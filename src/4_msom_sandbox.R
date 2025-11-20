library(jagsUI)

in_cache_detections     = "data/cache/1_preprocess_agg_pam_data/detections_calibrated_0.5.rds"

# Load species detection history data
detections = readRDS(in_cache_detections)
detections$long$common_name = tolower(detections$long$common_name)
detections$wide$common_name = tolower(detections$wide$common_name)

site_data_reach = readRDS("data/cache/3_calculate_vars/NEW_site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/NEW_site_data_5000m.rds")

site_data = site_data_reach

path_out = paste0("data/cache/models/reach_all.rds")

model_file = "src/msom.txt"

# Format detection data as multidimensional array (site x survey x species)
detections_long = detections$long %>% rename(site = site_id, survey = survey_num, species = common_name)
# TODO: manually remove any species or surveys needed
detections_long$count = ifelse(detections_long$n_detections > 0, 1, 0)
sites   = unique(detections_long$site)
surveys = sort(unique(detections_long$survey))
species = unique(detections_long$species)
y = xtabs(count ~ site + survey + species, data = detections_long)
y = as.array(y)

# Model data constants
J = dim(y)[1]
K = dim(y)[2]
Kmax = max(K)
I = dim(y)[3]

# Model data covariates
site_data = site_data %>% slice(match(sites, site_id)) # match site order in y

# Store alpha parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
site_data = site_data %>% rename(
  canopy = rast_usfs_canopycover_sum_proportion,
  imp = rast_nlcd_impervious_sum_proportion
)
param_alpha_names = c("bibi", "canopy", "imp")
param_alpha_data = tibble(param = paste0("alpha", 1:length(param_alpha_names)), name  = param_alpha_names)
param_alpha_data = param_alpha_data %>% rowwise() %>% mutate(scaled = list(scale(site_data[[name]]))) %>% ungroup()
n_alpha_params = nrow(param_alpha_data)

# NOTE: For now, just use survey number as a detection covariate (replace with yday eventually)
x_yday = matrix(
  rep(as.vector(scale(surveys)), each = J),
  nrow = J,
  ncol = K,
  byrow = FALSE
)
detect_data = list(
  yday = x_yday
)
# Store beta parameter ID, variable name, and standardize data to have mean 0, standard deviation 1
param_beta_data = tibble(param = paste0("beta", seq_along(detect_data)), name = names(detect_data))
param_beta_data$scaled = list(detect_data[["yday"]])
n_beta_params = nrow(param_beta_data)

# Package data for MSOM -----------------------------------------------------------------------------------

msom_data = list(
  J = J,       # number of sites sampled
  K = rep(K, J),       # number of secondary sampling periods (surveys) per site per season (site x season) # NOTE: same number of surveys per site here!
  Kmax = Kmax, # maximum number of surveys across all sites and seasons
  I = I,       # number of species observed
  y = y        # observed (detection-nondetection) data matrix
)

# Add covariates to msom_data
for (a in seq_len(n_alpha_params)) { # Add alpha covariates
  msom_data[[paste0("x_", param_alpha_data$param[a])]] <- as.vector(param_alpha_data$scaled[[a]])
}
for (b in seq_len(n_beta_params)) {
  vec <- as.vector(param_beta_data$scaled[[b]])
  msom_data[[paste0("x_", param_beta_data$param[b])]] <- array(
    vec,
    dim = c(J,K)
  )
}

str(msom_data)

# Initialize latent occupancy state z[i] as 1 if a detection occurred at site i, and 0 otherwise
z = array(NA, dim = c(J, I), dimnames = list(sites, species))
for (j in seq_along(sites)) {
  for (i in seq_along(species)) {
    z[j, i] = (sum(y[j, , i], na.rm = TRUE) > 0) * 1
  }
}

# Run JAGS -----------------------------------------------------------------------------------------------

message("\n", "System CPU: "); print(as.data.frame(t(benchmarkme::get_cpu())))
message("System RAM: "); print(benchmarkme::get_ram())

message("Running JAGS (current time ", time_start <- Sys.time(), ")")

msom = jags(data = msom_data,
            inits = function() { list( # initial values to avoid data/model conflicts
              z = z
              # v = rep(logit(0.70), length(species)), # informative priors are necessary to avoid invalid PPC log(0) values
              # w = rep(logit(0.05), length(species))
            ) },
            parameters.to.save = c( # monitored parameters
              "mu.u", "sigma.u", "u",
              "mu.v", "sigma.v", "v",
              # "mu.w", "sigma.w", "w",
              # "mu.b", "sigma.b", "b",
              paste0("mu.alpha", 1:n_alpha_params), paste0("sigma.alpha", 1:n_alpha_params), paste0("alpha", 1:n_alpha_params),
              paste0("mu.beta",  1:n_beta_params),  paste0("sigma.beta",  1:n_beta_params),  paste0("beta",  1:n_beta_params),
              "D_obs", "D_sim"
            ),
            model.file = model_file,
            n.chains = 3, n.adapt = 100, n.iter = 1000, n.burnin = 100, n.thin = 1, parallel = FALSE, # ETA: 
            # n.chains = 3, n.adapt = 1000, n.iter = 10000, n.burnin = 2000, n.thin = 1, parallel = TRUE, # TODO: ETA
            DIC = FALSE, verbose=TRUE)

message("Finished running JAGS (", round(msom$mcmc.info$elapsed.mins / 60, 2), " hr)")

message("MCMC information:")
print(data.frame(
  n.chains = msom$mcmc.info$n.chains,
  n.adapt = msom$mcmc.info$n.adapt[1],
  n.iter = msom$mcmc.info$n.iter,
  n.burnin = msom$mcmc.info$n.burnin,
  n.thin = msom$mcmc.info$n.thin,
  samples = msom$mcmc.info$n.samples
))

message("Model size: ", format(utils::object.size(msom), units = "MB"))

# Retrieve summary data and investigate goodness-of-fit ----------------------------------------------------

msom_summary = summary(msom)
msom_summary = msom_summary %>% as_tibble() %>%
  mutate(param = rownames(summary(msom)), overlap0 = as.factor(overlap0)) %>% relocate(param, .before = 1) %>%
  mutate(prob = plogis(mean), prob_lower95 = plogis(`2.5%`), prob_upper95 = plogis(`97.5%`))
rhat_threshold = 1.1
suspected_nonconvergence = msom_summary %>% filter(Rhat >= rhat_threshold) %>% filter(!str_starts(param, "z\\[") & !str_starts(param, "psi\\["))
suspected_nonconvergence = suspected_nonconvergence %>% mutate(
  index = str_extract(param, "(?<=\\[)\\d+(?=\\])"),
  index = as.integer(index),
  species = ifelse(!is.na(index), species[index], NA)
)
if (nrow(suspected_nonconvergence) > 1) {
  message("The following ", nrow(suspected_nonconvergence), " parameters may not have converged:")
  print(suspected_nonconvergence)
} else {
  message("All parameters appear to have converged (rhat < ", rhat_threshold, ")")
}

## Posterior predictive check - Bernoulli deviance contribution (Broms et al. 2016)
# "If the observed data are consistent with the model in question, then the Bayesian p-value should be close to 0.5. In practice, a p-value close to 0 or 1 indicates that the model is inadequate in some way -- close to 0 suggests a lack of fit and close to 1 suggests that the model over-fits the data, which may occur when it is too complex... A Bayesian p-value is calculated as the proportion of times the selected summary statistic calculated for the generated data is greater than the value calculated from the observed data" (MacKenzie et al. 2018)
# Is the overall likelihood of the observed detection histories under the fitted model about the same as the likelihood of new data generated from that model?
# This Bayesian p-value is the probability (proportion of iterations) that the simulated deviance is greater than the observed deviance.
p_val = mean(msom$sims.list$D_sim > msom$sims.list$D_obs)
message("Baysian p-value (deviance): ", round(p_val,3))

if (FALSE) {
  # "Examine trace plots for good mixing and convergence among chains. Each chain is displayed in a different colour. This means random paths exploring a lot of the parameter space on the y-axis without a clear pattern and each chain converging on the same value."
  MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'trace', Rhat = TRUE, n.eff = TRUE)
  
  # "Examine density plots for not super-wide or with irregular peaks. The more parameter space the density plots include, the higher the uncertainty in a parameter estimate. The density curves don’t have to be normal but shouldn’t have multiple peaks and each chain colour should have approximately the same peak."
  MCMCtrace(msom$samples, ISB = FALSE, pdf = F, exact = TRUE, post_zm = TRUE, type = 'density', Rhat = TRUE, n.eff = TRUE, ind = TRUE)
}

# Inspect the mean and 95% BCI of hyperparameter estimates
whiskerplot(msom, c(paste0('mu.', param_alpha_data$param)))
whiskerplot(msom, c(paste0('mu.', param_beta_data$param)))

# Write results to cache
msom_results = list(
  msom              = msom,
  msom_summary      = msom_summary,
  p_val             = p_val,
  param_alpha_data  = param_alpha_data,
  param_beta_data   = param_beta_data,
  sites             = sites,
  species           = species
)
if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
saveRDS(msom_results, file = path_out)
message(crayon::green("Cached model and results to", path_out))
