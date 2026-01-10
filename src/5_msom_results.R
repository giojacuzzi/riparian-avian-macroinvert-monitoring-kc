# 5_msom_results.R ===================================================================
# Visualize multi-species occupancy modeling results and estimate site-specific richness
#
# Input
msom_path = "data/cache/models/msom_migrant.rds"
# Output
path_out = paste0("data/cache/5_msom_results/msom_richness_estimates.rds")

source("src/global.R")

# Naive occupancy -----------

# in_cache_detections = "data/cache/1_preprocess_agg_pam_data/detections_calibrated_0.5.rds"
# detections = readRDS(in_cache_detections)
# detections$long$common_name = tolower(detections$long$common_name)

# Load data for multi-species occupancy model --------------------------------------------------

message("Loading data for multi-species occupancy model ", msom_path)
model_data = readRDS(msom_path)

msom_summary = model_data$msom_summary
msom = model_data$msom
groups = model_data$groups %>% arrange(common_name)
sites = model_data$sites
species = model_data$species

#  Calculate estimated group richness per site ---------------------------
# z = msom$sims.list$z
# group_idx = groups$group_idx
# group_names = groups %>% distinct(group_idx, .keep_all = TRUE) %>% arrange(group_idx) %>% pull(group)
# samples = dim(z)[1]
# J = dim(z)[2]
# I = dim(z)[3]
# G = max(group_idx)
# 
# rich_group = array(NA, c(samples, J, G))
# for (g in 1:G) {
#   species_in_g = which(group_idx == g)
#   rich_group[ , , g] = apply(z[ , , species_in_g, drop = FALSE], c(1,2), sum)
# }
# rich_group_mean  = apply(rich_group, c(2,3), mean)
# rich_group_lower = apply(rich_group, c(2,3), quantile, probs = 0.025)
# rich_group_upper = apply(rich_group, c(2,3), quantile, probs = 0.975)
# msom_richness_estimates = lapply(1:G, function(g) {
#   group = group_names[g]
#   tibble(
#     site_id = sites
#   ) %>%
#     dplyr::mutate(
#       !!paste0(group, "_rich_mean")  := rich_group_mean[, g],
#       !!paste0(group, "_rich_lower") := rich_group_lower[, g],
#       !!paste0(group, "_rich_upper") := rich_group_upper[, g]
#     )
# })
# names(msom_richness_estimates) = group_names
# 
# rich_total = apply(z, c(1, 2), sum)
# msom_richness_estimates$total = tibble(
#   site_id    = sites,
#   total_rich_mean  = apply(rich_total, 2, mean),
#   total_rich_lower = apply(rich_total, 2, quantile, probs = 0.025),
#   total_rich_upper = apply(rich_total, 2, quantile, probs = 0.975)
# )

# if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
# saveRDS(msom_richness_estimates, file = path_out)
# message(crayon::green("Cached site-specific richness estimates to", path_out))

# Visualize results ---------------------

## Inspect group membership

ggplot(groups %>%
         left_join(species_traits, by = "common_name") %>%
         mutate(migration = as.factor(migration)) %>%
         count(group, migration) %>% group_by(group) %>%
         mutate(prop = n / sum(n)) %>% ungroup(),
       aes(x = group, y = prop, fill = migration)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_discrete(
    name = "Migration",
    breaks = c("1", "2", "3"),
    labels = c("Resident", "Partial migrant", "Long-distance migrant")
  ) +
  labs(x = "Group", y = "Proportion", fill = "Migration")

ggplot(groups %>%
         left_join(species_traits, by = "common_name") %>%
         mutate(primary_lifestyle = as.factor(primary_lifestyle)) %>%
         count(group, primary_lifestyle) %>% group_by(group) %>%
         mutate(prop = n / sum(n)) %>% ungroup(),
       aes(x = group, y = prop, fill = primary_lifestyle)) +
  geom_col(position = "fill") + 
  scale_y_continuous(labels = scales::percent_format())

## Retrieve baseline responses

community_baselines = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.u", "sigma.u", "mu.v", "sigma.v",
      "mu.w", "sigma.w", "mu.b", "sigma.b"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, prob, prob_lower95, prob_upper95)

community_baselines = community_baselines %>% left_join(groups %>% distinct(group_idx, group), by = "group_idx")

message("Baseline occurrence probability:")
print(subset(community_baselines, param == "mu.u"))

message("Baseline detection probability:")
print(subset(community_baselines, param == "mu.v"))

ggplot(community_baselines, aes(x = prob, y = param, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = prob_lower95, xmax = prob_upper95), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Parameter",
       title = "Hyperparameter estimates across models (mean and 95% BCI)",
       subtitle = "
u - occupancy probability across sites
v - true positive detection probability given presence
") + theme(legend.position = "bottom")

species_baselines = model_data$msom_summary %>%
  filter(grepl("^[uvwb]\\[", param)) %>%
  select(param, prob, prob_lower95, prob_upper95) %>%
  mutate(species_idx = str_extract(param, "(?<=\\[)[0-9]+") %>% as.integer()) %>%
  mutate(common_name = species[species_idx]) %>%
  left_join(groups, by = "common_name")

message("Species occurrence probability range:")
species_baselines %>% filter(startsWith(param, "u")) %>% summarise(
  prob_min = min(prob),
  species_min = common_name[which.min(prob)],
  prob_max = max(prob),
  species_max = common_name[which.max(prob)],
  .groups = "drop"
)

ggplot(species_baselines %>% filter(startsWith(param, "u")), aes(x = prob, y = reorder(common_name, prob), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("orange", "gray20")) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific occurrence probability `u` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

ggplot(species_baselines %>% filter(startsWith(param, "v")), aes(x = prob, y = reorder(common_name, prob), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("orange", "gray20")) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific TP detection probability given presence `v` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

# TODO
# Nsite_posterior = model_data$msom_summary %>% filter(stringr::str_starts(param, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
#   mutate(site_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(site = sites[site_idx])
# Nsite_mean = mean(Nsite_posterior$mean)
# message("Mean estimated species richness across all sites: ", round(Nsite_mean,1), " (range ", round(min(Nsite_posterior$mean),1), "–", round(max(Nsite_posterior$mean),1), ")")
# ggplot(Nsite_posterior, aes(x = as.factor(plot_order), y = mean)) +
#   geom_hline(yintercept = Nsite_mean, linetype = "solid", color = "blue") +
#   geom_point() + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
#   scale_x_discrete(labels = Nsite_posterior$site) + 
#   labs(title = "Estimated species richness per site", x = "Site", y = "Estimated species richness") +
#   coord_flip()
# ggplot(Nsite_posterior, aes(x = mean)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(xintercept = Nsite_mean, color = "blue") +
#   labs(title = "Estimated species richness per site", x = "Number of species estimated", y = "Number of sites")

# Compare coefficients on occupancy -------------------------------------------------------------------------

message("Comparing model coefficients")

# occ_coef_summary = model_data$msom_summary %>%
#   filter(str_detect(param, "^(mu|sigma)\\.alpha")) %>% arrange(param) %>%
#   select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

occ_coef_summary = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.alpha", "sigma.alpha"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
  left_join(groups %>% distinct(group_idx, group), by = "group_idx")

# detect_coef_summary = model_data$msom_summary %>%
#   filter(str_detect(param, "^(mu|sigma)\\.beta")) %>% arrange(param) %>%
#   select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

detect_coef_summary = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.beta", "sigma.beta"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
  left_join(groups %>% distinct(group_idx, group), by = "group_idx")

param_alpha_data = model_data$param_alpha_data
param_beta_data = model_data$param_beta_data
  
param_occ_data = param_alpha_data

occ_effect_sizes = full_join(occ_coef_summary %>% filter(str_starts(param, "mu")), param_occ_data %>% mutate(param = paste0("mu.", param)), by='param')

# occ_effect_sizes = occ_effect_sizes %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

# Compare species level effects of each covariate on occurrence
ggplot(occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
    labs(title = "Community level effect sizes for occurrence covariates",
         subtitle = tools::file_path_sans_ext(basename(msom_path)),
         x = "Coefficient estimate", y = "Parameter")

species_effects = param_occ_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])

species_effects = species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")

# species_effects = species_effects %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))
# species_effects = left_join(species_effects, groups, by = "common_name")

mu_list = list(
  mu.alpha1 = model_data$msom$sims.list[["mu.alpha1"]],
  mu.alpha2 = model_data$msom$sims.list[["mu.alpha2"]],
  mu.alpha3 = model_data$msom$sims.list[["mu.alpha3"]],
  mu.alpha4 = model_data$msom$sims.list[["mu.alpha4"]],
  mu.alpha5 = model_data$msom$sims.list[["mu.alpha5"]]
)

# Convert each to long format and combine
mu_long = bind_rows(
  lapply(names(mu_list), function(mu_name) {
    mat = mu_list[[mu_name]]
    df = as.data.frame(mat)
    df$samp = 1:nrow(df)
    df_long = pivot_longer(df, cols = -samp, names_to = "species_group", values_to = "value")
    df_long$parameter = mu_name
    return(df_long)
  })
)
# mu_long$species_group = factor(mu_long$species_group, labels = c("invert_predator", "NA"))

# Plot posterior densities
ggplot(mu_long, aes(x = value, fill = species_group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free", labeller = labeller(parameter = c(mu.alpha1 = "B-IBI", mu.alpha2 = "Canopy %", mu.alpha3 = "Impervious %"))) +
  # scale_fill_manual(values = c("dodgerblue2", "gray20")) +
  labs(x = "Posterior value", y = "Density", fill = "Group")
  
# ===========================================================================================================
ggplot(species_effects, aes(x = coef_mean, y = name, group = group, color = group)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0.5, color = "gray90", linetype = "dashed") +
  geom_vline(xintercept = -0.5, color = "gray90", linetype = "dashed") +
  geom_beeswarm(aes(shape = coef_overlap0), dodge.width = 0.5, cex = 1, priority = "density", size = 2, alpha = 0.6) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), xmin = `2.5%`, xmax = `97.5%`), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), xmin = `25%`, xmax = `75%`), size = 1, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), group = as.factor(group)), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 1)) +
  # scale_color_manual(values = c("orange", "gray20")) +
  labs(x = "Occupancy coefficient mean", y = "", color = "Diet group") +
  theme_sleek() + guides(shape = "none")

p_bibi = ggplot(species_effects %>% filter(name == "bibi"), aes(x = coef_mean, y = reorder(common_name, coef_mean), color = group)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  # scale_color_manual(values = c("orange", "gray20")) +
  geom_point(); print(p_bibi)

stop()

# Some species without documented diets containing target taxa nevertheless exhibited positive trends in occupancy in response to B-IBI
# EX: Spiders are preferred prey of Chestnut-backed chickadee, Black-headed Grosbeak, which may benefit from B-IBI
# These (black-throated gray) species have shown trophic links to aquatic resources: https://doi.org/10.1002/ecs2.3148
# EX: Belted Kingfisher, which feeds supplementally but is not an invertivore

# TODO: Compare species rarity between groups -- do inverts have proportionally more rare species?

# “The most widespread species (PSF, GCK, CY) show negative mean responses to B-IBI, while less common species tend to show weak or uncertain positive responses.”
# The more common species provide more information, and disproportionately influence the hyperparameter

binwidth = 0.04
species_effects$coef_bin = round(species_effects$coef_mean / binwidth) * binwidth

p_occ = ggplot(species_effects, aes(x = coef_bin, y = name, group = group)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`, color = as.factor(group)), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `25%`, xmax = `75%`, color = as.factor(group)), size = 1.0, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), group = as.factor(group)), size = 3, position = position_dodge(width = 0.5)) +
  geom_beeswarm(aes(color = group, alpha = coef_f),
    dodge.width = 0.5, shape = 16, cex = 1, priority = "density", size = 1.1, side=1L) +
  # geom_text_repel(data = species_effects, aes(x = coef_mean, y = name, label = common_name, color = group), size = 1, direction = "y", hjust = 0.05, max.overlaps = 20, position = position_dodge(0.5)) +
  # scale_color_manual(values = c("orange", "gray20")) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  labs(x = "Occupancy coefficient mean", y = "", color = "Diet group", alpha = "coef_f") +
  theme_sleek() + guides(shape = "none") + theme(
    legend.position = "bottom"
  ); print(p_occ)


param_detect_data = param_beta_data

detect_effect_sizes = full_join(detect_coef_summary %>% filter(str_starts(param, "mu")), param_detect_data %>% mutate(param = paste0("mu.", param)), by='param')
detect_effect_sizes = detect_effect_sizes %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

detect_species_effects = param_detect_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])
detect_species_effects = detect_species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")
detect_species_effects = detect_species_effects %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

detect_species_effects$coef_bin = round(detect_species_effects$coef_mean / binwidth) * binwidth
p_detect = ggplot(detect_species_effects, aes(x = coef_bin, y = name, group = group_label)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`, color = as.factor(group_label)), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `25%`, xmax = `75%`, color = as.factor(group_label)), size = 1.0, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), group = as.factor(group_label)), size = 3, position = position_dodge(width = 0.5)) +
  geom_beeswarm(aes(color = group_label, alpha = coef_f),
                shape = 16, dodge.width = 0.5, cex = 1, priority = "density", size = 1.1, side=1L) +
  # geom_text_repel(data = detect_species_effects, aes(x = coef_mean, y = name, label = common_name, color = group_label), size = 1, direction = "y", hjust = 0.05, max.overlaps = 20, position = position_dodge(0.5)) +
  scale_color_manual(values = c("orange", "gray20")) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  labs(x = "Detection coefficient mean", y = "Predictor", color = "Diet group", alpha = "coef_f") +
  theme_sleek() + guides(shape = "none") + theme(
    legend.position = "none"
  ); print(p_detect)

# Marginal responses --------------------------------------------------------------------

param_name = "alpha1"
param_data   = param_alpha_data %>% filter(param == param_name)

# Mean marginal probabilities of occurrence for the metacommunity in relation to alpha
coefs = model_data$msom_summary %>% filter(stringr::str_starts(param, param_name)) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(common_name = species[species_idx])
param_scaled_data = param_data %>% pull(scaled)
p_mu = attr(param_scaled_data[[1]], "scaled:center") # to transform between z-scale and pred_range_original scale
p_sd = attr(param_scaled_data[[1]], "scaled:scale")
bound_low  = min(param_scaled_data[[1]]) * p_sd + p_mu
bound_high = max(param_scaled_data[[1]]) * p_sd + p_mu
pred_range_original = seq(bound_low, bound_high, by = 0.01) # range of possible alpha values
pred_range_standardized = (pred_range_original - p_mu) / p_sd
# species-specific occurrence intercepts u[i]
intercepts = model_data$msom_summary %>% filter(str_starts(param, "u")) %>%
  mutate(species_idx = as.integer(str_extract(param, "\\d+")), common_name = species[species_idx]) %>%
  select(common_name, u_i = mean)
# species-specific coefficients
intercepts_and_coeffs = coefs %>% rename(alpha6_i = mean) %>% select(common_name, alpha6_i) %>%
  left_join(intercepts, by = "common_name")
# predict species-specific occurrence probabilities
predictions = intercepts_and_coeffs %>%
  rowwise() %>% do({
    i <- .
    tibble(
      common_name = i$common_name,
      idx = pred_range_original,
      psi = plogis(i$u_i + i$alpha6_i * pred_range_standardized)
    )
  }) %>% bind_rows()
predictions = predictions %>% left_join(groups %>% select(common_name, group), by = "common_name")

mu_u_samples      = as.matrix(model_data$msom$sims.list[["mu.u"]])
mu_param_samples  = as.matrix(model_data$msom$sims.list[[paste0("mu.", param_name)]])
n_groups = ncol(mu_u_samples)

meta_preds = map_dfr(1:nrow(mu_u_samples), function(i) {
  map_dfr(1:n_groups, function(g) {
    tibble(
      group_idx = g,
      idx = pred_range_original,
      psi = plogis(
        mu_u_samples[i, g] +
          mu_param_samples[i, g] * pred_range_standardized
      ),
      draw = i
    )
  })
})
meta_summary = meta_preds %>% # calculate means and 95% BCIs
  group_by(group_idx, idx) %>% summarize(
    psi_mean = mean(psi),
    psi_lower = quantile(psi, 0.025),
    psi_upper = quantile(psi, 0.975)
  )
meta_summary = meta_summary %>%
  left_join(
    groups %>% select(group_idx, group) %>% distinct(),
    by = "group_idx"
  )

ggplot() +
  # geom_line(data = predictions, aes(x = idx, y = psi, group = common_name, color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  # scale_x_continuous(limits = c(bound_low, bound_high)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = param_data$name, y = "Occurrence probability") +
  theme_sleek()

# ggplot(predictions %>% filter(common_name == "belted kingfisher"), aes(x = idx, y = psi)) +
#   geom_line() +
#   scale_x_continuous(limits = c(bound_low, bound_high)) +
#   scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0))
# 
# ggplot(predictions %>% filter(common_name == "wilson's warbler"), aes(x = idx, y = psi)) +
#   geom_line() +
#   scale_x_continuous(limits = c(bound_low, bound_high)) +
#   scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0))
# 
# ggplot(predictions %>% filter(common_name == "bewick's wren"), aes(x = idx, y = psi)) +
#   geom_line() +
#   scale_x_continuous(limits = c(bound_low, bound_high)) +
#   scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0))

# ggplot(predictions, aes(x = idx, y = psi, group = common_name)) +
#   geom_line(aes(color = group)) +
#   # gghighlight::gghighlight(use_direct_label = FALSE, unhighlighted_params = list(color = alpha("black", 0.1))) +
#   facet_wrap(~reorder(common_name, psi)) +
#   scale_x_continuous(limits = c(bound_low, bound_high)) +
#   scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0)) +
#   labs(title = "", x = param_data$name, y = "Occurrence probability") +
#   theme(legend.position = "none")

p_bibi = ggplot(predictions, aes(x = idx, y = psi, group = common_name)) +
  geom_line(aes(color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group), alpha = 0.2, inherit.aes = FALSE) +
  geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  facet_wrap(~group) +
  scale_fill_manual(values = c("orange", "gray20")) +
  scale_color_manual(values = c("orange", "gray20")) +
  scale_x_continuous(limits = c(bound_low, bound_high)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0)) + theme(
    legend.position = "none"
  ); print(p_bibi)

# Combo plots
p_detect + p_occ + p_bibi

