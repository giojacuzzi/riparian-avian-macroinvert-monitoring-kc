# 4_sem.R ===================================================================
# Fit piecewise structural equation models
#
# Input
use_msom_richness_estimates = TRUE # use richness estimates from the msom instead of naive observed values
msom_path = "data/cache/models/msom_all.rds"
in_cache_detections = "data/cache/1_preprocess_agg_pam_data/detections_calibrated_0.5.rds" # detections_calibrated_0.75.rds
# Output
out_cache_dir = "data/cache/4_sem"

source("src/global.R")

# Load site variable data ------------------------------------------------------------

site_data_reach = readRDS("data/cache/3_calculate_vars/site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/site_data_5000m.rds")

site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)

# Urbanization x B-IBI gradient
ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = bibi)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.01) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "orange", alpha = 0.01) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "yellow", alpha = 0.01) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "green", alpha = 0.01) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "dodgerblue", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Canopy x B-IBI gradient
ggplot(site_data_reach, aes(x = rast_usfs_canopycover_sum_proportion, y = bibi)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.01) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "orange", alpha = 0.01) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "yellow", alpha = 0.01) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "green", alpha = 0.01) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "dodgerblue", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Urbanization x Canopy gradient
ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = rast_usfs_canopycover_sum_proportion)) +
  geom_rect(aes(xmin = 0.0, xmax = 0.20, ymin = -Inf, ymax = Inf), fill = "pink", alpha = 0.01) +
  geom_rect(aes(xmin = 0.20, xmax = 0.50, ymin = -Inf, ymax = Inf), fill = "tomato", alpha = 0.01) +
  geom_rect(aes(xmin = 0.50, xmax = 0.80, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Load species detection history data ------------------------------------------------
message("Loading species detection history data")

detections = readRDS(in_cache_detections)
detections$long$common_name = tolower(detections$long$common_name)
detections$wide$common_name = tolower(detections$wide$common_name)

# Exclude certain sites from analysis ------------------------------------------------

# Inspect total detections per site
# total_detections_by_site = detections$long %>% group_by(site_id) %>%
#   summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
#   arrange(desc(total_detections))
# ggplot(total_detections_by_site, aes(x = reorder(site_id, total_detections), y = total_detections)) + geom_col()

message("Excluding site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

message("Retaining ", length(unique(site_data_reach$site_id)), " sites")

# Visualize joint data ----------------------------------------------------------------
message("Visualizing joint data")

presence_absence = detections$long %>% group_by(site_id, common_name) %>%
  summarise(presence = if_else(sum(n_detections, na.rm = TRUE) > 0, 1, 0), .groups = "drop")

# Summarize richness of different groups by site
site_group_richness = presence_absence %>%
  group_by(site_id) %>%
  summarise(
    obs_all = sum(presence) # observed richness of all species
  )

# Incorporate mean msom richness estimates
if (use_msom_richness_estimates) {
  
  message("Loading data for multi-species occupancy model ", msom_path)
  msom_data = readRDS(msom_path)
  
  msom_summary = msom_data$msom_summary
  msom    = msom_data$msom
  groups  = msom_data$groups %>% arrange(common_name)
  sites   = msom_data$sites
  species = msom_data$species
  
  z = msom$sims.list$z
  # group_idx = groups$group_idx
  # group_names = groups %>% distinct(group_idx, .keep_all = TRUE) %>% arrange(group_idx) %>% pull(group)
  samples = dim(z)[1]
  J = dim(z)[2]
  I = dim(z)[3]
  # G = max(group_idx)
  
  # species_groupings = list(
  #   species %in% (species_traits %>% filter(group_all == "all")         %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_migrant == "migrant") %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_diet == "diet")       %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_forage == "aerial")   %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_forage == "gleaner")  %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_forage == "ground")   %>% pull(common_name)),
  #   species %in% (species_traits %>% filter(group_forage == "bark")     %>% pull(common_name))
  # )
  
  message("Calculating richness for each grouping")
  pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = 4, clear = FALSE)
  for (grouping in c("group_all", "group_migrant", "group_diet", "group_forage")) {
  
    species_group = tibble(
      common_name = species_traits$common_name,
      group = species_traits[[grouping]],
      group_idx = as.integer(as.factor(group))
    )
    
    G = length(unique(species_group$group))
    
    rich_group = array(NA, c(samples, J, G))
    for (g in unique(species_group$group)) {
      species_in_g = which(species %in% (species_group %>% filter(group == g) %>% pull(common_name)))
      g = (species_group %>% filter(group == g) %>% pull(group_idx))[1]
      rich_group[ , , g] = apply(z[ , , species_in_g, drop = FALSE], c(1,2), sum)
    }
    rich_group_mean  = apply(rich_group, c(2,3), mean)
    rich_group_lower = apply(rich_group, c(2,3), quantile, probs = 0.025)
    rich_group_upper = apply(rich_group, c(2,3), quantile, probs = 0.975)
    msom_richness_estimates = lapply(1:G, function(g) {
      g_name = (species_group %>% filter(group_idx == g) %>% pull(group))[1]
      gn = paste0(grouping, "_", g_name)
      tibble(
        site_id = sites
      ) %>%
        dplyr::mutate(
          !!paste0(gn, "_rich_mean")  := rich_group_mean[, g],
          !!paste0(gn, "_rich_lower") := rich_group_lower[, g],
          !!paste0(gn, "_rich_upper") := rich_group_upper[, g]
        )
    })
    msom_richness_estimates = Reduce(
      function(x, y) left_join(x, y, by = "site_id"),
      msom_richness_estimates
    )
    
    site_group_richness = left_join(site_group_richness, msom_richness_estimates, by = "site_id")
    pb$tick()
  }
}

# ## Indicator species analysis
# library(indicspecies)
# spmat = presence_absence %>% pivot_wider(names_from = common_name, values_from = presence)
# d = left_join(spmat, site_data_reach, by = "site_id")
# spmat = spmat %>% select(-site_id) %>% as.data.frame()
# bibi_excellent = ifelse(d$bibi >= 80, 1, 0)
# indval = multipatt(spmat, bibi_excellent, control = how(nperm=999)) 
# summary(indval)

# Arrange richness to match site data
stopifnot(site_data_reach$site_id == site_data_basin$site_id)
site_group_richness = site_group_richness %>% arrange(match(site_id, site_data_reach$site_id))

# Species-specific presence/absence as a function of BIBI
ggplot(left_join(presence_absence %>% filter(common_name == "wilson's warbler"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

# Structural equation modeling -----------------------------------------------------------------------------

# Calculate pairwise collinearity among predictors

# Urbanization variables
pairwise_collinearity(site_data_reach %>% st_drop_geometry() %>% select(
  rast_nlcd_impervious_sum_proportion, density_roads_paved))
pairwise_collinearity(site_data_basin %>% st_drop_geometry() %>% select(
  rast_nlcd_impervious_sum_proportion, density_roads_paved))

# Cover configuration variables
pairwise_collinearity(site_data_reach %>% st_drop_geometry() %>% select(
  pd_riphab, ed_riphab))
pairwise_collinearity(site_data_basin %>% st_drop_geometry() %>% select(
  pd_riphab, ed_riphab))

# Riparian vegetation variables
pairwise_collinearity(site_data_reach %>% st_drop_geometry() %>% select(
  rast_usfs_canopycover_sum_proportion, nlcd_forest, nlcd_forest_and_wetlands, rast_gedi_height_mean, rast_gedi_fhd_mean))
pairwise_collinearity(site_data_basin %>% st_drop_geometry() %>% select(
  rast_usfs_canopycover_sum_proportion, nlcd_forest, nlcd_forest_and_wetlands, rast_gedi_height_mean, rast_gedi_fhd_mean))

# Combinations
pairwise_collinearity(site_data_reach %>% st_drop_geometry() %>% select(
  rast_nlcd_impervious_sum_proportion, pd_riphab, rast_usfs_canopycover_sum_proportion))
pairwise_collinearity(site_data_basin %>% st_drop_geometry() %>% select(
  rast_nlcd_impervious_sum_proportion, pd_riphab, rast_usfs_canopycover_sum_proportion))

# 550 m represents riparian zone within the local reach, and the 90% dispersal distance
# 5 km represents the catchment landscape (roughly basin)
d_raw = data.frame(
  rich_mean_all            = site_group_richness$group_all_all_rich_mean,
  rich_mean_migrant        = site_group_richness$group_migrant_migrant_rich_mean,
  rich_mean_diet           = site_group_richness$group_diet_diet_rich_mean,
  rich_mean_forage_aerial  = site_group_richness$group_forage_aerial_rich_mean,
  rich_mean_forage_gleaner = site_group_richness$group_forage_gleaner_rich_mean,
  rich_mean_forage_ground  = site_group_richness$group_forage_ground_rich_mean,
  rich_mean_forage_bark    = site_group_richness$group_forage_bark_rich_mean,
  # B-IBI
  "bibi"            = site_data_reach$bibi,
  # Environmental variables
  "imp_reach"       = site_data_reach$rast_nlcd_impervious_sum_proportion,
  "imp_basin"       = site_data_basin$rast_nlcd_impervious_sum_proportion,

  "ed_reach"        = site_data_reach$ed,
  "ed_basin"        = site_data_basin$ed,
  "ed_riphab_reach" = site_data_reach$ed_riphab,
  "ed_riphab_basin" = site_data_basin$ed_riphab,
  "ed_dev_reach"    = site_data_reach$ed_dev,
  "ed_dev_basin"    = site_data_basin$ed_dev,
  "agg_reach"       = site_data_reach$agg,
  "agg_basin"       = site_data_basin$agg,
  "pd_reach"        = site_data_reach$pd,
  "pd_basin"        = site_data_basin$pd,
  "pd_riphab_reach" = site_data_reach$pd_riphab,
  "pd_riphab_basin" = site_data_basin$pd_riphab,
  "pd_dev_reach"    = site_data_reach$pd_dev,
  "pd_dev_basin"    = site_data_basin$pd_dev,
  
  "roads_reach"      = site_data_reach$density_roads_paved,
  "roads_basin"      = site_data_basin$density_roads_paved,
  
  "dev_reach"     = site_data_reach$nlcd_developed_variable_intensity + site_data_reach$nlcd_developed_open_space,
  "dev_basin"     = site_data_basin$nlcd_developed_variable_intensity + site_data_basin$nlcd_developed_open_space,
  "devvi_reach"   = site_data_reach$nlcd_developed_variable_intensity,
  "devvi_basin"   = site_data_basin$nlcd_developed_variable_intensity,
  "forest_reach"  = site_data_reach$nlcd_forest,
  "forest_basin"  = site_data_basin$nlcd_forest,
  "riphab_reach"  = site_data_reach$nlcd_forest_and_wetlands,
  "riphab_basin"  = site_data_basin$nlcd_forest_and_wetlands,
  
  "tcc_reach"  = site_data_reach$rast_usfs_canopycover_sum_proportion,
  "tcc_basin"  = site_data_basin$rast_usfs_canopycover_sum_proportion,
  
  "height_reach"  = site_data_reach$rast_gedi_height_mean,
  "height_basin"  = site_data_basin$rast_gedi_height_mean,
  
  "fhd_reach"     = site_data_reach$rast_gedi_fhd_mean,
  "fhd_basin"     = site_data_basin$rast_gedi_fhd_mean,
  "site_id"       = site_data_reach$site_id
)

# Print summary stats for environmental variables
num_cols = names(d_raw)[sapply(d_raw, is.numeric)]
for (col in num_cols) {
  cat("\n", col, "\n")
  cat("  Mean:", mean(d_raw[[col]], na.rm = TRUE), "\n")
  cat("  SD:  ", sd(d_raw[[col]], na.rm = TRUE), "\n")
  cat("  Min: ", min(d_raw[[col]], na.rm = TRUE), "\n")
  cat("  Max: ", max(d_raw[[col]], na.rm = TRUE), "\n")
}
pairwise_collinearity(d_raw %>% select(where(is.numeric)))

## SEM =====================================================================
d = d_raw
n = nrow(d)

# Canopy reach (logit transformation for proportional bounds)
d$tcc_reach = qlogis(d_raw$tcc_reach)

# Canopy basin (logit transformation for proportional bounds)
d$tcc_basin = qlogis(d_raw$tcc_basin)

# B-IBI (logit transformation for proportional bounds)
d$bibi = d_raw$bibi / 100
d$bibi = (d$bibi * (n - 1) + 0.5) / n # Smithson & Verkuilen (2006) squeeze
d$bibi = qlogis(d$bibi)

# NOTE: piecewiseSEM does not support glm Gamma distribution; must use normal distribution

## Model selection
options(na.action = "na.fail")

# Confirm the spatial scale by comparing candidate sets of predictors within each component model

# Model selection for bibi; optimal: bibi ~ imp_basin + tcc_reach OR imp_reach
# "Model selection indicated strong support for basin-scale impervious cover as a predictor of aquatic biotic integrity. Several competitive models (ΔAICc ≤ 2) included an additional local-scale predictor, but support was divided between reach-scale imperviousness and reach-scale canopy cover. We retained reach-scale canopy cover based on a priori hypotheses regarding local habitat mediation and to avoid redundant representation of urbanization at multiple scales."
m_bibi_global = lm(bibi ~ tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_bibi = dredge(m_bibi_global)

# Model selection for diet; optimal: diet ~ bibi + tcc_reach OR imp_reach
# “For the diet guild, AICc-based model selection indicated strong support for aquatic biotic integrity (BIBI) as a predictor of species richness, with all top-ranked models (ΔAICc ≤ 1) including BIBI. Support for additional local predictors was divided between reach-scale canopy cover and reach-scale imperviousness, indicating that these variables represent alternative expressions of the same local urbanization gradient. Because reach-scale canopy is largely determined by impervious cover, we retained imperviousness as the direct predictor and treated canopy as a mediating variable elsewhere in the SEM.”
m_diet_global   = lm(rich_mean_diet ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_diet = dredge(m_diet_global)

# Model selection for migrant; optimal: migrant ~ bibi + tcc_reach OR imp_reach
# “For neotropical migrant richness, AICc-based model selection also indicated strong and consistent support for aquatic biotic integrity (BIBI), with all top-ranked models (ΔAICc ≤ 1) including BIBI. Support for additional local predictors was divided between reach-scale canopy cover and reach-scale imperviousness, indicating that these variables represent alternative expressions of the same local urbanization gradient. We therefore retained imperviousness as the direct predictor and treated canopy as a mediating variable elsewhere in the SEM.”
m_migrant_global   = lm(rich_mean_migrant ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_migrant = dredge(m_migrant_global)

# Model selection for forage aerial
# “For aerial foragers, AICc-based model selection indicated substantial model uncertainty, with no one terrestrial or aquatic predictor clearly dominating at either spatial scale. The top-ranked model included BIBI alone, and additional predictors provided only marginal improvements (ΔAICc > 1). We therefore retained a single, theory-consistent path from aquatic condition to aerial forager richness, reflecting a weak but plausible trophic linkage.”
m_aerial_global   = lm(rich_mean_forage_aerial ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_aerial = dredge(m_aerial_global)

# Model selection for forage gleaner
# “For foliage gleaners, terrestrial predictors at multiple spatial scales (local imperviousness and basin-scale canopy) were among the top-supported variables, whereas aquatic condition (BIBI) received little support.”
m_gleaner_global   = lm(rich_mean_forage_gleaner ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_gleaner = dredge(m_gleaner_global)

# Model selection for forage ground
# “For ground-foraging birds, AICc-based model selection indicated that aquatic condition (BIBI) was the strongest predictor of richness, whereas terrestrial variables at either spatial scale had weak and inconsistent support across competitive models.”
m_ground_global   = lm(rich_mean_forage_ground ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_ground = dredge(m_ground_global)

# Model selection for forage bark
# “For bark-foraging birds, AICc-based model selection indicated substantial model uncertainty, with terrestrial habitat predictors at multiple spatial scales appearing among the top models, whereas aquatic condition (BIBI) had little support.”
m_bark_global   = lm(rich_mean_forage_bark ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_bark = dredge(m_bark_global)

# Model selection for all
# “For total bird richness, AICc-based model selection indicated substantial model uncertainty, with no single terrestrial or aquatic predictor clearly dominating at either spatial scale, and the top model containing no predictors besides the intercept. This suggests that overall richness may be less sensitive to landscape or aquatic variation than guild-specific richness patterns.”
m_all_global   = lm(rich_mean_all ~ bibi + tcc_reach + imp_reach + tcc_basin + imp_basin, d)
candidates_all = dredge(m_all_global)

## SEMs
coefs_sems = tibble()

# Diet group richness
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_diet      = lm(rich_mean_diet ~ bibi + tcc_reach, d)
sem_diet    = psem(m_tcc_reach, m_bibi, m_diet); plot(sem_diet); print(summary(sem_diet))
coefs_sems  = rbind(coefs_sems, coefs(sem_diet) %>% clean_names() %>% mutate(model = "sem_diet"))

# Migrant group richness
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_migrant   = lm(rich_mean_migrant ~ bibi + tcc_reach, d)
sem_migrant = psem(m_tcc_reach, m_bibi, m_migrant); plot(sem_migrant); print(summary(sem_migrant))
coefs_sems  = rbind(coefs_sems, coefs(sem_migrant) %>% clean_names() %>% mutate(model = "sem_migrant"))

# Aerial
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_aerial    = lm(rich_mean_forage_aerial ~ bibi + tcc_reach, d)
sem_aerial  = psem(m_tcc_reach, m_bibi, m_aerial); plot(sem_aerial); print(summary(sem_aerial))
coefs_sems  = rbind(coefs_sems, coefs(sem_aerial) %>% clean_names() %>% mutate(model = "sem_aerial"))

# Gleaner
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_gleaner   = lm(rich_mean_forage_gleaner ~ bibi + tcc_reach, d)
sem_gleaner = psem(m_tcc_reach, m_bibi, m_gleaner); plot(sem_gleaner); print(summary(sem_gleaner))
coefs_sems  = rbind(coefs_sems, coefs(sem_gleaner) %>% clean_names() %>% mutate(model = "sem_gleaner"))

# Ground
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_ground    = lm(rich_mean_forage_ground ~ bibi + tcc_reach, d)
sem_ground  = psem(m_tcc_reach, m_bibi, m_ground); plot(sem_ground); print(summary(sem_ground))
coefs_sems  = rbind(coefs_sems, coefs(sem_ground) %>% clean_names() %>% mutate(model = "sem_ground"))

# Bark
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_bark      = lm(rich_mean_forage_bark ~ bibi + tcc_reach, d)
sem_bark    = psem(m_tcc_reach, m_bibi, m_bark); plot(sem_bark); print(summary(sem_bark))
coefs_sems  = rbind(coefs_sems, coefs(sem_bark) %>% clean_names() %>% mutate(model = "sem_bark"))

# All (no group) richness
m_tcc_reach = lm(tcc_reach ~ imp_reach, d)
m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d)
m_all       = lm(rich_mean_all ~ bibi + tcc_reach, d)
sem_all     = psem(m_tcc_reach, m_bibi, m_all); plot(sem_all); print(summary(sem_all))
coefs_sems  = rbind(coefs_sems, coefs(sem_all) %>% clean_names() %>% mutate(model = "sem_all"))

# Save sem coefs to file
sem_coefs = coefs_sems %>%
  group_by(
    response, predictor, estimate, std_error, df,
    crit_value, p_value, std_estimate, x
  ) %>%
  summarise(model = paste(sort(unique(model)), collapse = ", "), .groups = "drop") %>% arrange(predictor) %>% rename(signif = x)
if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
out_filepath = file.path(out_cache_dir, paste0("sem_coefs.csv"))
write_csv(sem_coefs, out_filepath)
message(crayon::green("Cached", out_filepath))

# TODO: Test to support assumption of normality
# shapiro.test(residuals(m_predator)) # p > 0.05 => approximately normally distributed residuals

# Save SEMs for prediction
out_filepath = file.path(out_cache_dir, paste0("sem_migrant.rds"))
saveRDS(sem_migrant, out_filepath)
message(crayon::green("Cached", out_filepath))
out_filepath = file.path(out_cache_dir, paste0("sem_diet.rds"))
saveRDS(sem_diet, out_filepath)
message(crayon::green("Cached", out_filepath))

# Propagate uncertainty from MSOM --------------------------------------------------------
if (use_msom_richness_estimates) {
  
  coefs_final_stats = tibble()
  for (grouping in c("group_all", "group_migrant", "group_diet", "group_forage")) {
    
    species_group = tibble(
      common_name = species_traits$common_name,
      group = species_traits[[grouping]],
      group_idx = as.integer(as.factor(group))
    )
    
    G = length(unique(species_group$group))
    
    # Fit sem over all posterior draws to propagate uncertainty
    for (g in 1:(G+1)) {
      g_name = (species_group %>% filter(group_idx == g) %>% pull(group))[1]
      message("Propagating MSOM uncertainty for grouping '", grouping, "' group '", g_name, "' group richness response")
      z = msom$sims.list$z
      draws = dim(z)[1]
      stopifnot(all(groups$common_name == msom_data$species))
      coeffs_draws = tibble()
      pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = draws, clear = FALSE)
      for (draw in 1:draws) {
        # Calculate estimated group richness at each site
        rich_group_draw = data.frame(
          site_id = msom_data$sites,
          rich_group_draw = rich_group[draw, , g]
        )
        d_msom = d %>% left_join(rich_group_draw, by = "site_id")
        
        # Fit the component models and SEM
        m_tcc_reach = lm(tcc_reach ~ imp_reach, d_msom)
        m_bibi      = lm(bibi ~ tcc_reach + imp_basin, d_msom)
        m_rich_draw = lm(rich_group_draw ~ bibi + tcc_reach, d_msom)
        sem_draw    = psem(m_tcc_reach, m_bibi, m_rich_draw) #; plot(sem_draw); print(summary(sem_draw))
        
        # Extract SEM coefficients and store
        coeffs_draw = coefs(sem_draw) %>% clean_names() %>%
          # filter(response == "rich_group_draw") %>%
          mutate(draw = draw) %>% rename(signif = x)
        coeffs_draws = rbind(coeffs_draws, coeffs_draw)
        pb$tick()
      }
      
      rich_group_coefs = coeffs_draws %>%
        group_by(response, predictor) %>%
        summarise(
          mean  = mean(std_estimate),
          lower = quantile(std_estimate, 0.025),
          upper = quantile(std_estimate, 0.975),
          .groups = "drop"
        ) %>% rename(std_estimate = mean)
      
      rich_group_coefs$grouping = grouping
      rich_group_coefs$g_name = g_name

      p = ggplot(rich_group_coefs, aes(x = std_estimate, y = interaction(predictor, response))) +
        geom_vline(xintercept = 0, color = "gray") +
        geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
        geom_point() + labs(title = paste0("Grouping '", grouping, "' group '", g_name, "'")); print(p)
      
      message("SEM coefficients including propagated uncertainty:")
      print(rich_group_coefs)
      
      coefs_final_stats = rbind(coefs_final_stats, rich_group_coefs)
    }
  }
}

# TODO: print and save coefs_final_stats
stop("DEBUG")

# Check over/underdispersion -- if overdispersed, fit negative binomial
# simres_pois = simulateResiduals(m_predator, n = 1000); plot(simres_pois); testDispersion(simres_pois)

# TODO: Moran's I test for spatial autocorrelation
# library(spdep)
# # site coords
# coords <- cbind(stream_data$longitude, stream_data$latitude)
# # spatial weights matrix (inverse distance)
# dists <- as.matrix(dist(coords))
# inv_dists <- 1/dists
# diag(inv_dists) <- 0  # no self-weight
# listw <- mat2listw(inv_dists)
# # Moran's I for residuals
# resids <- residuals(your_model, type="response")  # e.g., BIBI residuals
# moran.test(resids, listw)

# print(presence_absence %>% group_by(common_name) %>% summarize(n_sites = sum(presence)) %>% arrange(n_sites) %>% filter(common_name %in% sp_invert), n = 100)
# 
# species_traits %>% filter(common_name %in% c(
#   "green-winged teal", "mallard", "gadwall", "wood duck", "common merganser", "blue-winged teal", "cackling goose", "caspian tern", "great blue heron", "red-winged blackbird", "american bittern", "green heron"
# )) %>% select(common_name, trophic_level, trophic_niche, diet_5cat)

# Marginal effect plots =======================================================================
{
  ## Canopy reach as a function of impervious reach
  # Generate prediction values
  imp_grid = seq(min(d$imp_reach), max(d$imp_reach), length.out = 300)
  d_pred_canopy = data.frame(imp_reach = imp_grid)
  pred_canopy = predict(m_canopy, newdata = d_pred_canopy, type = "response", se.fit = TRUE)
  d_pred_canopy = d_pred_canopy %>%
    mutate(fit  = pred_canopy$fit, se = pred_canopy$se.fit, low = fit - 1.96 * se, high = fit + 1.96 * se)
  
  # Back-transform predicted canopy to original proportional scale
  d_pred_canopy = d_pred_canopy %>%
    mutate(
      canopy      = plogis(fit)  * 100,
      canopy_low  = plogis(low)  * 100,
      canopy_high = plogis(high) * 100
    )
  # Back-transform observed canopy
  d$canopy_pct = plogis(d$canopy_reach) * 100
  
  # Plot
  p_canopy_imp = ggplot(d_pred_canopy %>% mutate(imp_reach), aes(x = imp_reach, y = canopy)) +
    geom_ribbon(aes(ymin = canopy_low, ymax = canopy_high), fill = "forestgreen", alpha = 0.25) +
    geom_line(color = "forestgreen", size = 1) +
    geom_point(data = d, aes(x = imp_reach, y = canopy_pct), color = "forestgreen", alpha = 0.4) +
    labs(x = "Reach impervious (%)", y = "Reach canopy cover (%)") +
    theme_sleek() + theme(aspect.ratio = 1)
  
  ## B-IBI as a function of impervious watershed
  # Generate prediction values, fixing canopy to mean value
  imp_grid = seq(min(d$imp_basin), max(d$imp_basin), length.out = 300)
  d_pred_bibi = data.frame(imp_basin = imp_grid, canopy_reach = mean(d$canopy_reach))
  pred_bibi = predict(m_bibi, newdata = d_pred_bibi, type = "response", se.fit = TRUE)
  d_pred_bibi = d_pred_bibi %>%
    mutate(fit = pred_bibi$fit, se = pred_bibi$se.fit, low = fit - 1.96 * se, high = fit + 1.96 * se)
  
  # Inverse Smithson & Verkuilen squeeze
  n = nrow(d)
  inverse_squeeze = function(p_sq) { ((p_sq * n) - 0.5) / (n - 1) }
  
  # Back-transform predicted B-IBI to original proportional scale
  d_pred_bibi = d_pred_bibi %>%
    mutate(
      bibi      = inverse_squeeze(plogis(fit))  * 100,
      bibi_low  = inverse_squeeze(plogis(low))  * 100,
      bibi_high = inverse_squeeze(plogis(high)) * 100
    )
  # Back-transform observed B-IBI
  d$bibi_pct_orig = inverse_squeeze(plogis(d$bibi)) * 100
  
  # Plot
  p_bibi_imp = ggplot(d_pred_bibi, aes(x = imp_basin, y = bibi)) +
    geom_ribbon(aes(ymin = bibi_low, ymax = bibi_high), fill = "royalblue", alpha = 0.25) +
    geom_line(color = "royalblue", size = 1) +
    geom_point(data = d, aes(x = imp_basin, y = bibi_pct_orig), color = "royalblue", alpha = 0.4) +
    labs(x = "Watershed impervious (%)", y = "B-IBI") +
    theme_sleek() + theme(aspect.ratio = 1)
  
  ## Predator richness as a function of B-IBI
  # Generate prediction values, fixing canopy and impervious at
  # mean values (start with original B-IBI proportional scale,
  # then transform to logit scale for the model)
  bibi_grid = seq(1, 99, length.out = 300)
  bibi_logit = qlogis(((bibi_grid / 100) * (n - 1) + 0.5) / n) # apply squeeze
  d_pred_rich = data.frame(bibi = bibi_logit, canopy_reach = mean(d$canopy_reach), imp_reach = mean(d$imp_reach, na.rm = TRUE))
  if (use_msom_richness_estimates) {
    pred = predict(m_predator, newdata = d_pred_rich, se.fit = TRUE) 
    d_pred_rich = d_pred_rich %>%
      mutate(
        fit = pred$fit, se = pred$se.fit, low  = fit - 1.96 * se, high = fit + 1.96 * se,
        # Back-transform Poisson response
        rich = fit, rich_low  = low, rich_high = high,
        # Include bibi on original proportional scale
        bibi = bibi_grid
      )
  } else {
    pred = predict(m_predator, newdata = d_pred_rich, type = "link", se.fit = TRUE) 
    d_pred_rich = d_pred_rich %>%
      mutate(
        fit = pred$fit, se  = pred$se.fit, low  = fit - 1.96 * se, high = fit + 1.96 * se,
        # Back-transform Poisson response
        rich = exp(fit), rich_low  = exp(low), rich_high = exp(high),
        # Include bibi on original proportional scale
        bibi = bibi_grid
      )
  }
  
  # Back-transform observed B-IBI
  d$bibi_pcnt = (plogis(d$bibi) * n - 0.5) / (n - 1) * 100
  
  # Plot
  if (use_msom_richness_estimates) {
    p_bird_bibi = ggplot(d_pred_rich, aes(x = bibi, y = rich)) +
      geom_ribbon(aes(ymin = rich_low, ymax = rich_high), fill = "orange", alpha = 0.2) +
      geom_line(color = "orange", size = 1) +
      geom_point(data = d, aes(x = bibi_pct_orig, y = invert_predator_rich_mean), color = "orange", alpha = 0.4) +
      labs(x = "B-IBI", y = "Estimated predator richness") +
      theme_sleek() + theme(aspect.ratio = 1)
  } else {
    p_bird_bibi = ggplot(d_pred_rich, aes(x = bibi, y = rich)) +
      geom_ribbon(aes(ymin = rich_low, ymax = rich_high), fill = "orange", alpha = 0.2) +
      geom_line(color = "orange", size = 1) +
      geom_point(data = d, aes(x = bibi_pct_orig, y = rich_predator), color = "orange", alpha = 0.4) +
      labs(x = "B-IBI", y = "Observed predator richness") +
      theme_sleek() + theme(aspect.ratio = 1)
  }
  
  ## Plot all
  print(p_canopy_imp + p_bibi_imp + p_bird_bibi)
}

