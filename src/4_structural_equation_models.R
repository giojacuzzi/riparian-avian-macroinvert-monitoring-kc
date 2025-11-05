library(tidyverse)
library(sf)
library(lme4)
library(DHARMa)
library(ggrepel)
library(patchwork)
library(piecewiseSEM)
theme_set(theme_minimal())

files = list.files("data/cache/3_structural_equation_model", pattern = "^site_data.*\\.rds$", full.names = TRUE)

site_data = files %>%
  lapply(function(f) {
    df <- readRDS(f)
    f <- tools::file_path_sans_ext(basename(f))
    df$variant <- f
    df$variant_m <- str_match(f, "_(\\d+)m")[,2]  # captures digits before "m_"
    df$variant_t <- str_match(f, "_(\\d*\\.?\\d+)t(?:_|$)")[,2]  # captures digits before "t_"
    df <- df %>% mutate(across(where(~inherits(.x, "units")), as.numeric))
    df
  }) %>%
  bind_rows()
site_data$variant_m = factor(site_data$variant_m, levels = sort(as.numeric(unique(site_data$variant_m))))
site_data$variant_t = factor(site_data$variant_t, levels = sort(as.numeric(unique(site_data$variant_t))))

# ggplot(site_data %>% filter(variant == "site_data_550m_0.75t_3d"),
#        aes(x = rast_usfs_canopycover_sum_proportion, y = richness_invert_avo)) +
#   geom_point() + geom_smooth(method = "lm")
{
  y_var = "richness_invert_avo" # richness_invert_avo, richness, bibi, 
  x_vars = c(
    "rast_nlcd_impervious_sum_proportion",
    "rast_usfs_canopycover_sum_proportion",
    "bibi",
    "nlcd_forest",
    "nlcd_wetlands",
    "rast_lemma_age_mean",
    "rast_lemma_denall_mean",
    "rast_lemma_qmd_mean",
    "rast_gedi_height_cv",
    "rast_gedi_fhd_mean",
    "rast_lemma_ba_mean",
    "density_roads_paved"
  )
  cor_results = map_dfr(x_vars, function(x) {
    site_data %>% st_drop_geometry() %>%
      group_by(variant_m, variant) %>%
      summarise(
        x_var = x,
        r = cor(.data[[x]], .data[[y_var]], use = "complete.obs"),
        r2 = r^2,
        .groups = "drop"
      )
  }) %>% arrange(desc(r2))
  cor_results$variant_t = factor(str_match(cor_results$variant, "_(\\d*\\.?\\d+)t(?:_|$)")[,2])
  cor_results$variant_m = factor(cor_results$variant_m)
  
  ggplot(cor_results, aes(x = variant_t, y = r, fill = factor(variant_m))) +
    geom_col(position = position_dodge(width = 0.8)) + scale_fill_viridis_d() +
    facet_wrap(~ reorder(x_var, -abs(r)), scales = "free_x") +
    labs(title = paste0("Univariate predictors for ", y_var),
      x = "Confidence threshold", y = "Value", fill = "Buffer (m)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

{
  plots = map(x_vars, function(x) {
    ggplot(site_data,
           aes(x = .data[[x]],
               y = .data[[y_var]],
               color = factor(variant_t))) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = "lm", se = FALSE) +
      facet_wrap(~ variant_m) +
      labs(x = x, y = y_var, color = "Threshold") 
  })
  names(plots) = x_vars
  plots["rast_nlcd_impervious_sum_proportion"]
}

# Data exploration -------------------------------------------------------------------

variant_to_use = "site_data_550m_0.75t_3d"

cor_results %>% filter(variant == variant_to_use) %>% select(x_var, r, r2)

data = site_data %>% filter(variant == variant_to_use)

y_var = "richness_invert_avo"

# Richness as a function of different predictors
p_cr = ggplot(data, aes(x = rast_usfs_canopycover_sum_proportion, y = !!sym(y_var))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, color = "forestgreen") + geom_text_repel(aes(label = site_id), color = "gray") +
  labs(title = "Canopy cover")

p_ir = ggplot(data, aes(x = rast_nlcd_impervious_sum_proportion, y = !!sym(y_var))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, color = "forestgreen") + geom_text_repel(aes(label = site_id), color = "gray") +
  labs(title = "Impervious % (interaction zone)")

p_br = ggplot(data, aes(x = bibi, y = !!sym(y_var))) +
  geom_point() + geom_smooth(method = "lm", se = FALSE, color = "forestgreen") + geom_text_repel(aes(label = site_id), color = "gray") +
  labs(title = "BIBI")

p_cb = ggplot(data, aes(x = rast_usfs_canopycover_sum_proportion, y = bibi)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + geom_text_repel(aes(label = site_id), color = "gray") +
  labs(title = "Canopy cover")

p_ib = ggplot(data, aes(x = rast_nlcd_impervious_sum_proportion, y = bibi)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) + geom_text_repel(aes(label = site_id), color = "gray") +
  labs(title = "Impervious % (interaction zone)")

(p_cr + p_ir + p_br) / (p_cb + p_ib + plot_spacer())


# Structural equation modeling -------------------------------------------------------

candidate_vars = c(
  ## Predictors
  "area_int"             = "area_interaction_zone",
  "area_rip"             = "area_riparian_usfs",
  "bibi"                 = "bibi",
  "site"                 = "site_id",
  "basin"                = "basin_name",
  # "elevation"            = "elevation",
  "imperviousness_basin" = "basin_impervious",
  "imperviousness"       = "rast_nlcd_impervious_sum_proportion",
  "abund_dev_varint"     = "nlcd_developed_variable_intensity",
  "abund_dev_opensp"     = "nlcd_developed_open_space",
  "abund_forest"         = "nlcd_forest",
  "abund_wetland"        = "nlcd_wetlands",
  "abund_openwater"      = "nlcd_open_water",
  "canopy_usfs"          = "rast_usfs_canopycover_sum_proportion",
  "canopy_gedi"          = "rast_gedi_cover_sum",
  # "height_landfire"      = "rast_landfire_treeheight_mean",
  "height_gedi"          = "rast_gedi_height_mean",
  # "height_cv_landfire"   = "rast_landfire_treeheight_cv",
  "height_cv_gedi"       = "rast_gedi_height_cv",
  "height_cv_gedi_mask"  = "rast_gedi_height_impmask_cv",
  "fhd"                  = "rast_gedi_fhd_mean",
  "pavd20m"              = "rast_gedi_pavd20m_mean",
  "pavd5to10m"           = "rast_gedi_pavd5to10m_mean",
  "stand_age_mean"       = "rast_lemma_age_mean",
  "stand_age_cv"         = "rast_lemma_age_cv",
  "stand_ba"             = "rast_lemma_ba_mean",
  "stand_denall"         = "rast_lemma_denall_mean",
  "stand_dencon"         = "rast_lemma_dencon_mean",
  "stand_denhw"          = "rast_lemma_denhw_mean",
  "stand_qmd"            = "rast_lemma_qmd_mean",
  "edge_density"         = "edge_density",
  "aggregation"          = "aggregation",
  "tree_richness"        = "tree_species_richness",
  # "dist_road_major"      = "dist_road_major",
  # "dist_road_paved"      = "dist_road_paved",
  "density_roads_major"  = "density_roads_major",
  "density_roads_paved"  = "density_roads_paved", 
  ## Responses
  "richness"             = "richness",
  "richness_invert_et"   = "richness_invert_et",
  "richness_invert_avo"  = "richness_invert_avo",
  "richness_ripdep"      = "richness_ripdep"
)
d = as.data.frame(data) %>% select(all_of(candidate_vars))

# Any NA?
sort(colnames(data)[colSums(is.na(data)) > 0])

# Calculate pairwise collinearity among predictors
pairwise_collinearity = function(vars, threshold = 0.7) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

pairwise_collinearity(d %>% select(
  bibi,
  # imperviousness_basin,
  imperviousness,
  canopy_usfs,
  # abund_forest,
  # abund_dev_varint,
  # height_gedi,
  fhd,# height_cv_gedi,
  # height_cv_gedi_mask,
  # pavd20m, pavd5to10m,
  # elevation,
  # edge_density, aggregation,
  stand_age_mean #, stand_age_cv, stand_ba, stand_denall, stand_dencon, stand_denhw, stand_qmd, tree_richness
))

# TODO: VIF analyses for multicollinearity, e.g. sort(car::vif(m_bibi))

# Create and inspect structural equation model claims and fit
# - A significant independence claim from a test of directed separation suggests that the path is missing
# or misspecified.
# - A significant global Fisher’s C p-value (< 0.05) suggests that the modeled structure is statistically
# significantly different than the structure implied by the data, and that alternative pathways or causal
# links with missing variables warrant further exploration
{
  # Fit component regressions
  m1 = lm(bibi ~ imperviousness + canopy_usfs,
          d)
  m2 = glm(richness_invert_avo ~ bibi + imperviousness + canopy_usfs,
           d, family = poisson(link = "log"))
  sem = psem(m1, m2); plot(sem); print(summary(sem))
  
  # Check overdispersion -- if overdispersed, fit negative binomial
  # simres_pois = simulateResiduals(m2, n = 1000)
  # plot(simres_pois)
  # testDispersion(simres_pois)
  
  # Check spatial autocorrelation
  # coords = cbind(site_data$Easting, site_data$Northing)
  # testSpatialAutocorrelation(simres, x = coords[,1], y = coords[,2], plot = TRUE)
  
  # Bootstrap CIs for indirect effects?
}

{
  m1 = lm(bibi ~ imperviousness + canopy_usfs,
          data)
  m2 = glm(richness ~ imperviousness,
           data, family = poisson(link = "log"))
  sem = psem(m1, m2); plot(sem); print(summary(sem))
}

m2 = glm(richness ~ imperviousness, data, family = poisson(link = "log"))
plot(data$imperviousness_basin, data$richness_invert_avo)

# Workshop --------------------------------------------------------------------

# TODO: Area-weight covariates?

# Roads
sem_roads = psem(
  lm(bibi ~ imperviousness_local + canopy_usfs,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi_mask + (density_roads_paved),
      data, family = poisson(link = "log"))
); plot(sem_roads,
        # node_attrs = list(fontsize = 8, fontcolor = "black", fillcolor = c("lightblue", "yellow", "gray", "green", "green", "gray"))
); print(summary(sem_roads))
mapview(site_data %>% select(site_id, richness, richness_invert_avo, dist_road_major, dist_road_paved), zcol = "dist_road_major")

# Local
sem_local = psem(
  lm(bibi ~ imperviousness_local + canopy_usfs,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
); plot(sem_local); print(summary(sem_local))

# Local ripdep
sem_local = psem(
  lm(bibi ~ imperviousness_local + canopy_usfs,
     data),
  glm(richness_ripdep ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
); plot(sem_local); print(summary(sem_local))

# Comparison with overall richness
plot(psem(
  lm(bibi ~ imperviousness_local + canopy_usfs,
     data),
  glm(richness ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
))

# Explore local stand and landscape structure (Gradient Nearest Neighbor vegetation maps from USFS and OSU)
# NOTE: OVERPARAMETRIZED
alt = psem( # using stand structure
  lm(bibi ~ imperviousness_basin + canopy_usfs + stand_qmd + stand_denall + tree_richness,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + stand_qmd + stand_denall + height_cv_gedi_mask + tree_richness,
      data, family = poisson(link = "log"))
); plot(alt); summary(alt)
pairwise_collinearity(data %>% select(
  bibi, imperviousness_basin, canopy_usfs, stand_qmd, stand_denall, tree_richness, richness_invert_avo, height_cv_gedi_mask, tree_richness
))

alt = psem( # parsimonious
  lm(bibi ~ imperviousness_basin + canopy_usfs + stand_qmd,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + stand_qmd + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
); plot(alt); summary(alt)

# TODO:
# - riparian zone area? floodplain (lateral) connectivity?
# - tree density? woody debris? snags? understory vegetation?
# - proximity to forest on landscape (basin) level?

stop("Arrived at workshop")

plot(psem( # using stand age
  lm(bibi ~ imperviousness_basin + canopy_usfs + stand_age_mean,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + stand_age_mean + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
))

# NOTE: There are not enough observations per basin to support inclusion of a random effect to account for the nested structure of site within basin
table(data$basin)

# TODO: Explore promising common vegetation vars from literature:
# - Canopy cover (USFS canopy cover seems more accurate than GEDI)
# - Canopy height (gedi_height_mean)
#    - canopy height classes (e.g. average cover height class 4-9m, average cover height class 0-4m)
# - Foliage height diversity / gedi_height_cv
# Other vars:
# - pattern / quality / connectivity of riparian areas (mean patch size, connectivity, edge density)
# - stream network connectivity
# - geomorphology
# - elevation
# Other responses:
# - occupancy probability
# - native vs nonnative species

# - Effects of urbanization at multiple spatial scales?

# - Interactions?

# TESTING
plot(psem(
  lm(bibi ~ imperviousness_basin + canopy_usfs,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + height_cv_gedi + edge_density,
      data, family = poisson(link = "log"))
))

plot(psem(
  lm(bibi ~ imperviousness_basin + abund_forest,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + abund_forest + height_cv_gedi,
      data, family = poisson(link = "log"))
))

# Interaction
plot(psem(
  lm(bibi ~ imperviousness_basin + canopy_usfs,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + (height_cv_gedi * canopy_usfs),
      data, family = poisson(link = "log"))
))

# Multi-scale imperviousness
plot(psem(
  lm(bibi ~ imperviousness_basin + canopy_usfs,
     data),
  lm(imperviousness_local ~ imperviousness_basin,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi,
      data, family = poisson(link = "log"))
))

plot(psem(
  lm(bibi ~ abund_dev_varint + abund_forest,
     data),
  glm(richness_invert_avo ~ bibi + abund_dev_varint + abund_forest + abund_dev_opensp,
      data, family = poisson(link = "log"))
))

# Indicator species analysis to find species associated with high B-IBI sites?
