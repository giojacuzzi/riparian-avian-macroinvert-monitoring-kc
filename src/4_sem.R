# 4_sem.R ===================================================================
# Fit piecewise structural equation models
#
# Input
exclude_agri_sites = TRUE # exclude outlier agricultural sites
use_msom_richness_estimates = TRUE # use richness estimates from the msom instead of naive observed values
msom_path = "data/cache/models/reach_invert_predator.rds"
in_cache_detections = "data/cache/1_preprocess_agg_pam_data/detections_calibrated_0.5.rds" # detections_calibrated_0.75.rds
# Output
out_cache_dir = "data/cache/4_sem"

source("src/global.R")

# Load site variable data ------------------------------------------------------------

site_data_reach = readRDS("data/cache/3_calculate_vars/site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/site_data_5000m.rds")

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

# NOTE: Very few aerial specialist foragers
eltontraits %>% filter(for_strat_aerial >= 10) %>% pull(common_name)

# Get AVONET invertivore community subset
# "Invertivore = species obtaining at least 60% of food resources from invertebrates in terrestrial systems, including insects, worms, arachnids, etc."
invertivores_avonet = species_traits %>% filter(trophic_niche == "Invertivore") %>% pull(common_name) %>% sort()
# "Aquatic Predator = species obtaining at least 60% of food resources from vertebrate and invertebrate animals in aquatic systems, including fish, crustacea, molluscs, etc."
species_traits %>% filter(trophic_niche == "Aquatic predator") %>% pull(common_name)

# Get Eltontraits invertivore community subset
# "Percent use of: Invertebrates-general, aquatic invertebrates, shrimp, krill, squid, crustacaeans, molluscs, cephalapod, polychaetes, gastropods, orthoptera, terrestrial Invertebrates, ground insects, insect larvae, worms, orthopterans, flying insects"
species_invert_ETdietGt10 = species_traits %>% filter(diet_inv >= 10) %>% pull(common_name) %>% sort() # nearly all species
# "Assignment to the dominant among five diet categories based on the summed scores of constituent individual diets."
invertivores_eltontraits = species_traits %>% filter(diet_5cat == "Invertebrate") %>% pull(common_name) %>% sort()

setdiff(invertivores_eltontraits, invertivores_avonet)
setdiff(invertivores_avonet, invertivores_eltontraits)

# Insectivores
sp_invert = species_traits %>% filter(diet_inv >= 10) %>% pull(common_name) %>% sort() # most species
sp_invert_primary = species_traits %>% filter(diet_5cat == "Invertebrate") %>% pull(common_name) %>% sort()

# Foraging guild: Aerial insectivores (e.g. swallows, swifts, flycatchers)
# primarily capture prey while they are flying in the air
sp_g_aerial_invert = species_traits %>% filter(foraging_guild_cornell %in% c("aerial forager", "flycatching")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

sp_g_aerial_invert_primary = species_traits %>% filter(foraging_guild_cornell %in% c("aerial forager", "flycatching")) %>%
  filter(common_name %in% sp_invert_primary) %>% pull(common_name)

sp_g_foliage_invert_primary = species_traits %>% filter(foraging_guild_cornell %in% c("foliage gleaner")) %>%
  filter(common_name %in% sp_invert_primary) %>% pull(common_name)

sp_g_ground_invert_primary = species_traits %>% filter(foraging_guild_cornell %in% c("ground forager")) %>%
  filter(common_name %in% sp_invert_primary) %>% pull(common_name)

sp_g_bark_invert_primary = species_traits %>% filter(foraging_guild_cornell %in% c("bark forager")) %>%
  filter(common_name %in% sp_invert_primary) %>% pull(common_name)

# Foraging guild: Foliage gleaners (e.g. warblers, vireos)
# typically capture insects located on vegetation or woody substrate
sp_g_foliage_invert = species_traits %>% filter(foraging_guild_cornell %in% c("foliage gleaner")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Foraging guild: Ground foragers (e.g. american robin)
# procure prey within forest leaf-litter and at the soil surface
sp_g_ground_invert = species_traits %>% filter(foraging_guild_cornell %in% c("ground forager")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Foraging guild: Bark-probers (e.g. brown creeper, red-breated nuthatch, woodpeckers and sapsuckers)
# extract prey from under bark or by boring into wood
sp_g_bark_invert = species_traits %>% filter(foraging_guild_cornell %in% c("bark forager")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Riparian associates
sp_ripasso = species_traits %>% filter(rip_asso_rich2002 == "X") %>% pull(common_name) %>% sort()
# Riparian obligats
sp_ripobl  = species_traits %>% filter(rip_obl_rich2002 == "X")  %>% pull(common_name) %>% sort()

# A priori list of species that:
# - Are primarily invertivores
# - Are reported to forage on aquatic macroinvertebrate taxa measured in B-IBI
#
# https://pugetsoundstreambenthos.org/About-BIBI.aspx
# B-IBI is calculated from stream site samples of the following insect orders...
# - Ephemeroptera, mayflies
# - Plecoptera, stoneflies
# - Trichoptera, caddisflies
# - Diptera, true flies and midges (with aquatic larvae e.g. Chironomidae, Simuliidae, Tipulidae)
# ... and classes of mollusks:
# - Bivalvia, freshwater clams (e.g. Unionidae, Sphaeriidae)
# - Gastropoda, freshwater snails (e.g. Lymnaeidae, Planorbidae)
# - aquatic worms
#
# Diets for each species obtained from Birds of the World, unless otherwise noted in raw data
sp_predator = species_traits %>% filter(invert_predator == "invert_predator") %>% pull(common_name) %>% sort()

# TODO: FINALIZE VALIDATIONS
sp_predator = setdiff(sp_predator, c(
  "black-capped chickadee", 
  "chestnut-backed chickadee"
))
#####

# Riparian associate invertivores
sp_ripasso_inv = c(sp_invert[sp_invert %in% c(sp_ripasso)])
sp_ripobl_inv  = c(sp_invert[sp_invert %in% c(sp_ripobl)])

# Load data for multi-species occupancy model --------------------------------------------------
if (use_msom_richness_estimates) {
  message("Loading data for multi-species occupancy model ", msom_path)
  msom_data = readRDS(msom_path)
  msom_summary = msom_data$msom_summary
  msom         = msom_data$msom
  groups       = msom_data$groups %>% arrange(common_name)
  sites        = msom_data$sites
  species      = msom_data$species
  
  ## Calculate estimated group richness per site
  z = msom_data$msom$sims.list$z
  group_idx = groups$group_idx
  group_names = groups %>% distinct(group_idx, .keep_all = TRUE) %>% arrange(group_idx) %>% pull(group)
  samples = dim(z)[1]
  J = dim(z)[2]
  I = dim(z)[3]
  G = max(group_idx)
  
  rich_group = array(NA, c(samples, J, G))
  for (g in 1:G) {
    species_in_g = which(group_idx == g)
    rich_group[ , , g] = apply(z[ , , species_in_g, drop = FALSE], c(1,2), sum)
  }
  rich_group_mean  = apply(rich_group, c(2,3), mean)
  rich_group_lower = apply(rich_group, c(2,3), quantile, probs = 0.025)
  rich_group_upper = apply(rich_group, c(2,3), quantile, probs = 0.975)
  msom_richness_estimates <- lapply(1:G, function(g) {
    group = group_names[g]
    tibble(
      site_id = sites
    ) %>%
      dplyr::mutate(
        !!paste0(group, "_rich_mean")  := rich_group_mean[, g],
        !!paste0(group, "_rich_lower") := rich_group_lower[, g],
        !!paste0(group, "_rich_upper") := rich_group_upper[, g]
      )
  })
  names(msom_richness_estimates) = group_names
  
  rich_total = apply(z, c(1, 2), sum)
  msom_richness_estimates$total = tibble(
    site_id    = sites,
    total_rich_mean  = apply(rich_total, 2, mean),
    total_rich_lower = apply(rich_total, 2, quantile, probs = 0.025),
    total_rich_upper = apply(rich_total, 2, quantile, probs = 0.975)
  )
}

# if (!dir.exists(dirname(path_out))) dir.create(dirname(path_out), recursive = TRUE)
# saveRDS(msom_richness_estimates, file = path_out)
# message(crayon::green("Cached model and results to", path_out))

# Exclude certain sites from analysis ------------------------------------------------

# Inspect total detections per site
# total_detections_by_site = detections$long %>% group_by(site_id) %>%
#   summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
#   arrange(desc(total_detections))
# ggplot(total_detections_by_site, aes(x = reorder(site_id, total_detections), y = total_detections)) + geom_col()

# Exclude sites 257 and 259 that are dominated by agriculture
if (exclude_agri_sites) {
  sites_to_exclude = c("257", "259")
  message("Excluding agricultural site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
  site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
  site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
  detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)
}

# Exclude sites with incomplete surveys
sites_to_exclude = setdiff(unique(site_data_reach$site_id), unique(detections$long$site_id))
sites_to_exclude = c(sites_to_exclude, "150", "3097") # ARU at sites 150 and 3097 destroyed mid-survey by water damage
message("Excluding incomplete survey site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

# Exclude one of site 155 / 159 to prevent spatial autocorrelation
sites_to_exclude = c(sites_to_exclude, "155") # ARU at sites 150 and 3097 destroyed mid-survey by water damage
message("Excluding spatially correlated survey site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

message("Retaining ", length(unique(site_data_reach$site_id)), " sites")

# Visualize joint data ----------------------------------------------------------------
message("Visualizing joint data")

presence_absence = detections$long %>% group_by(site_id, common_name) %>%
  summarise(presence = if_else(sum(n_detections, na.rm = TRUE) > 0, 1, 0), .groups = "drop")

# Summarize richness of different groups by site
site_group_richness =  presence_absence %>%
  group_by(site_id) %>%
  summarise(
    rich_all         = sum(presence),
    rich_inv         = sum(presence[common_name %in% sp_invert]),
    rich_inv_primary = sum(presence[common_name %in% sp_invert_primary]),
    # Foraging guilds
    rich_aerial_inv  = sum(presence[common_name %in% sp_g_aerial_invert]),
    rich_foliage_inv = sum(presence[common_name %in% sp_g_foliage_invert]),
    rich_ground_inv  = sum(presence[common_name %in% sp_g_ground_invert]),
    rich_bark_inv    = sum(presence[common_name %in% sp_g_bark_invert]),
    rich_aerial_inv_primary  = sum(presence[common_name %in% sp_g_aerial_invert_primary]),
    rich_foliage_inv_primary  = sum(presence[common_name %in% sp_g_foliage_invert_primary]),
    rich_ground_inv_primary  = sum(presence[common_name %in% sp_g_ground_invert_primary]),
    rich_bark_inv_primary  = sum(presence[common_name %in% sp_g_bark_invert_primary]),
    rich_predator     = sum(presence[common_name %in% sp_predator]),
    # Riparian habitat association
    rich_ripasso_inv = sum(presence[common_name %in% sp_ripasso_inv]),
    rich_ripasso     = sum(presence[common_name %in% sp_ripasso]),
    rich_ripobl_inv  = sum(presence[common_name %in% sp_ripobl_inv]),
    rich_ripobl      = sum(presence[common_name %in% sp_ripobl])
  )

# Incorporate mean msom richness estimates
if (use_msom_richness_estimates) {
  site_group_richness = left_join(site_group_richness, msom_richness_estimates$invert_predator, by = "site_id")
  site_group_richness = left_join(site_group_richness, msom_richness_estimates$`NA`, by = "site_id")
  site_group_richness = left_join(site_group_richness, msom_richness_estimates$total, by = "site_id")
}

# TODO: Do all exploratory analyses below with mean richness estimates

## Indicator species analysis
library(indicspecies)
spmat = presence_absence %>% pivot_wider(names_from = common_name, values_from = presence)
d = left_join(spmat, site_data_reach, by = "site_id")
spmat = spmat %>% select(-site_id) %>% as.data.frame()
bibi_excellent = ifelse(d$bibi >= 80, 1, 0)
indval = multipatt(spmat, bibi_excellent, control = how(nperm=999)) 
summary(indval)

# Richness across AVONET guilds per site
richness_by_trophic_niche = presence_absence %>% left_join(species_traits, by = "common_name") %>%
  group_by(site_id, trophic_niche) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_primary_lifestyle = presence_absence %>% left_join(species_traits, by = "common_name") %>%
  group_by(site_id, primary_lifestyle) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_ripasso = presence_absence %>% mutate(ripasso = common_name %in% c(sp_ripasso)) %>% group_by(site_id, ripasso) %>% summarize(count = sum(presence))
richness_by_predator = presence_absence %>% mutate(predator = common_name %in% c(sp_predator)) %>% group_by(site_id, predator) %>% summarize(count = sum(presence))

ggplot(richness_by_trophic_niche, aes(x = count, y = reorder(site_id, count), fill = trophic_niche)) + geom_col()
ggplot(richness_by_primary_lifestyle, aes(x = count, y = reorder(site_id, count), fill = primary_lifestyle)) + geom_col()
ggplot(richness_by_ripasso, aes(x = count, y = reorder(site_id, count), fill = ripasso)) + geom_col()
ggplot(richness_by_predator, aes(x = count, y = reorder(site_id, count), fill = predator)) + geom_col()

# Join with site data
site_data_reach = full_join(site_data_reach, site_group_richness, by = "site_id")
site_data_basin = full_join(site_data_basin, site_group_richness, by = "site_id")

# Richness as a function of different predictors
ggplot(site_data_reach, aes(x = rast_usfs_canopycover_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Canopy cover") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Reach imperviousness") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_basin, aes(x = rast_nlcd_impervious_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Basin imperviousness") + geom_text_repel(aes(label = site_id))

ggplot(site_data_reach, aes(x = bibi, y = rich_inv)) + geom_point() + geom_smooth(method = "lm") +
  labs(title = "BIBI") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_reach, aes(x = bibi, y = rich_aerial_inv)) + geom_point() + geom_smooth(method = "lm") +
  labs(title = "BIBI") + geom_text_repel(aes(label = site_id))

ggplot(site_data_reach, aes(x = (density_roads_paved), y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Density paved roads") + geom_text_repel(aes(label = site_id)) 

mapview(site_data_reach %>% select(site_id, rich_inv), zcol = "rich_inv")

# Richness of guilds as a function of...
ggplot(left_join(richness_by_trophic_niche, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = trophic_niche, fill = trophic_niche)) +
  geom_point() + geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_ripasso, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = ripasso, fill = ripasso)) +
  geom_point() + geom_smooth(aes(group = ripasso), method = "lm", se = FALSE)
ggplot(left_join(richness_by_ripasso, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = ripasso, fill = ripasso)) +
  geom_point() + geom_smooth(aes(group = ripasso), method = "lm", se = FALSE)

ggplot(left_join(richness_by_predator, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = predator, fill = predator)) +
  geom_point() + geom_smooth(aes(group = predator), method = "lm", se = FALSE)
ggplot(left_join(richness_by_predator, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = predator, fill = predator)) +
  geom_point() + geom_smooth(aes(group = predator), method = "lm", se = FALSE)

site_group_richness_long = site_group_richness %>% left_join(site_data_reach) %>% pivot_longer(
  cols = starts_with("rich_"),   # select all richness columns
  names_to = "richness_type",
  values_to = "richness_value"
)
ggplot(site_group_richness_long, aes(x = bibi, y = richness_value, color = richness_type)) +
  geom_point() +
  facet_wrap(~richness_type) +
  geom_smooth(aes(group = richness_type), method = "lm", se = FALSE)

# Species-specific presence/absence as a function of BIBI
ggplot(left_join(presence_absence %>% filter(common_name == "wilson's warbler"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "black-throated gray warbler"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "western wood-pewee"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

ggplot(left_join(presence_absence %>% filter(common_name == "violet-green swallow"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "pacific wren"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "swainson's thrush"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "black-headed grosbeak"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "golden-crowned kinglet"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

# Urban adapted species
ggplot(left_join(presence_absence %>% filter(common_name == "bewick's wren"), site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

ggplot(left_join(presence_absence %>% filter(common_name == "song sparrow"), site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

# Structural equation modeling -----------------------------------------------------------------------------

# Calculate pairwise collinearity among predictors
pairwise_collinearity = function(vars, threshold = 0.7) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

{
  # 550 m represents riparian zone within the local reach, and the 90% dispersal distance
  data_reach   = site_data_reach
  # 5 km represents the catchment landscape (roughly basin)
  data_basin = site_data_basin
  d = data.frame(
    "rich_all"         = data_reach$rich_all,
    "rich_inv"         = data_reach$rich_inv,
    "rich_inv_primary" = data_reach$rich_inv_primary,
    "rich_aerial_inv"  = data_reach$rich_aerial_inv,
    "rich_foliage_inv" = data_reach$rich_foliage_inv,
    "rich_ground_inv"  = data_reach$rich_ground_inv,
    "rich_bark_inv"    = data_reach$rich_bark_inv,
    "rich_aerial_inv_primary"  = data_reach$rich_aerial_inv_primary,
    "rich_foliage_inv_primary"  = data_reach$rich_foliage_inv_primary,
    "rich_ground_inv_primary"  = data_reach$rich_ground_inv_primary,
    "rich_bark_inv_primary"  = data_reach$rich_bark_inv_primary,
    "rich_predator"        = data_reach$rich_predator,
    "rich_ripasso"     = data_reach$rich_ripasso,
    "rich_ripobl"      = data_reach$rich_ripobl,
    "rich_ripasso_inv"     = data_reach$rich_ripasso_inv,
    "rich_ripobl_inv"      = data_reach$rich_ripobl_inv,
    "bibi"          = data_reach$bibi,
    "imp_reach"     = data_reach$rast_nlcd_impervious_sum_proportion,
    "imp_basin"     = data_basin$rast_nlcd_impervious_sum_proportion,
    "fhd_reach"     = data_reach$rast_gedi_fhd_mean,
    "fhd_basin"     = data_basin$rast_gedi_fhd_mean,
    "canopy_reach"  = data_reach$rast_usfs_canopycover_sum_proportion,
    "canopy_basin"  = data_basin$rast_usfs_canopycover_sum_proportion,
    "ed_reach"      = data_reach$edge_density,
    "ed_basin"      = data_basin$edge_density,
    "forest_reach"  = data_reach$nlcd_forest,
    "forest_basin"  = data_basin$nlcd_forest,
    "riphab_reach"  = data_reach$nlcd_forest_and_wetlands,
    "riphab_basin"  = data_basin$nlcd_forest_and_wetlands,
    "site_id"       = data_reach$site_id
  )
  if (use_msom_richness_estimates) {
    d = left_join(d, msom_richness_estimates$invert_predator, by = "site_id")
    d = left_join(d, msom_richness_estimates$`NA`, by = "site_id")
    d = left_join(d, msom_richness_estimates$total, by = "site_id")
  }
  pairwise_collinearity(d %>% select(where(is.numeric)))

  # Exploratory models
  if (FALSE) {
    
    m_bibi = lm(bibi ~ canopy_reach + imp_basin, d)
    
    # All species
    m_all = glm(rich_all ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_all); plot(sem); print(summary(sem))

    # Riparian associates/obligates
    m_ripasso_inv = glm(rich_ripasso_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_ripasso_inv); plot(sem); print(summary(sem))
    
    m_ripobl_inv = glm(rich_ripobl_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_ripobl_inv); plot(sem); print(summary(sem))
    
    # All invertivores (primary and >10% diet)
    m_inv = glm(rich_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_inv); plot(sem); print(summary(sem))
    
    m_inv_primary = glm(rich_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_inv_primary); plot(sem); print(summary(sem))
    
    # Invertivorous foraging guilds (>= 10% invertivores)
    m_aerial_inv = glm(rich_aerial_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_aerial_inv); plot(sem); print(summary(sem))
  
    m_foliage_inv = glm(rich_foliage_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_foliage_inv); plot(sem); print(summary(sem))
    
    m_ground_inv = glm(rich_ground_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_ground_inv); plot(sem); print(summary(sem))
    
    m_bark_inv = glm(rich_bark_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_bark_inv); plot(sem); print(summary(sem))
    
    # Primary invertivores
    m_aerial_inv_primary = glm(rich_aerial_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_aerial_inv_primary); plot(sem); print(summary(sem))
    
    m_foliage_inv_primary = glm(rich_foliage_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_foliage_inv_primary); plot(sem); print(summary(sem))
    
    m_ground_inv_primary = glm(rich_ground_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_ground_inv_primary); plot(sem); print(summary(sem))
    
    m_bark_inv_primary = glm(rich_bark_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_bark_inv_primary); plot(sem); print(summary(sem))
  }
  
  ## SEM for invertivore guild ==========================================================
  # (species that are primarily insectivorous and are reported to predate on aquatic inverts)
  n = nrow(d)
  
  # Canopy reach (logit transformation for proportional bounds) as a function of reach-scale impervious %
  d$canopy_reach = qlogis(d$canopy_reach)
  m_canopy = lm(canopy_reach ~ imp_reach, d)
  
  # B-IBI (logit transformation for proportional bounds) as a function of watershed-scale impervious % and reach-scale canopy %
  d$bibi = d$bibi / 100
  d$bibi = (d$bibi * (n - 1) + 0.5) / n # Smithson & Verkuilen (2006) squeeze
  d$bibi = qlogis(d$bibi)
  m_bibi = lm(bibi ~ canopy_reach + imp_basin, data = d)
  
  # Predator group richness as a function of bibi, reach canopy, and reach impervious
  if (!use_msom_richness_estimates) {
    message("Using naive observed richness as response variable")
    m_predator = glm(rich_predator ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    # Fit piecewise SEM with component models
    sem_predator = psem(m_canopy, m_bibi, m_predator); plot(sem_predator); print(summary(sem_predator))
    
  } else {
    message("Using msom richness mean as response variable")
    d_msom = d
    d_msom$rich_predator = d$invert_predator_rich_mean
    # NOTE: piecewiseSEM does not support glm Gamma distribution; must use normal distribution
    m_predator = lm(rich_predator ~ bibi + canopy_reach + imp_reach, d_msom)
    sem_predator = psem(m_canopy, m_bibi, m_predator); plot(sem_predator); print(summary(sem_predator))
    coefs_predator = coefs(sem_predator) %>% clean_names() %>% rename(signif = x)
    # Test to support assumption of normality
    shapiro.test(residuals(m_predator)) # p > 0.05 => approximately normally distributed residuals
    
    # "Other" group richness
    d_msom$rich_NA = d$NA_rich_mean
    m_NA = lm(rich_NA ~ bibi + canopy_reach + imp_reach, d_msom)
    sem_NA = psem(m_canopy, m_bibi, m_NA); plot(sem_NA); print(summary(sem_NA))
    coefs_NA = coefs(sem_NA) %>% clean_names() %>% rename(signif = x)
    
    # Total group richness
    d_msom$rich_total = d$total_rich_mean
    m_total = lm(rich_total ~ bibi + canopy_reach + imp_reach, d_msom)
    sem_total = psem(m_canopy, m_bibi, m_total); plot(sem_total); print(summary(sem_total))
    coefs_total = coefs(sem_total) %>% clean_names() %>% rename(signif = x)
  }
  
  if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
  out_filepath = file.path(out_cache_dir, paste0("sem_predator.rds"))
  saveRDS(sem_predator, out_filepath)
  message(crayon::green("Cached", out_filepath))
  
  # TODO: Also do full community model
  
  # Propagate uncertainty from MSOM --------------------------------------------------------
  if (use_msom_richness_estimates) {
    # Fit sem over all posterior draws to propagate uncertainty
    for (g in 1:(G+1)) {
      if (g > G) {
        group_name = "total"
      } else {
        group_name = group_names[g]
      }
      message("Propagating MSOM uncertainty for '", group_name, "' group richness response")
      z = msom$sims.list$z
      draws = dim(z)[1]
      stopifnot(all(groups$common_name == msom_data$species))
      coeffs_draws = tibble()
      pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = draws, clear = FALSE)
      for (draw in 1:draws) {
        # Calculate estimated group richness at each site
        if (group_name == "total") {
          rich_group_draw = data.frame(
            site_id = msom_data$sites,
            rich_group_draw = rich_group[draw, , 1] + rich_group[draw, , 2]
          )
        } else {
          rich_group_draw = data.frame(
            site_id = msom_data$sites,
            rich_group_draw = rich_group[draw, , g]
          )
        }
        d_msom = d %>% left_join(rich_group_draw, by = "site_id")
        
        # Fit the component model and SEM
        m_rich_draw = lm(rich_group_draw ~ bibi + canopy_reach + imp_reach, d_msom)
        sem_draw = psem(m_canopy, m_bibi, m_rich_draw)
        
        # Extract SEM coefficients and store
        coeffs_draw = coefs(sem_draw) %>% clean_names() %>%
          filter(response == "rich_group_draw") %>%
          mutate(draw = draw) %>% rename(signif = x)
        coeffs_draws = rbind(coeffs_draws, coeffs_draw)
        pb$tick()
      }
      rich_group_coefs = coeffs_draws %>%
        group_by(predictor) %>%
        summarise(
          mean  = mean(std_estimate),
          lower = quantile(std_estimate, 0.025),
          upper = quantile(std_estimate, 0.975),
          response = "rich_group_draw",
          .groups = "drop"
        ) %>% rename(std_estimate = mean)
      if (group_name == "invert_predator") {
        coefs_to_compare = coefs_predator
      } else if (group_name == "NA") {
        coefs_to_compare = coefs_NA
      } else if (group_name == "total") {
        coefs_to_compare = coefs_total
      }
      all_coefs = bind_rows(rich_group_coefs, coefs_to_compare)
      all_coefs$group = group_name
      p = ggplot(all_coefs, aes(x = std_estimate, y = interaction(predictor, response))) +
        geom_vline(xintercept = 0, color = "gray") +
        geom_errorbar(aes(xmin = lower, xmax = upper), width = 0) +
        geom_point() + labs(title = group_name); print(p)
      
      message("SEM coefficients including propagated uncertainty:")
      print(all_coefs)
    }
  }
  
  if (FALSE) { # TODO: Alternative scale hypotheses validations (compare AIC and R.squared)
    m_reach_canopy = lm(canopy_reach ~ imp_reach, d)
    m_reach_bibi   = lm(bibi ~ canopy_reach + imp_reach, d)
    m_reach_rich   = glm(rich_predator ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem_reach = psem(m_reach_canopy, m_reach_bibi, m_reach_rich); plot(sem_reach); print(summary(sem_reach))
    
    d$canopy_basin = qlogis(d$canopy_basin)
    m_basin_canopy = lm(canopy_basin ~ imp_basin, d)
    m_basin_bibi   = lm(bibi ~ canopy_basin + imp_basin, d)
    m_basin_rich   = glm(rich_predator ~ bibi + canopy_basin + imp_basin, d, family = poisson)
    sem_basin = psem(m_basin_canopy, m_basin_bibi, m_basin_rich); plot(sem_basin); print(summary(sem_basin)) 
  }
  
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
}

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
