# 3_explore_data.R ===================================================================
# Load ARU and PSSB site locations, calculate all site covariates from land cover
# data, and do some exploratory structural equation modeling
#
# Inputs:
in_path_site_metadata   = "data/site_metadata.csv"
in_path_pssb_data       = "data/raw/pssb/ScoresByYear.csv"
in_path_nlcd_metadata   = "data/raw/nlcd_metadata.csv"
in_cache_detections     = "data/cache/1_preprocess_agg_pam_data/detections.rds"
in_cache_geospatial_dir = "data/cache/2_preprocess_geospatial_data"
in_path_species_list    = "data/pam/species_list.txt"
in_path_avonet_traits   = "data/traits/AVONET Supplementary dataset 1.xlsx"
in_path_elton_traits    = "data/traits/EltonTraits 1.0/BirdFuncDat.csv"

# Load required packages (automatically install any missing) -------------------------
pkgs = c(
  "tidyverse",   # data manipulation
  "janitor",     # data cleaning
  "sf",          # vector data
  "terra",       # raster data
  "tidyterra",   # raster data manipulation
  "tigris",      # political boundaries
  "geosphere",   # distance metrics
  "mapview",     # interactive geospatial visualization
  "leafsync",    # synchronized mapview panels
  "progress",    # dynamic progress bar
  "ggrepel",     # plot annotations
  "piecewiseSEM" # structural equation modeling
)
sapply(pkgs, function(pkg) {
  if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
  as.character(packageVersion(pkg))
})

# Load ARU and PSSB site locations and define study area -----------------------------
message("Loading ARU and PSSB site locations")

site_metadata = read_csv(in_path_site_metadata, show_col_types = FALSE) %>% clean_names() %>% mutate(site_id = as.character(site_id))

crs_standard = "EPSG:32610" # shared coordinate reference system (metric)

# Create sf points for site locations (reference original crs 4326)
sites_aru = site_metadata %>%
  st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326) %>% st_transform(crs = crs_standard)
sites_pssb = site_metadata %>%
  st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326) %>% st_transform(crs = crs_standard)

# Define study area by King County geopolitical boundary
study_area = counties(state = "WA", cb = TRUE) %>% filter(NAME == "King") %>% st_transform(crs = crs_standard)

# View study area and sites
mapview(study_area, alpha.regions = 0, lwd = 2, layer.name = "King County") +
  mapview(sites_aru, col.region = "green", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

# Calculate distance between paired ARU and PSSB sites
summary(site_metadata %>% rowwise() %>%
          mutate(dist_m = geosphere::distHaversine(c(long_aru,  lat_aru), c(long_pssb, lat_pssb))) %>% pull(dist_m))

# Load Puget Sound Stream Benthos monitoring data ------------------------------------
message("Loading PSSB BIBI data")

pssb_data = read_csv(in_path_pssb_data, show_col_types = FALSE) %>% clean_names() %>%
  mutate(site_id = as.character(site_id)) %>% select(site_id, x2024)

# Get the BIBI for the year in which the site was sampled
sites_aru = left_join(sites_aru, pssb_data, by = c("site_id"))
sites_aru = sites_aru %>% mutate(year = year(date_start)) %>% rowwise() %>% mutate(bibi = NA) %>%
  mutate(bibi = case_when(
    year == 2024 ~ x2024,
    year == 2025 ~ x2024, # TODO: 2025 data not yet available
    TRUE ~ bibi
  ))

# View sites by BIBI value
mapview(sites_aru, zcol = "bibi")

# Classify values into discrete categories
m = mapview(sites_aru, zcol = "bibi", cex = "bibi", at = c(0, 20, 40, 60, 80, 100),
        col.regions = c("red", "orange", "yellow", "green", "blue")); m

# Load and process geospatial data ---------------------------------------------------
message("Loading geospatial data")
theme_set(theme_minimal())

# Load cached raster data
rast_filepaths = list.files(in_cache_geospatial_dir, pattern = "^rast_.*\\.tif$", full.names = TRUE)
rast_data = lapply(rast_filepaths, rast)
names(rast_data) = gsub("\\.tif$", "", basename(rast_filepaths))

# Load cached vector
sf_ripfb = st_read(paste0(in_cache_geospatial_dir, "/sf_ripfb.gpkg"), quiet = TRUE)

# Load cached basin sf objects and retain only those sampled
sf_basins12d = st_read(paste0(in_cache_geospatial_dir, "/sf_basins12d.gpkg"), quiet = TRUE) %>%
  clean_names() %>% select(huc12, name, area_sq_km) %>% mutate(basin_name = name, basin_area = area_sq_km)
sf_basins12d = sf_basins12d %>% filter(lengths(st_intersects(., sites_aru)) > 0) # TODO: exclude lakes?

# Store subsequent data calculated per site in `site_data`
site_data = sites_aru

# Calculate imperviousness at the basin (landscape) scale ------------------------------------
message("Calculating imperviousness at the basin scale")

# Calculate mean imperviousness per basin
basin_impervious = terra::extract(rast_data$rast_nlcd_impervious, vect(sf_basins12d), fun = mean, na.rm = TRUE)
sf_basins12d = sf_basins12d %>% mutate(basin_impervious = basin_impervious[[2]])

# Store results
site_data = st_intersection(site_data, sf_basins12d)

# Calculate geospatial variables for all sites at local scale --------------------------------
buffer_size = 500 # ~500 m insect emergence 90% flux range falloff
message("Calculating geospatial variables for all sites at local scale (buffer size ", buffer_size, " m)")

# Load NLCD metadata for cover class codes
nlcd_metadata = read.csv(in_path_nlcd_metadata)

# Helper function to crop and mask raster to a buffer
crop_and_mask = function(rast, buff) { mask(crop(rast, buff), buff) }

# Helper function to calculate summary statistics for a continuous raster
rast_stats = function(rast, na.rm = TRUE) {
  stats = global(rast, fun = c("max", "min", "sum", "mean", "sd"), na.rm = na.rm)
  stats$cv = stats$sd / stats$mean
  return(stats)
}

geospatial_site_data = list()
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
for (s in 1:nrow(sites_aru)) {

  # Create a buffer around the site
  site = st_as_sf(sites_aru[s,])
  site_buffer = st_buffer(site, buffer_size)
  # mapview(site) + mapview(site_buffer)
  
  ## Calculate stats for continuous raster data
  
  # NLCD imperviousness
  site_nlcd_impervious  = crop_and_mask(rast_data$rast_nlcd_impervious, site_buffer)
  stats_nlcd_impervious = rast_stats(site_nlcd_impervious)
  # mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_nlcd_impervious)
  
  # NASA GEDI foliage height diversity
  site_gedi_fhd  = crop_and_mask(rast_data$rast_gedi_fhd, site_buffer)
  stats_gedi_fhd = rast_stats(site_gedi_fhd)
  # mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_gedi_fhd)
  
  # NASA GEDI canopy cover
  site_gedi_cover  = crop_and_mask(rast_data$rast_gedi_cover, site_buffer)
  stats_gedi_cover = rast_stats(site_gedi_cover)
  
  # NASA GEDI canopy height
  site_gedi_height  = crop_and_mask(rast_data$rast_gedi_height, site_buffer)
  stats_gedi_height = rast_stats(site_gedi_height)
  
  # NASA GEDI proportion of vegetation density (mature upper canopy)
  site_gedi_pavd20m  = crop_and_mask(rast_data$rast_gedi_pavd20m, site_buffer)
  stats_gedi_pavd20m = rast_stats(site_gedi_pavd20m)
  
  # NASA GEDI proportion of vegetation density
  site_gedi_pavd5to10m  = crop_and_mask(rast_data$rast_gedi_pavd5to10m, site_buffer)
  stats_gedi_pavd5to10m = rast_stats(site_gedi_pavd5to10m)
  
  # USFS canopy cover
  site_usfs_canopycover  = crop_and_mask(rast_data$rast_usfs_canopycover, site_buffer)
  stats_usfs_canopycover = rast_stats(site_usfs_canopycover)
  
  # Landfire tree height
  site_landfire_treeheight  = crop_and_mask(rast_data$rast_landfire_treeheight, site_buffer)
  stats_landfire_treeheight = rast_stats(site_landfire_treeheight)
  
  ## Calculate composition for categorical NLCD land cover raster data
  
  site_nlcd_landcover = crop_and_mask(rast_data$rast_nlcd_landcover, site_buffer)

  # Calculate relative abundance of each cover class
  freq_table = freq(site_nlcd_landcover) %>% select(value, count)
  freq_table$percent = (freq_table$count / sum(freq_table$count)) 
  cover_detail = full_join(freq_table %>% rename(class = value),
                           nlcd_metadata %>% select(class), by = "class")
  cover_detail = cover_detail %>% mutate(
    count = ifelse(is.na(count), 0, count),
    percent = ifelse(is.na(percent), 0, percent)
  ) %>% arrange(desc(percent))
  
  # Summarize by combining similar classes into broader groups
  cover_summary = cover_detail %>%
    mutate(group = case_when(
      grepl("Intensity", class)          ~ "Developed, Variable Intensity",
      grepl("Forest", class)             ~ "Forest",
      grepl("Wetlands", class)           ~ "Wetlands",
      grepl("Pasture|Cultivated", class) ~ "Agriculture",
      TRUE ~ class  # retain other classes as they are
    )) %>% group_by(group) %>%
    summarise(
      count = sum(count),
      percent = sum(percent)
    ) %>% arrange(desc(percent))
  
  # Store all data for the site
  rasters = list(
    site_nlcd_impervious,
    site_gedi_fhd,
    site_gedi_cover,
    site_gedi_height,
    site_gedi_pavd20m,
    site_gedi_pavd5to10m,
    site_usfs_canopycover,
    site_landfire_treeheight,
    site_nlcd_landcover
  )
  names(rasters) = unlist(lapply(rasters, names))
  raster_stats = do.call(rbind, list(
    stats_nlcd_impervious,
    stats_gedi_fhd,
    stats_gedi_cover,
    stats_gedi_height,
    stats_gedi_pavd20m,
    stats_gedi_pavd5to10m,
    stats_usfs_canopycover,
    stats_landfire_treeheight
  )) %>% rownames_to_column(var = "name")
  
  geospatial_site_data[[as.character(site$site_id)]] = list(
    rasters = rasters,
    raster_stats = raster_stats,
    nlcd_cover = list(
      "nlcd_cover_detail" = cover_detail,
      "nlcd_cover_summary" = cover_summary
    )
  )
  pb$tick()
}

###
# Explore landfire vs gedi height
debug_buffers = st_buffer(sites_aru, buffer_size)
debug_landfire = crop_and_mask(rast_data$rast_landfire_treeheight, debug_buffers)
debug_gedi = crop_and_mask(rast_data$rast_gedi_height, debug_buffers)
debug_forest = crop_and_mask(rast_data$rast_nlcd_landcover, debug_buffers)
debug_canopy = crop_and_mask(rast_data$rast_usfs_canopycover, debug_buffers)
mapview(debug_landfire) + mapview(debug_gedi) + mapview(debug_canopy) + mapview(debug_forest)
###

# Join land cover data with sites into a single table --------------------------------
message("Joining data for all sites")

# Convert geospatial_site_data into a single table and pivot wider
nlcd_summary = lapply(names(geospatial_site_data), function(id) {
  df = geospatial_site_data[[id]]$nlcd_cover$nlcd_cover_summary
  df$site_id = id
  return(df)
}) %>% bind_rows()
nlcd_summary = nlcd_summary %>%
  select(site_id, group, percent) %>%
  pivot_wider(names_from = group, values_from = percent, names_prefix = "nlcd_") %>%
  janitor::clean_names()

# Join land cover data
site_data = site_data %>% left_join(nlcd_summary, by = "site_id")

# Extract raster stats for each site in wide format
stats = lapply(names(geospatial_site_data), function(i) {
  df = geospatial_site_data[[i]]$raster_stats
  df$site_id = i
  return(df)
}) %>% bind_rows() %>% pivot_wider(
  id_cols = site_id,
  names_from = name,
  values_from = c(max, min, sum, mean, sd, cv),
  names_glue = "{name}_{.value}"
)

# Join raster stats
site_data = site_data %>% left_join(stats, by = "site_id")

# Visualize geospatial data ---------------------------------------------------------------------
message("Visualizing geospatial data")

# Order sites by development intensity
site_order = unique(site_data %>% arrange(desc(nlcd_developed_variable_intensity)) %>% pull(site_id))

nlcd_cols = grep("^nlcd_", names(site_data), value = TRUE)
nlcd_long = site_data %>% st_drop_geometry() %>% select(site_id, all_of(nlcd_cols)) %>%
  pivot_longer(
    cols = -site_id,
    names_to = "landcover",
    values_to = "percent"
  ) %>%
  mutate(
    landcover = gsub("^nlcd_", "", landcover),
    site_id = factor(site_id, levels = site_order)
  )

ggplot(nlcd_long, aes(y = site_id, x = percent, fill = landcover)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Percent cover",
    y = "Site ID",
    fill = "Summary class",
    title = paste0("NLCD land cover composition (", buffer_size, " m buffer)")
  ) +
  scale_fill_manual(values = c(
    "developed_variable_intensity"             = "#eb0000",
    "developed_open_space"             = "#eb9999",
    "forest"                = "#1c5f2c",
    "wetlands"              = "#6c9fb8",
    "agriculture"           = "#ead963",
    "barren_land"           = "#b3ac9f",
    "grassland_herbaceous"  = "#dde9af",
    "open_water"            = "#466b9f",
    "perennial_ice_snow"    = "#d1defa",
    "shrub_scrub"           = "#af963c"
  )) +
  theme_minimal()

# Visualize land cover for all sites

sites_lc_df = bind_rows(
  lapply(names(geospatial_site_data), function(site_id) {
    df = as.data.frame(geospatial_site_data[[site_id]]$rasters$class, xy = TRUE)
    value_col = names(df)[3]
    df = df %>% rename(value = all_of(value_col))  # rename for ggplot
    df$site <- site_id
    return(df)
  })
)
present_levels = nlcd_metadata %>% filter(class %in% unique(sites_lc_df$value))

# Order sites by decreasing amount of development
site_order = sites_lc_df %>% group_by(site) %>%
  summarise(developed_n = sum(grepl("^Developed", value)), .groups = "drop") %>%
  arrange(desc(developed_n)) %>% pull(site)
sites_lc_df = sites_lc_df %>%
  mutate(site = factor(site, levels = site_order))

ggplot(sites_lc_df, aes(x = x, y = y, fill = value)) + geom_tile() +
  scale_fill_manual(values = setNames(present_levels$color, present_levels$class), name = "Class") +
  facet_wrap(~site, scales = "free") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD land cover configuration (", buffer_size, " m buffer)"), x = "", y = "")

# ggplot expects data.frame objects for plotting
raster_to_df = function(r, site_name) {
  df = as.data.frame(r, xy = TRUE)
  colnames(df)[3] = "class"
  df$site = site_name
  return(df)
}

# Visualize imperviousness for all sites
sites_imp_df = bind_rows(
  lapply(seq_along(geospatial_site_data), function(i) {
    raster_to_df(geospatial_site_data[[i]]$rasters$rast_nlcd_impervious, names(geospatial_site_data)[i])
  })
) %>% mutate(site = factor(site, levels = site_order))

ggplot(sites_imp_df, aes(x = x, y = y, fill = class)) + geom_tile() +
  facet_wrap(~ site, scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD fractional impervious surface (", buffer_size, " m buffer)"), x = "", y = "", fill = "Impervious %")

# Visualize canopy cover for all sites
sites_tcc_df = bind_rows(
  lapply(seq_along(geospatial_site_data), function(i) {
    raster_to_df(geospatial_site_data[[i]]$rasters$rast_usfs_canopycover, names(geospatial_site_data)[i])
  })
) %>% mutate(site = factor(site, levels = site_order))

ggplot(sites_tcc_df, aes(x = x, y = y, fill = class)) + geom_tile() +
  facet_wrap(~ site, scales = "free") +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("USFS tree canopy cover (", buffer_size, " m buffer)"), x = "", y = "", fill = "Cover %")

## Visualize relationships between raster summary stats

# Site attributes (e.g. BIBI)
ggplot(site_data, aes(x = rast_nlcd_impervious_sum, y = bibi)) +
  geom_point() + geom_smooth(method = "lm") + geom_text_repel(aes(label = site_id)) 

# Comparing different rasters
ggplot(site_data, aes(x = rast_landfire_treeheight_sum, y = rast_gedi_height_sum)) +
  geom_point() + geom_abline(intercept = 0, slope = 1) +
  geom_smooth(method = "lm") + geom_text_repel(aes(label = site_id))

# Load species detection history data ------------------------------------------------
message("Loading species detection history data")

# Load species detection history data
detections = readRDS(in_cache_detections)
detections$long$common_name = tolower(detections$long$common_name)
detections$wide$common_name = tolower(detections$wide$common_name)

species_names = read_lines(in_path_species_list) %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>% mutate(
    common_name = tolower(common_name), scientific_name = tolower(scientific_name)
  ) %>% filter(common_name %in% sort(unique(detections$long$common_name)))

# Load AVONET species trait metadata
avonet = readxl::read_xlsx(in_path_avonet_traits, sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  mutate(scientific_name = tolower(scientific_name)) %>%
  filter(scientific_name %in% species_names$scientific_name) %>%
  select(scientific_name, family, order, mass, habitat, habitat_density, migration, trophic_level, trophic_niche, primary_lifestyle)

# Load EltonTraits 1.0 species trait metadata
eltontraits = read_csv(in_path_elton_traits) %>% clean_names() %>% rename(scientific_name = scientific, common_name = english) %>%
  mutate(
    common_name = str_to_lower(common_name), scientific_name = str_to_lower(scientific_name)
  ) %>% mutate(
    common_name = case_when(
      common_name == "american treecreeper" ~ "brown creeper",
      common_name == "winter wren"          ~ "pacific wren",
      common_name == "black-throated grey warbler" ~ "black-throated gray warbler",
      common_name == "common teal" ~ "green-winged teal",
      TRUE ~ common_name   # keep all others unchanged
    )
  ) %>% filter(common_name %in% species_names$common_name) %>%
  select(common_name, diet_inv, diet_5cat, starts_with("for_strat"), body_mass_value)

setdiff(species_names$common_name, eltontraits %>% filter(common_name %in% species_names$common_name) %>% pull(common_name))

# Add both trait data sources to species metadata
species_metadata = left_join(species_names, avonet, by = "scientific_name")
species_metadata = left_join(species_metadata, eltontraits, by = "common_name")

# NOTE: Very few aerial specialist foragers
eltontraits %>% filter(for_strat_aerial >= 10) %>% pull(common_name)

# Get AVONET invertivore community subset
# "Invertivore = species obtaining at least 60% of food resources from invertebrates in terrestrial systems, including insects, worms, arachnids, etc."
invertivores_avonet = species_metadata %>% filter(trophic_niche == "Invertivore") %>% pull(common_name)
# "Aquatic Predator = species obtaining at least 60% of food resources from vertebrate and invertebrate animals in aquatic systems, including fish, crustacea, molluscs, etc."
species_metadata %>% filter(trophic_niche == "Aquatic predator") %>% pull(common_name)

# Get Eltontraits invertivore community subset
# "Percent use of: Invertebrates-general, aquatic invertebrates, shrimp, krill, squid, crustacaeans, molluscs, cephalapod, polychaetes, gastropods, orthoptera, terrestrial Invertebrates, ground insects, insect larvae, worms, orthopterans, flying insects"
invertivores_eltontraits_dietGt10 = species_metadata %>% filter(diet_inv >= 10) %>% pull(common_name) # nearly all species
# "Assignment to the dominant among five diet categories based on the summed scores of constituent individual diets."
invertivores_eltontraits = species_metadata %>% filter(diet_5cat == "Invertebrate") %>% pull(common_name)

setdiff(invertivores_eltontraits, invertivores_avonet)
setdiff(invertivores_avonet, invertivores_eltontraits)

# Exclude certain sites from analysis ------------------------------------------------

# Exclude sites 257 and 259 that are dominated by agriculture
sites_to_exclude = c(257, 259)
message("Excluding agricultural site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data = site_data %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)
sites_to_exclude = setdiff(unique(site_data$site_id), unique(detections$long$site_id))
message("Excluding incomplete survey site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data = site_data %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

# Visualize joint data ----------------------------------------------------------------
message("Visualizing joint data")

presence_absence = detections$long %>% group_by(site_id, common_name) %>%
  summarise(presence = if_else(sum(n_detections, na.rm = TRUE) > 0, 1, 0), .groups = "drop")

# Calculate richness of different groups
richness = presence_absence %>% group_by(site_id) %>% summarise(richness = sum(presence))
richness_invert_et = presence_absence %>% filter(common_name %in% invertivores_eltontraits) %>%
  group_by(site_id) %>% summarise(richness_invert_et = sum(presence))
richness_invert_avo = presence_absence %>% filter(common_name %in% invertivores_avonet) %>%
  group_by(site_id) %>% summarise(richness_invert_avo = sum(presence))

# Join with site data
site_data = left_join(site_data, richness, by = "site_id")
site_data = left_join(site_data, richness_invert_et, by = "site_id")
site_data = left_join(site_data, richness_invert_avo, by = "site_id")

# Richness across AVONET guilds per site
richness_by_trophic_niche = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, trophic_niche) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_primary_lifestyle = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, primary_lifestyle) %>% summarise(count = sum(presence), .groups = "drop")

ggplot(richness_by_trophic_niche, aes(x = count, y = reorder(site_id, count), fill = trophic_niche)) + geom_col()
ggplot(richness_by_primary_lifestyle, aes(x = count, y = reorder(site_id, count), fill = primary_lifestyle)) + geom_col()

# Richness of guilds as a function of BIBI
ggplot(left_join(richness_by_trophic_niche, site_data, by = "site_id"),
       aes(x = bibi, y = count, color = trophic_niche, fill = trophic_niche)) +
  geom_point() + geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data, by = "site_id"),
       aes(x = bibi, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

# Wilson's Warbler presence/absence as a function of BIBI
ggplot(left_join(presence_absence %>% filter(common_name == "wilson's warbler"), site_data, by = "site_id"),
       aes(x = bibi, y = presence)) +
  geom_point() +geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

# Structural equation modeling -------------------------------------------------------

candidate_vars = c(
  # Predictors
  "bibi" = "bibi",
  "imperviousness_basin" = "basin_impervious",
  "imperviousness_local" = "rast_nlcd_impervious_mean",
  "abund_dev_varint" = "nlcd_developed_variable_intensity",
  "abund_dev_opensp" = "nlcd_developed_open_space",
  "abund_forest"     = "nlcd_forest",
  "abund_wetland"    = "nlcd_wetlands",
  "canopy_usfs"      = "rast_usfs_canopycover_sum",
  "canopy_gedi"      = "rast_gedi_cover_sum",
  "height_landfire"  = "rast_landfire_treeheight_mean",
  "height_gedi"      = "rast_gedi_height_mean",
  "height_cv_landfire"  = "rast_landfire_treeheight_cv",
  "height_cv_gedi"     = "rast_gedi_height_cv",
  "fhd"              = "rast_gedi_fhd_mean",
  # Responses
  "richness" = "richness",
  "richness_invert_et" = "richness_invert_et",
  "richness_invert_avo" = "richness_invert_avo"
)
data = as.data.frame(site_data) %>% select(all_of(candidate_vars))

# Explore pairwise collinearity among predictors
pairwise_collinearity = function(vars, threshold = 0.7) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

pairwise_collinearity(data %>% select(bibi, imperviousness_basin, canopy_usfs, height_landfire))

# Fit component regressions
# TODO: Include random intercepts in component models to account for nested structure of site within basin
m_bibi = lm(bibi ~ imperviousness_basin + canopy_usfs,
                data)

m_richness_invert = glm(richness_invert_avo ~ bibi + imperviousness_basin + canopy_usfs + height_landfire,
                        data, family = poisson(link = "log"))

# VIF analyses for multicollinearity
sort(car::vif(m_bibi))
sort(car::vif(m_richness_invert))

# Create structural equation model
sem_alpha_total = psem(
  m_bibi,
  m_richness_invert
)

# Inspect model structure
plot(sem_alpha_total)

# Inspect model claims and fit
# - A significant independence claim from a test of directed separation suggests that the path is missing
# or misspecified.
# - A significant global Fisher’s C p-value (< 0.05) suggests that the modeled structure is statistically
# significantly different than the structure implied by the data, and that alternative pathways or causal
# links with missing variables warrant further exploration
print(summary(sem_alpha_total))
