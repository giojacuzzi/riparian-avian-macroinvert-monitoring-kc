# Load ARU and PSSB site locations, define study area,
# and calculate all site covariates from land cover data

library(tidyverse)
library(terra)
library(mapview)
library(sf)
library(geosphere)
library(progress)

buffer_size = 500 # ~500m insect emergence flux range falloff; try 1km also
threshold = 0.9
days_threshold = 3

standard_crs_number = 32610
standard_crs_code = "EPSG:32610"

############################################################
# Load ARU and PSSB site locations and define study area
message("Loading ARU and PSSB site locations")

site_metadata = read_csv("data/site_metadata.csv", show_col_types = FALSE)
site_metadata$year = year(site_metadata$date_start)

# Calculate distance between paired ARU and PSSB sites
site_metadata = site_metadata %>% rowwise() %>%
  mutate(dist_m = geosphere::distHaversine(
    c(long_aru, lat_aru),  # c(lon, lat) for ARU
    c(long_pssb, lat_pssb) # c(lon, lat) for PSSB
  )) %>%
  ungroup()

sites_aru = site_metadata %>% st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326)
sites_pssb = site_metadata %>% st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326)
study_area = st_as_sfc(st_bbox(st_buffer(sites_aru, 10000)))

mapview(study_area, alpha.regions = 0, lwd = 2) +
  mapview(sites_aru, zcol = "dist_m", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

############################################################
# Load PSSB B-IBI data
message("Loading PSSB B-IBI data")

pssb_data = read_csv("data/raw/benthos/ScoresByYear.csv", show_col_types = FALSE) %>%
  janitor::clean_names() %>% select(site_id, x2020, x2021, x2022, x2023, x2024)

sites_aru = left_join(sites_aru, pssb_data, by = c("site_id"))
sites_aru = sites_aru %>% rowwise() %>% mutate(bibi_mean = mean(c_across(x2020:x2024), na.rm = TRUE)) %>% ungroup()
sites_aru = sites_aru %>% rowwise() %>% mutate(bibi = NA) %>%
  mutate(bibi = case_when( # Get the BIBI for the year in which the site was sampled
    year == 2024 ~ x2024,
    year == 2025 ~ x2024, # NOTE: 2025 not yet available
    TRUE ~ bibi
  ))

mapview(sites_aru, zcol = "bibi")

############################################################
# Load lc land cover and impervious surface data
message("Loading land cover and impervious surface data")

lc_raw  = rast("data/raw/environment/NLCD/Annual_NLCD_LndCov_2023_CU_C1V0.tif")
imp_raw = rast('data/raw/environment/NLCD/Annual_NLCD_FctImp_2023_CU_C1V0.tif')
stopifnot(crs(lc_raw) == crs(imp_raw))

study_area = project(vect(study_area), crs(lc_raw))
sites_aru  = project(vect(sites_aru), crs(lc_raw))
stopifnot(crs(lc_raw) == crs(study_area) && crs(lc_raw) == crs(sites_aru))

lc  = mask(crop(lc_raw, study_area), study_area)
imp = mask(crop(imp_raw, study_area), study_area)

# Factor land cover data
# https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
lc = as.factor(lc)
nlcd_levels <- data.frame(
  id = c(11,12,21,22,23,24,31,41,42,43,52,71,81,82,90,95),
  class = c(
    "Open Water", # areas of open water, generally with less than 25% cover of vegetation or soil.
    "Perennial Ice/Snow", # areas characterized by a perennial cover of ice and/or snow, generally greater than 25% of total cover.
    "Developed, Open Space", # areas with a mixture of some constructed materials, but mostly vegetation in the form of lawn grasses. Impervious surfaces account for less than 20% of total cover. These areas most commonly include large-lot single-family housing units, parks, golf courses, and vegetation planted in developed settings for recreation, erosion control, or aesthetic purposes.
    "Developed, Low Intensity", # areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 20% to 49% percent of total cover. These areas most commonly include single-family housing units.
    "Developed, Medium Intensity", # areas with a mixture of constructed materials and vegetation. Impervious surfaces account for 50% to 79% of the total cover. These areas most commonly include single-family housing units.
    "Developed, High Intensity", # highly developed areas where people reside or work in high numbers. Examples include apartment complexes, row houses and commercial/industrial. Impervious surfaces account for 80% to 100% of the total cover.
    "Barren Land", # areas of bedrock, desert pavement, scarps, talus, slides, volcanic material, glacial debris, sand dunes, strip mines, gravel pits and other accumulations of earthen material. Generally, vegetation accounts for less than 15% of total cover.
    "Deciduous Forest", # areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species shed foliage simultaneously in response to seasonal change.
    "Evergreen Forest", # areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. More than 75% of the tree species maintain their leaves all year. Canopy is never without green foliage.
    "Mixed Forest", # areas dominated by trees generally greater than 5 meters tall, and greater than 20% of total vegetation cover. Neither deciduous nor evergreen species are greater than 75% of total tree cover.
    "Shrub/Scrub", # areas dominated by shrubs; less than 5 meters tall with shrub canopy typically greater than 20% of total vegetation. This class includes true shrubs, young trees in an early successional stage or trees stunted from environmental conditions.
    "Grassland/Herbaceous", # areas dominated by gramanoid or herbaceous vegetation, generally greater than 80% of total vegetation. These areas are not subject to intensive management such as tilling, but can be utilized for grazing.
    "Pasture/Hay", # areas of grasses, legumes, or grass-legume mixtures planted for livestock grazing or the production of seed or hay crops, typically on a perennial cycle. Pasture/hay vegetation accounts for greater than 20% of total vegetation.
    "Cultivated Crops", # areas used for the production of annual crops, such as corn, soybeans, vegetables, tobacco, and cotton, and also perennial woody crops such as orchards and vineyards. Crop vegetation accounts for greater than 20% of total vegetation. This class also includes all land being actively tilled.
    "Woody Wetlands", # areas where forest or shrubland vegetation accounts for greater than 20% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
    "Emergent Herbaceous Wetlands" # areas where perennial herbaceous vegetation accounts for greater than 80% of vegetative cover and the soil or substrate is periodically saturated with or covered with water.
  ),
  color = c(
    "#466b9f",
    "#d1defa",
    "#dec5c5",
    "#d99282",
    "#eb0000",
    "#ab0000",
    "#b3ac9f",
    "#68ab5f",
    "#1c5f2c",
    "#b5c58f",
    "#af963c",
    "#dde9af",
    "#ead963",
    "#ab6c28",
    "#b3d1d1",
    "#6c9fb8"
  )
)
levels(lc) <- nlcd_levels
levels(lc)[[1]] # Verify

plot(lc, col = nlcd_levels$color)
plot(imp)

# Project all data to EPSG:32610
sites_aru  = project(sites_aru, standard_crs_code)
study_area = project(study_area, standard_crs_code)
lc  = project(lc, standard_crs_code)
imp = project(imp, standard_crs_code)

############################################################
# Load tree canopy cover
message("Loading tree canopy cover data")
tcc_raw = rast("data/raw/environment/Forest Service Science TCC/science_tcc_conus_wgs84_v2023-5_20230101_20231231.tif")

template = project(study_area, crs(tcc_raw))
message("Cropping tree canopy cover data")
tcc = mask(crop(tcc_raw, template), template)
message("Projecting tree canopy cover data")
tcc = project(tcc, standard_crs_code)
# Clean data
tcc[tcc > 100] = NA
tcc[is.nan(tcc)] <- NA

############################################################
# Calculate land cover data for all sites
message("Calculating land cover composition with buffer size ", buffer_size)

nlcd_data = list()
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
for (s in 1:nrow(sites_aru)) {
  
  # TODO: Exclude specific sites from analysis?
  
  site = st_as_sf(sites_aru[s])
  site_buffer = st_buffer(site, buffer_size)

  site_lc = mask(crop(lc, site_buffer), site_buffer)
  site_lc[] = as.integer(factor(site_lc[], levels = nlcd_levels$id))
  levels(site_lc) = data.frame(id = 1:nrow(nlcd_levels), class = nlcd_levels$class)
  # mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_lc) + mapview(site_imp)
  
  freq_table <- freq(site_lc)
  total_cells <- sum(freq_table$count)
  freq_percent <- data.frame(
    class = freq_table$value,
    count = freq_table$count,
    percent = (freq_table$count / total_cells) * 100
  )
  cover_detail = full_join(freq_percent, nlcd_levels %>% select(class), by = "class")
  cover_detail = cover_detail %>% mutate(
    count = ifelse(is.na(count), 0, count),
    percent = ifelse(is.na(percent), 0, percent)
  ) %>% arrange(desc(percent))
  # print(cover_detail)
  
  cover_summary <- cover_detail %>%
    mutate(
      group = case_when(
        grepl("Intensity", class) ~ "Developed, Variable Intensity",
        grepl("Forest", class) ~ "Forest",
        # grepl("Shrub|Grass", class) ~ "Grass/Shrub",
        grepl("Wetlands", class) ~ "Wetlands",
        grepl("Pasture|Cultivated", class) ~ "Agriculture",
        TRUE ~ class  # keep other classes as-is
      )
    ) %>%
    group_by(group) %>%
    summarise(
      count = sum(count),
      percent = sum(percent)
    ) %>%
    arrange(desc(percent))
  # print(cover_summary)
  
  present_ids = sort(unique(site_lc[]))
  present_levels = nlcd_levels[present_ids, ]
  # plot(site_lc, col = present_levels$color)
  
  site_imp = mask(crop(imp, site_buffer), site_buffer)
  imp_mean = mean(values(site_imp), na.rm = TRUE)
  imp_sum  = sum(values(site_imp), na.rm = TRUE)
  
  site_tcc = mask(crop(tcc, site_buffer), site_buffer)
  tcc_mean = mean(values(site_tcc), na.rm = TRUE)
  tcc_sum  = sum(values(site_tcc), na.rm = TRUE)
  
  nlcd_data[[as.character(site$site_id)]] = list(
    rast_lc       = site_lc,
    rast_imp      = site_imp,
    rast_tcc      = site_tcc,
    cover_summary = cover_summary,
    cover_detail  = cover_detail,
    imp_mean      = imp_mean,
    imp_sum       = imp_sum,
    tcc_mean      = tcc_mean,
    tcc_sum       = tcc_sum
  )
  pb$tick()
}

# Inspect a particular site
if (FALSE) {
  site_id = "189"
  site_lc = nlcd_data[[site_id]]$rast_lc
  present_ids = sort(unique(site_lc[]))
  present_levels = nlcd_levels[present_ids, ]
  site_lc_df = as.data.frame(site_lc, xy = TRUE) %>% left_join(nlcd_levels, by = "class")
  ggplot(site_lc_df, aes(x = x, y = y, fill = class)) + geom_raster() +
    scale_fill_manual(values = setNames(present_levels$color, present_levels$class)) +
    coord_equal() + theme_minimal()
  
  site_imp = nlcd_data[[1]]$rast_imp
  site_imp_df = as.data.frame(site_imp, xy = TRUE)
  ggplot(site_imp_df, aes(x = x, y = y, fill = Annual_NLCD_FctImp_2023_CU_C1V0)) + geom_raster() +
    scale_fill_continuous(limits = c(0, 100)) +
    coord_equal() + theme_minimal()
}

############################################################
# Join land cover data with sites into a single table
message("Joining data")

# Convert nlcd_data into a single table and pivot wider
nlcd_summary <- lapply(names(nlcd_data), function(id) {
  df <- nlcd_data[[id]]$cover_summary
  df$site_id <- id
  df
}) %>% bind_rows()
nlcd_summary <- nlcd_summary %>%
  select(site_id, group, percent) %>%
  pivot_wider(names_from = group, values_from = percent, names_prefix = "nlcd_") %>%
  janitor::clean_names()

# Join land cover data
sites_df <- as.data.frame(sites_aru)  # extract attributes
sites_df$site_id = as.character(sites_df$site_id)
site_data <- sites_df %>% left_join(nlcd_summary, by = "site_id")

# Extract imp_mean for each site
imp_mean_df <- lapply(names(nlcd_data), function(id) {
  data.frame(
    site_id = id,
    imp_mean = nlcd_data[[id]]$imp_mean,
    imp_sum  = nlcd_data[[id]]$imp_sum,
    tcc_mean = nlcd_data[[id]]$tcc_mean,
    tcc_sum  = nlcd_data[[id]]$tcc_sum
  )
}) %>% bind_rows()

# TODO: for other imp_sum etc.

# Join impervious surface
site_data <- site_data %>%
  left_join(imp_mean_df, by = "site_id")

############################################################
# Visualize data
message("Visualizing data")

# Step 1: Order site_id by nlcd_Developed
site_order <- site_data %>% arrange(desc(nlcd_developed_variable_intensity)) %>% pull(site_id)

# Step 2: Convert to long format (same as before)
nlcd_cols <- grep("^nlcd_", names(site_data), value = TRUE)

nlcd_long <- site_data %>% select(site_id, all_of(nlcd_cols)) %>%
  pivot_longer(
    cols = -site_id,
    names_to = "landcover",
    values_to = "percent"
  ) %>%
  mutate(
    landcover = gsub("^nlcd_", "", landcover),
    site_id = factor(site_id, levels = site_order)  # reorder sites
  )

# Step 3: Plot
p = ggplot(nlcd_long, aes(y = site_id, x = percent, fill = landcover)) +
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
  theme_minimal(); print(p)

# Visualize land cover rasters for all sites

raster_to_df <- function(r, site_name) {
  df <- as.data.frame(r, xy = TRUE)
  colnames(df)[3] <- "class"
  df$site <- site_name
  return(df)
}

all_sites_df <- bind_rows(
  lapply(seq_along(nlcd_data), function(i) {
    raster_to_df(nlcd_data[[i]]$rast_lc, names(nlcd_data)[i])
  })
)
all_sites_df$class <- factor(all_sites_df$class, levels = nlcd_levels$class)
present_levels <- nlcd_levels %>% filter(class %in% unique(all_sites_df$class))

# Order sites by decreasing amount of development
site_order = all_sites_df %>%
  group_by(site) %>%
  summarise(developed_n = sum(grepl("^Developed", class)), .groups = "drop") %>%
  arrange(desc(developed_n)) %>% pull(site)
all_sites_df = all_sites_df %>%
  mutate(site = factor(site, levels = site_order))

# Plot
p = ggplot(all_sites_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  scale_fill_manual(values = setNames(present_levels$color, present_levels$class), name = "Class") +
  facet_wrap(~ site, scales = "free") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD land cover configuration (", buffer_size, " m buffer)"), x = "", y = ""); print(p)

imp_df <- bind_rows(
  lapply(seq_along(nlcd_data), function(i) {
    raster_to_df(nlcd_data[[i]]$rast_imp, names(nlcd_data)[i])
  })
) %>%
  mutate(site = factor(site, levels = site_order))

p = ggplot(imp_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  facet_wrap(~ site, scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD fractional impervious surface (", buffer_size, " m buffer)"), x = "", y = "", fill = "Impervious %"); print(p)

tcc_df <- bind_rows(
  lapply(seq_along(nlcd_data), function(i) {
    raster_to_df(nlcd_data[[i]]$rast_tcc, names(nlcd_data)[i])
  })
) %>% mutate(class = as.numeric(class)) %>%
  mutate(site = factor(site, levels = site_order))

p = ggplot(tcc_df, aes(x = x, y = y, fill = class)) +
  geom_tile() +
  facet_wrap(~ site, scales = "free") +
  # scale_fill_viridis_c(option = "D") +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("USFS tree canopy cover (", buffer_size, " m buffer)"), x = "", y = "", fill = "Cover %"); print(p)

ggplot(site_data, aes(x = imp_sum, y = bibi)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_minimal()

ggplot(site_data, aes(x = imp_mean, y = nlcd_developed_variable_intensity)) +
  geom_point() + theme_minimal()

ggplot(site_data, aes(x = tcc_sum, y = bibi)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_minimal()

############################################################
# Derive putative detection histories and diversity metrics

# Load classifier prediction data
message("Loading classifier prediction data")
path_prediction_data = "data/cache/aggregate_pam_data/prediction_data.feather"
prediction_data = arrow::read_feather(path_prediction_data) %>% rename(site_id = site)

# Discard sites 257 and 259 that are dominated by agriculture
sites_to_discard = as.character(c(257, 259))
message("Discarding sites ", paste(sites_to_discard, collapse = ", "))
prediction_data = prediction_data %>% filter(!(site_id %in% sites_to_discard))

# Obtain putative detections with naive threshold
message("Obtaining putative detections with naive threshold")
detections = prediction_data %>% filter(confidence >= threshold)

start_date_2024 <- min(detections$date[detections$season == "2024"])
end_date_2024 <- max(detections$date[detections$season == "2024"])
start_date_2025 <- min(detections$date[detections$season == "2025"])
end_date_2025 <- max(detections$date[detections$season == "2025"])

species = read_lines("data/audio/Full_Species_List.txt") %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>%
  filter(common_name %in% sort(unique(detections$common_name)))

# Break out insectivore guild into aerial insectivores via AVONET
avonet = readxl::read_xlsx('data/traits/AVONET Supplementary dataset 1.xlsx', sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  filter(scientific_name %in% species$scientific_name) %>%
  select(scientific_name, family, order, mass, habitat, habitat_density, migration, trophic_level, trophic_niche, primary_lifestyle)
species_metadata = left_join(species, avonet, by = "scientific_name")

insectivore_species = species_metadata %>%
  filter(trophic_niche %in% c("Invertivore"))

# insectivore_species = species_metadata %>%
#   filter(trophic_niche %in% c("Invertivore", "Aquatic predator")) %>%
#   filter(primary_lifestyle %in% c("Aerial", "Insessorial")) %>%
#   pull(common_name)

alt_insectivores = read_csv("data/processed/Species_Habitat_List.csv", show_col_types = FALSE) %>% rename(common_name = species) %>% filter(common_name %in% species$common_name) %>%
  filter(insectivore == "Yes") %>%
  pull(common_name)

# site_species_matrix_surveys <- detections %>%
#   distinct(site_id, common_name) %>%
#   mutate(present = 1) %>%
#   pivot_wider(
#     names_from = common_name,
#     values_from = present,
#     values_fill = 0
#   )
# site_species_matrix_surveys

detections_with_survey <- detections %>%
  group_by(site_id) %>%
  mutate(survey_num = as.integer(difftime(date, min(date), units = "days")) + 1) %>%
  ungroup()

detections_presence <- detections_with_survey %>%
  group_by(site_id, survey_num, common_name) %>%
  summarise(present = 1, .groups = "drop")  # 1 = detected at least once that day

species_site_survey_matrix <- detections_presence %>%
  pivot_wider(
    names_from = survey_num,
    values_from = present,
    values_fill = 0  # 0 = not detected that day
  )

species_matrices <- species_site_survey_matrix %>%
  group_split(common_name) %>%
  setNames(unique(detections_presence$common_name))

species_site_survey_long <- species_site_survey_matrix %>%
  pivot_longer(
    cols = -c(site_id, common_name),
    names_to = "survey_num",
    values_to = "detected"
  ) %>%
  mutate(survey_num = as.integer(survey_num))

detections_summary <- species_site_survey_long %>%
  group_by(common_name, survey_num) %>%
  summarise(total_detections = sum(detected), .groups = "drop")

ggplot(detections_summary, aes(x = survey_num, y = as.integer(total_detections))) +
  geom_line() +
  facet_wrap(~ common_name, scales = "free_y") +
  theme_minimal() +
  labs(
    title = "Species detections across all sites and surveys",
    x = "Survey number",
    y = "Total detections across sites"
  )

# TODO: Species accumulation curves?

site_species_matrix_days_detected <- detections %>%
  group_by(site_id, common_name) %>%
  summarise(n_dates = n_distinct(date), .groups = "drop") %>%  # count unique survey dates
  pivot_wider(
    names_from = common_name,
    values_from = n_dates,
    values_fill = 0   # species not detected on any date = 0
  )
site_species_matrix <- site_species_matrix_days_detected %>%
  mutate(across(-site_id, ~ if_else(. > days_threshold, 1, 0)))

site_insectivores_matrix <- site_species_matrix %>%
  select(site_id, all_of(insectivore_species$common_name))

# Visualize species composition at each site
site_species_long <- site_species_matrix %>%
  pivot_longer(cols = -site_id, names_to = "common_name", values_to = "detected") %>%
  filter(detected == 1) %>%
  group_by(site_id) %>%
  mutate(richness = n()) %>%
  ungroup() %>%
  arrange(richness)

ggplot(left_join(site_species_long, species_metadata, by = "common_name"), aes(y = reorder(site_id, richness), fill = trophic_niche)) +
  geom_bar() +
  labs(title = "Species presence by trophic niche") +
  theme_minimal()

ggplot(left_join(site_species_long, species_metadata, by = "common_name"), aes(y = reorder(site_id, richness), fill = primary_lifestyle)) +
  geom_bar() +
  labs(title = "Species presence by primary foraging lifestyle") +
  theme_minimal()

ggplot(left_join(site_species_matrix, site_data, by = "site_id"),
       aes(x = bibi, y = `Wilson's Warbler`)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  geom_point() + theme_minimal()

ggplot(left_join(site_species_matrix, site_data, by = "site_id"),
       aes(x = bibi, y = `Belted Kingfisher`)) +
  geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) +
  geom_point() + theme_minimal()
  
# Calculate richness of different groups

richness = site_species_matrix %>% mutate(richness = rowSums(across(-site_id))) %>% select(site_id, richness)
richness_insectivore = site_insectivores_matrix %>% mutate(richness_insectivore = rowSums(across(-site_id))) %>% select(site_id, richness_insectivore)

# Join with other site data
site_data_bird = right_join(site_data, richness, by = "site_id")
site_data_bird = right_join(site_data_bird, richness_insectivore, by = "site_id")

ggplot(site_data_bird, aes(x = bibi, y = richness)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_minimal()

ggplot(site_data_bird, aes(x = bibi, y = richness_insectivore)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_minimal()

ggplot(site_data_bird, aes(x = tcc_sum, y = richness_insectivore)) +
  geom_point() + geom_smooth(method = "lm") +
  theme_minimal()

site_data_bird = left_join(site_data_bird %>% mutate(site_id = as.double(site_id)), site_metadata %>% select(site_id, lat_aru, long_aru), by = "site_id")
site_data_bird = site_data_bird %>% st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326)
mapview(site_data_bird, zcol = "richness_insectivore")

##############

# library(vegan)
# 
# species_mat <- site_species_matrix %>%
#   column_to_rownames("site_id") %>%
#   as.matrix()
# 
# nmds <- metaMDS(species_mat, distance = "jaccard", binary = TRUE, k = 2, trymax = 100)
# nmds
# 
# nmds_sites <- as.data.frame(scores(nmds, display = "sites")) %>%
#   rownames_to_column("site_id")
# 
# env_data <- site_data_bird %>% select(site_id, bibi, tcc_sum, imp_sum) %>% mutate(site_id = as.character(site_id))
# 
# nmds_sites <- nmds_sites %>%
#   left_join(env_data, by = "site_id")
# 
# nmds_species <- as.data.frame(scores(nmds, display = "species")) %>%
#   rownames_to_column("species")
# 
# ggplot(nmds_sites, aes(x = NMDS1, y = NMDS2)) +
#   geom_point(size = 3) +   # color by environmental variable
#   stat_ellipse(level = 0.95) +
#   geom_text(aes(label = site_id), vjust = -1, size = 3) +
#   scale_color_viridis_c(option = "plasma") +
#   theme_minimal(base_size = 12) +
#   labs(title = "NMDS of species composition")

############################################################
# SEM

library(piecewiseSEM)

ggplot(site_data_bird %>% st_drop_geometry() %>% select(imp_sum, tcc_sum) %>% pivot_longer(cols = everything(), names_to = "variable", values_to = "value"), aes(x = value, fill = variable)) +
  geom_histogram(show.legend = FALSE) +
  facet_wrap(~ variable, scales = "free") +
  theme_minimal()

candidate_vars = c("bibi",
                   "imp_sum", "nlcd_developed_variable_intensity", "nlcd_developed_open_space",
                   "richness", "richness_insectivore",
                   "nlcd_forest", "tcc_sum",
                   "nlcd_shrub_scrub", "nlcd_wetlands")

data = as.data.frame(site_data_bird) %>% janitor::clean_names() %>% select(
  all_of(candidate_vars)
)

# Explore pairwise collinearity
pairwise_collinearity = function(vars, threshold = 0.8) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

pairwise_collinearity(data)
data_subset = data %>% select(bibi, imp_sum, tcc_sum)
pairwise_collinearity(data_subset)

# VIF analysis for multicollinearity (consider dropping variable(s) with high VIF values (> 10))
model = lm(rnorm(nrow(data_subset)) ~ ., data = data_subset)
sort(car::vif(model))

# Fit component regressions
model_bibi = lm(
  bibi ~ imp_sum + tcc_sum,
  data
)
model_alpha_total = glm(
  richness_insectivore ~ bibi + imp_sum + tcc_sum,
  data,
  family = poisson(link = "log")
)

# Create structural equation model
sem_alpha_total = psem(
  model_bibi,
  model_alpha_total
)

# Inspect model structure
sem_alpha_total
plot(sem_alpha_total)

# Conduct tests of directed separation (for each missing path)
# Establish the basis set & evaluate independence claims
# Use `dsep` function to perform the tests automagically dSep(sem_alpha_total)
# Use `fisherC` function to evaluate claims
# A significant global Fisher’s C p-value (< 0.05) suggests that the modeled structure is statistically significantly different than the structure implied by the data, and that alternative pathways or causal links with missing variables warrant further exploration
# Nagelkerke R2 describes proportion of variance explained by the model
print(summary(sem_alpha_total))
