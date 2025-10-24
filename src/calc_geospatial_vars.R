############################################################################################################
# Load ARU and PSSB site locations and calculate all site covariates from land cover data
#
# Inputs:
in_cache_geospatial = "data/cache/preprocess_geospatial_data/geospatial_data.rds"
############################################################################################################

# Load required packages (install any missing)
pkgs = c(
  "tidyverse",  # data manipulation
  "sf",         # vector data
  "terra",      # raster data
  "tidyterra",  # raster data manipulation
  "tigris",     # political boundaries
  "geosphere",  # distance metrics
  "mapview",    # interactive geospatial visualization
  "progress"    # dynamic progress bar
)
sapply(pkgs, function(pkg) {
  if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
  as.character(packageVersion(pkg))
})

buffer_size = 500 # ~500m insect emergence flux range falloff
threshold = 0.9
days_threshold = 3

crs_standard = "EPSG:32610" # shared coordinate reference system (metric)

############################################################################################################
# Load ARU and PSSB site locations and define study area
message("Loading ARU and PSSB site locations")

site_metadata = read_csv("data/site_metadata.csv", show_col_types = FALSE)
site_metadata$year = year(site_metadata$date_start)

# Calculate distance between paired ARU and PSSB sites
site_metadata = site_metadata %>% rowwise() %>%
  mutate(dist_m = geosphere::distHaversine(
    c(long_aru,  lat_aru),
    c(long_pssb, lat_pssb)
  )) %>%
  ungroup()

sites_aru = site_metadata %>%
  st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326) %>% st_transform(crs = crs_standard)
sites_pssb = site_metadata %>%
  st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326) %>% st_transform(crs = crs_standard)

# Define study area by a buffered bounding box around the sites...
study_area = st_as_sfc(st_bbox(st_buffer(sites_aru, 1500)))
# ... or by a geopolitical boundary (e.g. King County)
study_area = counties(state = "WA", cb = TRUE) %>% filter(NAME == "King") %>% st_transform(crs = crs_standard)

mapview(study_area, alpha.regions = 0, lwd = 2) +
  mapview(sites_aru, zcol = "dist_m", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

############################################################################################################
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
# TODO: Load and process geospatial data
message("Loading geospatial data")

theme_set(theme_minimal())

# Load cached data
rast_filepaths = list.files(out_dir, pattern = "^rast_.*\\.tif$", full.names = TRUE)
rast_data = lapply(rast_filepaths, rast)
names(rast_data) = gsub("\\.tif$", "", basename(rast_filepaths))

sf_ripfb = st_read("data/cache/preprocess_geospatial_data/sf_ripfb.gpkg")

# Raster static plotting options
plot(rast_data$rast_landcover)
ggplot() + geom_spatraster(data = rast_data$rast_landcover)

plot(rast_data$rast_landcover)
plot(rast_data$rast_impervious)
plot(rast_data$rast_canopycover)
plot(rast_data$rast_vegcover)
plot(rast_data$rast_treeheight)
plot(rast_data$rast_shrubheight)
plot(rast_data$rast_herbheight)
plot(rast_data$rast_vegtype)
plot(rast_data$rast_sclass)
plot(rast_data$rast_fhd)
plot(rast_data$rast_cover)
plot(rast_data$rast_height)
plot(rast_data$rast_pavd55o10m)
plot(rast_data$rast_pavd20m)

# If your continuous raster object is very large
plot(rast_data$rast_canopycover) # plot statically
mapview(aggregate(rast_data$rast_canopycover, fact = 10, fun = "mean")) # ...plot dynamically at reduced (mean) resolution

# If your categorical raster object is very large
plot(rast_data$rast_landcover) # plot statically
mapview(aggregate(rast_data$rast_landcover, fact = 10, fun = "modal")) # ...plot dynamically at reduced (modal) resolution

# If your sf object is very large
ggplot(st_simplify(sf_ripfb, dTolerance = 250)) + geom_sf() # ...plot statically at reduced resolution
mapview(st_simplify(sf_ripfb, dTolerance = 250)) # ...plot dynamically at reduced resolution

############################################################
# Calculate land cover data for all sites
message("Calculating land cover composition with buffer size ", buffer_size)

nlcd_data = list()
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
for (s in 1:nrow(sites_aru)) {
  
  # TODO: Exclude specific sites from analysis?
  
  site = st_as_sf(sites_aru[s,])
  site_buffer = st_buffer(site, buffer_size)
  
  site_lc = mask(crop(rast_data$rast_landcover, site_buffer), site_buffer)
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
  
  site_imp = mask(crop(rast_data$rast_impervious, site_buffer), site_buffer)
  imp_sum  = sum(values(site_imp), na.rm = TRUE)
  
  site_tcc = mask(crop(rast_data$rast_canopycover, site_buffer), site_buffer)
  tcc_sum  = sum(values(site_tcc), na.rm = TRUE)
  
  site_treeheight = mask(crop(rast_data$rast_treeheight, site_buffer), site_buffer)
  treeheight_mean  = mean(values(site_treeheight), na.rm = TRUE)
  treeheight_sd   = sd(values(site_treeheight), na.rm = TRUE)
  
  # site_shrubheight = mask(crop(rast_shrubheight, site_buffer), site_buffer)
  # site_herbheight = mask(crop(rast_herbheight, site_buffer), site_buffer)
  
  # site_vegcover = mask(crop(rast_vegcover, site_buffer), site_buffer)
  
  site_ripfb = st_intersection(sf_ripfb, site_buffer)
  site_ripfb = st_union(site_ripfb)
  ripfb_area = st_area(site_ripfb)
  # mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_ripfb)
  
  nlcd_data[[as.character(site$site_id)]] = list(
    rast_lc       = site_lc,
    rast_impervious = site_imp,
    rast_tcc      = site_tcc,
    cover_summary = cover_summary,
    cover_detail  = cover_detail,
    imp_sum       = imp_sum,
    tcc_sum       = tcc_sum,
    treeheight_mean = treeheight_mean,
    treeheight_sd  = treeheight_sd,
    ripfb_area    = ripfb_area
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
  
  site_imp = nlcd_data[[1]]$rast_impervious
  site_imp_df = as.data.frame(site_imp, xy = TRUE)
  ggplot(site_imp_df, aes(x = x, y = y, fill = impervious)) + geom_raster() +
    scale_fill_continuous(limits = c(0, 100)) +
    coord_equal() + theme_minimal()
  
  mapview(site_lc) + mapview(site_imp)
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

# Extract sums for each site
imp_mean_df <- lapply(names(nlcd_data), function(id) {
  data.frame(
    site_id = id,
    imp_sum  = nlcd_data[[id]]$imp_sum,
    tcc_sum  = nlcd_data[[id]]$tcc_sum,
    treeheight_mean = nlcd_data[[id]]$treeheight_mean,
    treeheight_sd = nlcd_data[[id]]$treeheight_sd
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
    raster_to_df(nlcd_data[[i]]$rast_impervious, names(nlcd_data)[i])
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
  geom_point() + geom_smooth(method = "lm")

ggplot(site_data, aes(x = tcc_sum, y = bibi)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(site_data, aes(x = imp_sum, y = nlcd_developed_variable_intensity)) +
  geom_point()

ggplot(site_data, aes(x = imp_sum, y = tcc_sum)) +
  geom_point()