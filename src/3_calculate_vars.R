# 3_calculate_vars.R ===================================================================
# Load ARU and PSSB site locations, calculate all site covariates from land cover data
#
# Inputs:
overwrite_geospatial_site_data = FALSE
in_path_pssb_data       = "data/raw/pssb/ScoresByYear.csv"
in_path_nlcd_metadata   = "data/raw/nlcd_metadata.csv"
in_cache_geospatial_dir = "data/cache/2_preprocess_geospatial_data"
# Outputs:
out_cache_dir = "data/cache/3_calculate_vars"

source("src/global.R")

# site buffer:
# ~550 m insect emergence 90% flux range falloff
# 5000 m (5 km) basin/watershed
buffer_radius = set_units(550, m)

# Load ARU and PSSB site locations and define study area -----------------------------
message("Loading ARU and PSSB site locations")

# Retrieve study area
study_area = st_read(paste0(in_cache_geospatial_dir, "/sf_studyarea.gpkg"), quiet = TRUE)

# View study area and sites
mapview(study_area, alpha.regions = 0, lwd = 2, layer.name = "King County") +
  mapview(sites_aru, col.region = "green", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

# Visualize alternative buffers
mapview(sites_aru) +
  mapview(st_buffer(sites_aru, 250), alpha.regions = 0, lwd = 3, color = "black", col.regions = "black", layer.name = "250 m - homerange ballpark") +
  mapview(st_buffer(sites_aru, 550), alpha.regions = 0, lwd = 3, color = "green", col.regions = "green", layer.name = "550 m - 90% dispersal") +
  mapview(st_buffer(sites_aru, 1000), alpha.regions = 0, lwd = 3, color = "blue", col.regions = "blue", layer.name = "1 km") +
  mapview(st_buffer(sites_aru, 5000), alpha.regions = 0, lwd = 3, color = "red", col.regions = "red", layer.name = "5 km")

# Calculate distance between paired ARU and PSSB sites
summary(site_metadata %>% rowwise() %>%
          mutate(dist_m = geosphere::distHaversine(c(long_aru,  lat_aru), c(long_pssb, lat_pssb))) %>% pull(dist_m))

# Load Puget Sound Stream Benthos monitoring data ------------------------------------
message("Loading PSSB BIBI data")

pssb_data_by_year = read_csv(in_path_pssb_data, show_col_types = FALSE) %>% clean_names() %>%
  mutate(site_id = as.character(site_id)) %>% select(site_id, x2024)

# Get the BIBI for the year in which the site was sampled
sites_aru = left_join(sites_aru, pssb_data_by_year, by = c("site_id"))
sites_aru = sites_aru %>% mutate(year = year(date_start)) %>% rowwise() %>% mutate(bibi = NA) %>%
  mutate(bibi = case_when(
    year == 2024 ~ x2024,
    year == 2025 ~ x2024, # TODO: 2025 data not yet available
    TRUE ~ bibi
  ))

# View sites by BIBI value
mapview(sites_aru, zcol = "bibi")

# Classify values into discrete categories
m = mapview(sites_aru, zcol = "bibi", at = c(0, 20, 40, 60, 80, 100),
        col.regions = c("red", "orange", "yellow", "green", "blue")); m

# Compare B-IBI to EPT Richness
pssb_data_scores_all = read_tsv("data/raw/pssb/Scores.txt") %>% clean_names() %>% mutate(site_id = as.character(site_id))
message("Pearson's R (all): ", round(cor(pssb_data_scores_all$ept_richness_quantity, pssb_data_scores_all$overall_score), 3))
ggplot(pssb_data_scores_all, aes(x = ept_richness_quantity, y = overall_score)) +
  geom_point()
pssb_data_scores_sites = pssb_data_scores_all %>% filter(site_id %in% sites_aru$site_id)
message("Pearson's R (sites): ", round(cor(pssb_data_scores_sites$ept_richness_quantity, pssb_data_scores_sites$overall_score), 3))
ggplot(pssb_data_scores_sites, aes(x = ept_richness_quantity, y = overall_score)) +
  geom_point()

cor.test(pssb_data_scores_sites$ept_richness_quantity,
         pssb_data_scores_sites$overall_score,
         method = "pearson")

# Load and process geospatial data ---------------------------------------------------
message("Loading geospatial data")

# Load ICLUS impervious surface projections
message("- ICLUS projection rasters")
iclus_rast_filepaths = list.files(file.path(in_cache_geospatial_dir, "iclus"), pattern = "^is.*\\.tif$", full.names = TRUE, recursive = TRUE)

rast_imp_iclus_projections = list()
for (f in iclus_rast_filepaths) {
  filename = basename(f)
  year = str_extract(filename, "(?<=is)\\d{4}")
  scenario = str_extract(filename, "(?<=\\d{4})[a-z0-9]+")
  r = rast(f)
  if (is.null(rast_imp_iclus_projections[[scenario]])) {
    rast_imp_iclus_projections[[scenario]] = list()
  }
  rast_imp_iclus_projections[[scenario]][[year]] = r
}

# Load cached raster data
message("- Rasters")
rast_filepaths = list.files(in_cache_geospatial_dir, pattern = "^rast_.*\\.tif$", full.names = TRUE)
rast_data = lapply(rast_filepaths, rast)
names(rast_data) = gsub("\\.tif$", "", basename(rast_filepaths))

# Load NLCD metadata for cover class codes
nlcd_metadata = read.csv(in_path_nlcd_metadata)

# Load cached vector data
message("- King County riparian zones")
sf_ripfb = st_read(paste0(in_cache_geospatial_dir, "/sf_ripfb.gpkg"), quiet = TRUE)

message("- NHD flowlines")
sf_flowline = st_read(paste0(in_cache_geospatial_dir, "/sf_flowline.gpkg"), quiet = TRUE) %>% clean_names() %>%
  filter(f_code %in% c(
    46000, 46003, 46006, 46007, # include only stream/river flowlines
    55800 # NOTE: Sammamish River is designated as arificial (55800) but is covered by USFS riparian area
))
sf_flowline$f_code = factor(sf_flowline$f_code)

message("- NHD waterbodies")
sf_waterbody = st_read(paste0(in_cache_geospatial_dir, "/sf_waterbody.gpkg"), quiet = TRUE) %>% clean_names()

# Load cached basin sf objects and retain only those sampled
message("- HUC 12 basins")
sf_basins12d = st_read(paste0(in_cache_geospatial_dir, "/sf_basins12d.gpkg"), quiet = TRUE) %>%
  clean_names() %>% select(huc12, name, area_sq_km) %>% mutate(basin_name = name, basin_area = area_sq_km)
sf_basins12d = sf_basins12d %>% filter(lengths(st_intersects(., st_transform(sites_aru, st_crs(sf_basins12d)))) > 0) # TODO: exclude lakes?

# Load road sf data
message("- Roads")
roads = st_read("data/raw/tl_2024_53033_roads/tl_2024_53033_roads.shp", quiet = TRUE) %>% clean_names() %>% select(rttyp, mtfcc)
roads_paved = roads %>% filter(mtfcc %in% c("S1100", "S1200", "S1400", "S1630", "S1640", "S1740", "S1780"))
roads_major = roads %>% filter(mtfcc %in% c("S1100", "S1200"))

# Store subsequent data calculated per site in `site_data`
site_data = sites_aru

# Calculate elevation
# site_data = elevatr::get_elev_point(site_data)

# Calculate imperviousness at the basin (landscape) scale ------------------------------------
message("Calculating imperviousness at the basin scale")

# Calculate mean imperviousness per basin
basin_impervious = terra::extract(rast_data$rast_nlcd_impervious, vect(sf_basins12d), fun = mean, na.rm = TRUE)
sf_basins12d = sf_basins12d %>% mutate(basin_impervious = basin_impervious[[2]])

# Store results
site_data = st_intersection(st_transform(site_data, st_crs(sf_basins12d)), sf_basins12d)

# basin_canopy = terra::extract(rast_data$rast_usfs_canopycover, vect(sf_basins12d), fun = mean, na.rm = TRUE)
# sf_basins12d = sf_basins12d %>% mutate(basin_canopy = basin_canopy[[2]])

# Calculate geospatial variables for all sites --------------------------------
message("Calculating geospatial variables for all sites (site buffer size ", buffer_radius, " m)")

# Helper function to project a buffer to a raster's
# crs and crop and mask the raster to that buffer
transform_crop_and_mask = function(rast, buff_crop, buff_mask) {
  buff_crop = st_transform(buff_crop, st_crs(rast))
  buff_crop = vect(buff_crop)
  buff_mask = st_transform(buff_mask, st_crs(rast))
  buff_mask = vect(buff_mask)
  mask(crop(rast, buff_crop), buff_mask)
}

# Helper function to calculate summary statistics for a continuous raster
rast_stats = function(rast, na.rm = TRUE) {
  stats = global(rast, fun = c("max", "min", "sum", "mean", "sd"), na.rm = na.rm)
  stats$cv = stats$sd / stats$mean
  stats$sum_proportion = stats$sum / (global(!is.na(rast), "sum", na.rm=TRUE)[[1]] * 100)
  return(stats)
}

# TODO: SITE IDXs TO REVIEW
# Questionable: 37, 38, 45, 47
# Good: 39, 220, 50

site_zones = list()
## Compute site buffer and interaction zones
message("Computing site buffers and interaction zones")
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
for (s in 1:nrow(sites_aru)) {

  s_id = sites_aru$site_id[s]
  # print(s)
  # print(s_id)
  
  # Create a buffer around the ARU site
  site = site_metadata %>% filter(site_id == s_id) %>% st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326) #st_as_sf(sites_aru[s,])
  # site_buffer = st_buffer(site, buffer_radius)
  site_buffer_proj = st_buffer(st_transform(site, crs_standard), dist = buffer_radius)
  site_buffer = st_transform(site_buffer_proj, st_crs(site))
  
  site_pssb = site_metadata %>% filter(site_id == s_id) %>% st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326)
  # mapview(site) + mapview(site_pssb) + mapview(site_buffer) + mapview(site_buffer_proj)
  
  ## Buffers
  
  # 1. Get flowline and waterbodies within buffer of site
  flowlines = st_crop(sf_flowline, st_transform(site_buffer, crs = st_crs(sf_flowline))) %>% st_transform(crs_standard)
  flowlines = st_intersection(st_transform(flowlines, st_crs(site_buffer)), site_buffer)
  
  # mapview(flowlines, zcol = "f_code")
  
  f_code = flowlines$f_code[st_nearest_feature(site, flowlines)]
  
  # NOTE: Intermittent streams (46003) are sampled
  
  # Buffer waterbodies by 100m
  waterbodies = st_crop(sf_waterbody, st_transform(site_buffer, crs = st_crs(sf_waterbody))) %>% st_transform(crs_standard)
  waterbodies = st_intersection(st_transform(waterbodies, st_crs(site_buffer)), site_buffer) %>% st_make_valid()
  waterbodies_buffer = st_union(st_buffer(st_transform(waterbodies, crs_standard), 100))
  waterbodies_buffer = st_difference(waterbodies_buffer, st_union(st_transform(waterbodies, crs_standard)))
  # Remove waterbodies (without buffer) from flowlines
  flowlines_trimmed = st_difference(flowlines, st_union(waterbodies)) %>%             filter(st_geometry_type(.) %in% c("LINESTRING", "MULTILINESTRING"))

  # mapview(flowlines_trimmed) +
  #   mapview(flowlines, zcol = "f_code") +
  #   mapview(waterbodies, col.regions = "red")
  
  flowlines = flowlines_trimmed
  
  # 2. Buffer 50 m core emergence zone minimum from flowlines, and buffer 100m to include adjacent edge habitats
  emergence_zone = st_make_valid(st_union(st_buffer(flowlines %>% st_transform(crs_standard), 100)))
  
  # Remove waterbodies from emergence zone
  emergence_zone = st_difference(emergence_zone, st_union(waterbodies) %>% st_transform(st_crs(emergence_zone)))
  
  interaction_zone = emergence_zone
  
  # 3. Merge with USFS riparian areas (if they exist) buffered by 100 m to include adjacent edge habitats
  site_riparian_usfs = transform_crop_and_mask(rast_data$rast_riparian, site_buffer, site_buffer)
  riparian_zone = as.polygons(site_riparian_usfs == 1, dissolve = TRUE)
  
  if (nrow(riparian_zone) > 0) {
    area_usfs_riparian_zone = st_area(st_as_sf(riparian_zone))
    riparian_zone_buffered = st_buffer(st_as_sf(riparian_zone), 100)
    # Remove waterbodies from riparian zone
    riparian_zone_no_waterbodies = st_difference(riparian_zone_buffered, st_union(waterbodies) %>% st_transform(st_crs(riparian_zone)))
    # Merge with interaction zone
    interaction_zone = st_union(interaction_zone, st_transform(riparian_zone_no_waterbodies, crs_standard))
  } else {
    area_usfs_riparian_zone = 0
  }
  
  interaction_zone = st_union(interaction_zone, waterbodies_buffer)
  
  # EDIT: Buffer the flowlines by 550, not 50 or 100
  # buffy = st_make_valid(st_union(st_buffer(flowlines %>% st_transform(crs_standard), 550)))
  # interaction_zone = st_union(interaction_zone, buffy)
  # interaction_zone = st_difference(interaction_zone, st_union(waterbodies) %>% st_transform(st_crs(interaction_zone)))
  
  # Only select nearest interaction zone
  interaction_zone = st_cast(interaction_zone, "POLYGON")
  nearest_interaction_zone = st_nearest_feature(st_transform(site, crs_standard), interaction_zone)
  interaction_zone = interaction_zone[nearest_interaction_zone, ]
  
  # 4. Clip interaction zone by 550m buffer from flowline (90% dispersal distance)
  max_dispersal_zone = st_union(st_buffer(flowlines, 550))
  interaction_zone = st_intersection(interaction_zone, st_transform(max_dispersal_zone, crs_standard))
  # EDIT: Retain the full riparian zone
  interaction_zone = st_union(interaction_zone, st_transform(riparian_zone_no_waterbodies, crs_standard))
  
  # 5. Clip interaction zone by site buffer
  interaction_zone = st_intersection(interaction_zone, st_transform(site_buffer, crs_standard))
  interaction_zone = st_cast(st_union(interaction_zone), "POLYGON")
  interaction_zone = interaction_zone[which.max(st_area(interaction_zone))] # clean up small fragments
  area_interaction_zone = st_area(interaction_zone)
  
  # Inspect interaction zone and component layers
  # mapview(site, col.regions = "purple", layer.name = paste0("site_", s_id)) +
  #   mapview(site_pssb, col.regions = "orange") +
  #   mapview(site_buffer, alpha.regions = 0, lwd = 3, col.regions = "black", layer.name = "buffer_reach") +
  #   # mapview(waterbodies, col.regions = "red") +
  #   # mapview(waterbodies_buffer, col.regions = "pink") +
  #   mapview(flowlines, zcol = "f_code") +
  #   # mapview(max_dispersal_zone, col.regions = "gray", alpha = 0.2) +
  #   mapview(riparian_zone, col.regions = "yellow") +
  #   # mapview(emergence_zone, col.regions = "blue") +
  #   mapview(interaction_zone, col.regions = "green")

  if (s_id == "121") { # Figure
    emergence_zone = st_intersection(emergence_zone, site_buffer %>% st_transform(st_crs(emergence_zone)))
    
    imp = transform_crop_and_mask(rast_data$rast_nlcd_impervious, site_buffer, site_buffer)
    imp = project(imp, crs(site_buffer))
    imp = as.data.frame(imp, xy = TRUE)
    
    canopy = transform_crop_and_mask(rast_data$rast_usfs_canopycover, site_buffer, site_buffer)
    canopy = project(canopy, crs(site_buffer))
    canopy = as.data.frame(canopy, xy = TRUE)
    
    p_site = ggplot() +
      geom_sf(data = site_buffer, color = "grey10", fill = "transparent", linewidth = 0.5) +
      
      # geom_raster(data = imp, aes(x = x, y = y, fill = rast_nlcd_impervious)) +
      # scale_fill_gradient(low = alpha("white", 0.1), high = alpha("red", 1.0)) +
      # new_scale_fill() +
      # geom_raster(data = canopy, aes(x = x, y = y, fill = rast_usfs_canopycover)) +
      # scale_fill_gradient(low = alpha("white", 0.1), high = alpha("forestgreen", 1.0)) +
      # geom_sf(data = interaction_zone, fill = "white", color = "white", linetype = "dashed", linewidth = 0.5, alpha = 0.5) +
      
      geom_sf(data = riparian_zone, fill = "dodgerblue", color = "transparent", alpha = 0.5) +
      geom_sf(data = interaction_zone, fill = "transparent", color = "dodgerblue", linewidth = 0.5, alpha = 0) +
      geom_sf(data = flowlines, color = "dodgerblue", linewidth = 1, alpha = 0.75) +
      geom_sf(data = waterbodies, fill = "dodgerblue", color = "transparent", linewidth = 0, alpha = 0.75) +
      geom_sf(data = site, fill = "red", color = "grey10", shape = 21, size = 3, alpha = 0.75) +
      # geom_sf(data = site_pssb, fill = "blue", color = "grey10", shape = 21, size = 3, alpha = 0.75) +
      geom_sf(data = site_buffer, color = "grey70", fill = "transparent", linewidth = 0.5) +
      geom_sf(data = st_buffer(site, 550), color = "grey70", fill = "transparent", linewidth = 0.5) +
      coord_sf(expand = FALSE) +
      labs(x = "", y = "") +
      theme_sleek() + theme(
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ); print(p_site)

    ggsave(paste0("p_site_", buffer_radius, ".pdf"), p_site, width = 8, height = 6, dpi = 300)
  }
  
  vectors = list(
    "site_buffer"        = site_buffer,
    "interaction_zone"   = interaction_zone,
    "riparian_zone"      = riparian_zone
    # "emergence_zone"     = emergence_zone,
    # "max_dispersal_zone" = max_dispersal_zone,
    # "flowlines"          = flowlines
  )
  
  site_zones[[s_id]] = vectors
  
  pb$tick()
}

# Calculate impervious projections --------------------------------------------------
if (TRUE) {
  
  message("Calculating ICLUS impervious projections")
  # Baselines
  imp_projection_results = data.frame()
  i = 1
  for (s_id in names(site_zones)) {
    
    if (s_id %in% sites_to_exclude) next
    
    site_buffer = site_zones[[s_id]][["site_buffer"]]
    interaction_zone = site_zones[[s_id]][["interaction_zone"]]
    
    # NLCD "ground-truth"
    site_nlcd_impervious  = transform_crop_and_mask(rast_data$rast_nlcd_impervious, site_buffer, interaction_zone)
    stats_nlcd_impervious = rast_stats(site_nlcd_impervious)
    
    # Projections
    message("Calculating projections for site ", s_id)
    imp_projection_baselines = rbind(cbind(data.frame(site_id = s_id, scenario = "NLCD", year = "2023", stats_nlcd_impervious)))
    for (scenario in names(rast_imp_iclus_projections)) {
      for (year in as.character(seq(2020, 2100, 10))) {
        r = rast_imp_iclus_projections[[scenario]][[year]]
        r = disagg(r, fact = round(res(r)[1] / 30), method = "near")
        r = transform_crop_and_mask(r, site_buffer, interaction_zone)
        stats_imp_iclus_projections = rast_stats(r)
        imp_projection_baselines = rbind(
          imp_projection_baselines,
          cbind(data.frame(site_id = s_id, scenario = scenario, year = year), stats_imp_iclus_projections)
        )
      }
    }
    imp_projection_results = rbind(imp_projection_results, imp_projection_baselines)
    # print(imp_projection_baselines)
    message(round(i / length(names(site_zones)), 2) * 100, "%")
    i = i + 1
  }
  
  # Difference in observed NLCD imperviousness (2023) and ICLUS projections
  p = imp_projection_results %>%
    filter(year %in% c("2020", "2023", "2030", "2050", "2100")) %>%
    group_by(site_id) %>%
    mutate(nlcd_sum_proportion = sum_proportion[scenario == "NLCD"]) %>%
    mutate(diff_from_nlcd = sum_proportion - nlcd_sum_proportion) %>%
    ungroup() %>%
    filter(year %in% c("2020", "2030", "2050", "2100")) %>%
    ggplot(aes(x = factor(year), y = diff_from_nlcd, fill = scenario)) +
    geom_hline(yintercept = 0) +
    geom_boxplot(position = "dodge") +
    geom_text_repel(aes(label = site_id)) +
    facet_wrap(~scenario) +
    labs(title = "Difference in observed NLCD imperviousness (2023) and ICLUS projections"); print(p)
  
  # Remove outliers
  s_id_outliers = c(216,213,207,210)
  imp_diff = imp_projection_results %>%
    filter(!(site_id %in% s_id_outliers)) %>%
    filter(year %in% c("2020", "2023")) %>%
    group_by(site_id) %>%
    mutate(nlcd_sum_proportion = sum_proportion[scenario == "NLCD"]) %>%
    mutate(diff_from_nlcd = sum_proportion - nlcd_sum_proportion) %>%
    ungroup() %>%
    filter(year %in% c("2020"))
  
  message("Difference in observed NLCD imperviousness (2023) and ICLUS projections after removing outliers")
  print(tapply(imp_diff$diff_from_nlcd, imp_diff$scenario, summary))
  
  # NOTE: Sites with extremely high imp projections are also those that have underestimated imp
  # coverage in 2020/30 projections compared to NLCD baseline, all along Renton WA 515 corridor
  mapview(sites_aru %>% filter(site_id %in% s_id_outliers))
  # Exclude these outliers from analysis
  imp_projection_results_clean = imp_projection_results %>% filter(!(site_id %in% s_id_outliers))
  
  # Treat 2020 as baseline to compute % change per scenario, per site
  imp_projections = imp_projection_results_clean %>%
    select(site_id, scenario, year, sum_proportion) %>%
    group_by(site_id, scenario) %>%
    mutate(percent_change_sum_prop = ifelse(
      scenario != "NLCD",
      sum_proportion / sum_proportion[year == 2020],
      NA_real_
    )) %>%
    ungroup() %>% arrange(site_id, scenario, year) %>%
    filter(year != 2020, scenario != "NLCD") %>% select(-sum_proportion)
  imp_projections$sum_proportion = NA
  
  nlcd_baselines = imp_projection_results_clean %>%
    filter(scenario == "NLCD") %>%
    select(site_id, scenario, year, sum_proportion) %>%
    rename(original_scenario = scenario) %>% 
    crossing(scenario = unique(imp_projections$scenario)) %>%
    select(-original_scenario) %>% mutate(percent_change_sum_prop = 1.0)
  
  imp_projection_estimates = imp_projections %>%
    left_join(
      nlcd_baselines %>% select(site_id, scenario, baseline_sum = sum_proportion),
      by = c("site_id", "scenario")
    ) %>%
    mutate(sum_proportion = ifelse(
        is.na(sum_proportion),
        baseline_sum * percent_change_sum_prop,
        sum_proportion
      )) %>%
    select(-baseline_sum) %>% rbind(nlcd_baselines) %>% arrange(site_id, scenario, year)
  
  {
    # dbg_site_id = "216"
    # site_buffer = site_zones[[dbg_site_id]][["site_buffer"]]
    # interaction_zone = site_zones[[dbg_site_id]][["interaction_zone"]]
    # dbg_rast = transform_crop_and_mask(rast_data$rast_nlcd_impervious,
    #                                    site_buffer,
    #                                    interaction_zone)
    # r_2020 = rast_imp_iclus_projections[["a1"]][["2020"]]
    # r_2020 = disagg(r_2020, fact = round(res(r_2020)[1] / 30), method = "near")
    # r_2020 = transform_crop_and_mask(r_2020, site_buffer, interaction_zone)
    # r_2100 = rast_imp_iclus_projections[["a1"]][["2100"]]
    # r_2100 = disagg(r_2100, fact = round(res(r_2100)[1] / 30), method = "near")
    # r_2100 = transform_crop_and_mask(r_2100, site_buffer, interaction_zone)
    # mapview(dbg_rast) + mapview(r_2020) + mapview(r_2100)
  }
  
  imp_projection_estimates %>% arrange(desc(sum_proportion))
  
  p = imp_projection_estimates %>%
    ggplot(aes(x = year, y = sum_proportion, color = scenario, group = interaction(site_id, scenario))) +
    # geom_hline(yintercept = 1.0) +
    geom_line(alpha = 0.25) +
    facet_wrap(~scenario) +
    stat_summary(aes(group = scenario, color = scenario), fun = mean, geom = "line", size = 1.2) +
    labs(title = "Predicted impervious surface % per site from % change relative to ICLUS 2020"); print(p)
  
  # Cache ICLUS impervious projections ----------------------------------------
  
  if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
  out_filepath = file.path(out_cache_dir, paste0("imp_projections_", buffer_radius, "m.rds"))
  saveRDS(imp_projection_estimates, out_filepath)
  message(crayon::green("Cached", out_filepath))
}

## Calculate stats for continuous raster data
geospatial_site_data = list()
message("Calculating stats for continuous raster data within interaction zones")
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
for (s_id in names(site_zones)) {
  
  site_buffer = site_zones[[s_id]][["site_buffer"]]
  interaction_zone = site_zones[[s_id]][["interaction_zone"]]
  riparian_zone = site_zones[[s_id]][["riparian_zone"]]
  area_usfs_riparian_zone = st_area(st_as_sf(riparian_zone))
  if (length(area_usfs_riparian_zone) == 0) area_usfs_riparian_zone = 0
  area_interaction_zone = st_area(interaction_zone)
  
  # NLCD imperviousness
  site_nlcd_impervious  = transform_crop_and_mask(rast_data$rast_nlcd_impervious, site_buffer, interaction_zone)
  stats_nlcd_impervious = rast_stats(site_nlcd_impervious)
  mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_nlcd_impervious)
  
  # NASA GEDI foliage height diversity
  site_gedi_fhd  = transform_crop_and_mask(rast_data$rast_gedi_fhd, site_buffer, interaction_zone)
  stats_gedi_fhd = rast_stats(site_gedi_fhd)
  # mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_gedi_fhd)
  
  # NASA GEDI canopy cover
  site_gedi_cover  = transform_crop_and_mask(rast_data$rast_gedi_cover, site_buffer, interaction_zone)
  stats_gedi_cover = rast_stats(site_gedi_cover)
  
  # NASA GEDI canopy height
  site_gedi_height  = transform_crop_and_mask(rast_data$rast_gedi_height, site_buffer, interaction_zone)
  stats_gedi_height = rast_stats(site_gedi_height)
  
  # NASA GEDI proportion of vegetation density (mature upper canopy)
  site_gedi_pavd20m  = transform_crop_and_mask(rast_data$rast_gedi_pavd20m, site_buffer, interaction_zone)
  stats_gedi_pavd20m = rast_stats(site_gedi_pavd20m)
  
  # NASA GEDI proportion of vegetation density (mid canopy)
  site_gedi_pavd5to10m  = transform_crop_and_mask(rast_data$rast_gedi_pavd5to10m, site_buffer, interaction_zone)
  stats_gedi_pavd5to10m = rast_stats(site_gedi_pavd5to10m)
  
  # USFS canopy cover
  site_usfs_canopycover  = transform_crop_and_mask(rast_data$rast_usfs_canopycover, site_buffer, interaction_zone)
  stats_usfs_canopycover = rast_stats(site_usfs_canopycover)
  
  # Landfire tree height
  site_landfire_treeheight  = transform_crop_and_mask(rast_data$rast_landfire_treeheight, site_buffer, interaction_zone)
  stats_landfire_treeheight = rast_stats(site_landfire_treeheight)
  
  site_landfire_shrubheight  = transform_crop_and_mask(rast_data$rast_landfire_shrubheight, site_buffer, interaction_zone)
  stats_landfire_shrubheight = rast_stats(site_landfire_shrubheight)
  
  site_landfire_herbheight  = transform_crop_and_mask(rast_data$rast_landfire_herbheight, site_buffer, interaction_zone)
  stats_landfire_herbheight = rast_stats(site_landfire_herbheight)
  
  # LEMMA forest structure
  site_lemma_age  = transform_crop_and_mask(rast_data$rast_lemma_age, site_buffer, interaction_zone)
  site_lemma_age[site_lemma_age < 0] = NA
  stats_lemma_age = rast_stats(site_lemma_age)
  
  site_lemma_ba  = transform_crop_and_mask(rast_data$rast_lemma_ba, site_buffer, interaction_zone)
  site_lemma_ba[site_lemma_ba < 0] = NA
  stats_lemma_ba = rast_stats(site_lemma_ba)

  site_lemma_qmd  = transform_crop_and_mask(rast_data$rast_lemma_qmd, site_buffer, interaction_zone)
  site_lemma_qmd[site_lemma_qmd < 0] = NA
  stats_lemma_qmd = rast_stats(site_lemma_qmd)
  
  site_lemma_denall  = transform_crop_and_mask(rast_data$rast_lemma_denall, site_buffer, interaction_zone)
  site_lemma_denall[site_lemma_denall < 0] = NA
  stats_lemma_denall = rast_stats(site_lemma_denall)
  
  site_lemma_dencon  = transform_crop_and_mask(rast_data$rast_lemma_dencon, site_buffer, interaction_zone)
  site_lemma_dencon[site_lemma_dencon < 0] = NA
  stats_lemma_dencon = rast_stats(site_lemma_dencon)
  
  site_lemma_denhw  = transform_crop_and_mask(rast_data$rast_lemma_denhw, site_buffer, interaction_zone)
  site_lemma_denhw[site_lemma_denhw < 0] = NA
  stats_lemma_denhw = rast_stats(site_lemma_denhw)
  
  site_lemma_domts = transform_crop_and_mask(rast_data$rast_lemma_domts, site_buffer, interaction_zone)
  site_lemma_domts[site_lemma_domts < 0] = NA
  tree_species_richness = length(na.omit(unique(values(site_lemma_domts))))
  
  # Explore masking structural vars to non-urban areas
  site_gedi_height_impmask = mask(project(site_gedi_height, site_nlcd_impervious), site_nlcd_impervious <= 20, maskvalues = TRUE, inverse = TRUE)
  names(site_gedi_height_impmask) = "rast_gedi_height_impmask"
  stats_gedi_height_impmask = rast_stats(site_gedi_height_impmask)
  
  ## Calculate composition for categorical NLCD land cover raster data
  
  site_nlcd_landcover_site_buffer = transform_crop_and_mask(rast_data$rast_nlcd_landcover, site_buffer, site_buffer)
  names(site_nlcd_landcover_site_buffer) = "rast_nlcd_landcover_site_buffer"
  site_nlcd_landcover = transform_crop_and_mask(rast_data$rast_nlcd_landcover, site_buffer, interaction_zone)
  
  # mapview(as.factor(site_nlcd_landcover), col.region = nlcd_metadata %>%
  #           filter(id %in% unique(values(site_nlcd_landcover))) %>% select(id, class, color) %>% pull(color))

  # Calculate relative abundance of each cover class
  freq_table = freq(site_nlcd_landcover) %>% select(value, count)
  freq_table$percent = (freq_table$count / sum(freq_table$count)) 
  cover_detail = full_join(freq_table %>% rename(id = value), nlcd_metadata %>% select(id, class), by = "id")
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
  
  cover_summary = cover_summary %>% # Forest AND wetlands
    bind_rows(
      cover_summary %>%
        filter(group %in% c("Forest", "Wetlands")) %>%
        summarise(
          group = "Forest and Wetlands",
          count = sum(count),
          percent = sum(percent)
        )
    )
  
  # Major landcover classes
  r = site_nlcd_landcover
  # Reclassify summary group "Developed, Variable Intensity"
  r[r == 22] = 2
  r[r == 23] = 2
  r[r == 24] = 2
  # Reclassify summary group "Forest"
  r[r == 41] = 4
  r[r == 42] = 4
  r[r == 43] = 4
  # Reclassify summary group "Wetlands"
  r[r == 90] = 9
  r[r == 95] = 9
  # Reclassify summary group "Agriculture"
  r[r == 81] = 8
  r[r == 82] = 8
  r[] = as.integer(r[])
  site_nlcd_landcover_major = r
  
  # Reclassify simple landcover classes for landscape configuration metrics
  r = site_nlcd_landcover
  # Reclassify specific values
  r[r %in% c(41, 42, 43, 90, 95)] = 1 # "Riparian habitat" (Forested or Wetland)
  r[r %in% c(21, 22, 23, 24)]     = 2 # Developed
  r[r == 11]                      = 3 # Open water
  r[!r %in% c(1, 2, 3, NA)]       = 4 # Other
  
  ## Landscape configuration metrics
  
  ed = lsm_l_ed(r)$value # Edge density
  units(ed) = "m/ha"
  edge_density = lsm_c_ed(r)
  ed_riphab = edge_density %>% filter(class == 1) %>% pull(value)
  ed_dev    = edge_density %>% filter(class == 2) %>% pull(value)
  if (length(ed_riphab) == 0) ed_riphab = 0
  if (length(ed_dev) == 0)    ed_dev = 0
  units(ed_riphab) = "m/ha"
  units(ed_dev)    = "m/ha"

  agg = lsm_l_ai(r)$value
  aggregation = lsm_c_ai(r) # Aggregation of a class, 0 maximally disaggregated, 100 maximally aggregated
  agg_riphab = aggregation %>% filter(class == 1) %>% pull(value) # Riparian habitat
  agg_dev    = aggregation %>% filter(class == 2) %>% pull(value) # Developed
  if (length(agg_riphab) == 0) agg_riphab = NA
  if (length(agg_dev) == 0)    agg_dev = NA
  
  pd = lsm_l_pd(r)$value
  patch_density = lsm_c_pd(r) # Fragmentation of a class, increases as landscape gets more patchy
  pd_riphab = patch_density %>% filter(class == 1) %>% pull(value)
  pd_dev    = patch_density %>% filter(class == 2) %>% pull(value)
  if (length(pd_riphab) == 0) pd_riphab = 0
  if (length(pd_dev) == 0)    pd_dev = 0

  # mapview(as.factor(site_nlcd_landcover), col.region = nlcd_metadata %>%
  #           filter(id %in% unique(values(site_nlcd_landcover))) %>% select(id, class, color) %>% pull(color))
  
  # Roads
  # S1100 - Primary road 
  # Primary roads are limited-access highways that connect to other roads only at interchanges and not at at-grade intersections
  # S1200 - Secondary Road
  # Secondary roads are main arteries that are not limited access, usually in the U.S. highway, state highway, or county highway systems
  # S1400 - Local Neighborhood Road, Rural Road, City Street
  # Generally a paved non-arterial street, road, or byway that usually has a single lane of traffic in each direction.
  # site_roads = st_intersection(roads_major, st_transform(site_buffer, crs = st_crs(roads_major)))
  # dist_road_major = st_distance(st_transform(site, crs = st_crs(roads_major)), roads_major)
  # dist_road_major = min(dist_road_major)
  # dist_road_paved = st_distance(st_transform(site, crs = st_crs(roads_paved)), roads_paved) # NOTE: time intensive
  # dist_road_paved = min(dist_road_paved)
  
  site_roads_paved = st_crop(roads_paved, st_transform(interaction_zone, crs = st_crs(roads_paved)))
  site_roads_paved = st_intersection(st_transform(site_roads_paved, st_crs(interaction_zone)), interaction_zone) %>% st_transform(crs_standard)
  site_roads_major = st_crop(roads_major, st_transform(interaction_zone, crs = st_crs(roads_major)))
  site_roads_major = st_intersection(st_transform(site_roads_major, st_crs(interaction_zone)), interaction_zone) %>% st_transform(crs_standard)
  # mapview(site, col.regions = "green") + 
  #   mapview(site_pssb, col.regions = "blue") +
  #   mapview(site_buffer, alpha.regions = 0, lwd = 3) +
  #   mapview(flowlines, zcol = "f_code") +
  #   mapview(site_roads_paved, zcol = "mtfcc")
  
  # Density (m per square km) = total road length / buffer area (in m2) / 1e6
  density_roads_paved = sum(st_length(site_roads_paved)) / area_interaction_zone
  density_roads_major = sum(st_length(site_roads_major)) / area_interaction_zone
  density_roads_paved = set_units(density_roads_paved, "m/km^2")
  density_roads_major = set_units(density_roads_major, "m/km^2")
  
  # Store all data for the site
  rasters = list(
    site_riparian_usfs,
    site_nlcd_impervious,
    site_gedi_fhd,
    site_gedi_cover,
    site_gedi_height,
    site_gedi_pavd20m,
    site_gedi_pavd5to10m,
    site_usfs_canopycover,
    site_landfire_treeheight,
    site_gedi_height_impmask,
    site_nlcd_landcover,
    site_nlcd_landcover_major,
    site_nlcd_landcover_site_buffer,
    site_lemma_age,
    site_lemma_ba,
    site_lemma_denall,
    site_lemma_dencon,
    site_lemma_denhw,
    site_lemma_domts,
    site_lemma_qmd
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
    stats_landfire_treeheight,
    stats_gedi_height_impmask,
    stats_lemma_age,
    stats_lemma_ba,
    stats_lemma_denall,
    stats_lemma_dencon,
    stats_lemma_denhw,
    stats_lemma_qmd
  )) %>% rownames_to_column(var = "name")
  
  vectors = list(
    # "site_riparian_kc" = site_riparian_kc,
    "interaction_zone"   = interaction_zone,
    # "riparian_zone"      = riparian_zone,
    # "emergence_zone"     = emergence_zone,
    # "max_dispersal_zone" = max_dispersal_zone,
    # "flowlines"          = flowlines,
    "roads_paved"        = site_roads_paved,
    "roads_major"        = site_roads_major
  )
  
  geospatial_site_data[[s_id]] = list(
    rasters = rasters,
    raster_stats = raster_stats,
    nlcd_cover = list(
      "nlcd_cover_detail" = cover_detail,
      "nlcd_cover_summary" = cover_summary
    ),
    vectors = vectors,
    f_code = f_code,
    ed        = ed,
    ed_riphab = ed_riphab,
    ed_dev    = ed_dev,
    agg        = agg,
    agg_riphab = agg_riphab,
    agg_dev    = agg_dev,
    pd        = pd,
    pd_riphab = pd_riphab,
    pd_dev    = pd_dev,
    tree_species_richness = tree_species_richness,
    # dist_road_major = dist_road_major,
    # dist_road_paved = dist_road_paved,
    density_roads_paved = density_roads_paved,
    density_roads_major = density_roads_major,
    area_riparian_usfs = area_usfs_riparian_zone,
    # area_riparian_kc = area_riparian_kc,
    area_interaction_zone = area_interaction_zone
  )
  pb$tick()
}

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

additional_metrics = do.call(rbind, lapply(names(geospatial_site_data), function(i) {
  print(i)
  data.frame(
    site_id               = i,
    f_code                = geospatial_site_data[[i]]$f_code,
    ed                    = geospatial_site_data[[i]]$ed,
    ed_riphab             = geospatial_site_data[[i]]$ed_riphab,
    ed_dev                = geospatial_site_data[[i]]$ed_dev,
    agg                   = geospatial_site_data[[i]]$agg,
    agg_riphab            = geospatial_site_data[[i]]$agg_riphab,
    agg_dev               = geospatial_site_data[[i]]$agg_dev,
    pd                    = geospatial_site_data[[i]]$pd,
    pd_riphab             = geospatial_site_data[[i]]$pd_riphab,
    pd_dev                = geospatial_site_data[[i]]$pd_dev,
    tree_species_richness = geospatial_site_data[[i]]$tree_species_richness,
    density_roads_paved   = geospatial_site_data[[i]]$density_roads_paved,
    density_roads_major   = geospatial_site_data[[i]]$density_roads_major,
    area_riparian_usfs    = geospatial_site_data[[i]]$area_riparian_usfs,
    area_interaction_zone = geospatial_site_data[[i]]$area_interaction_zone,
    stringsAsFactors = FALSE
  )
}))

# Join land cover data
site_data = site_data %>% left_join(nlcd_summary, by = "site_id")
site_data = site_data %>% left_join(additional_metrics, by = "site_id")

# Extract raster stats for each site in wide format
stats = lapply(names(geospatial_site_data), function(i) {
  df = geospatial_site_data[[i]]$raster_stats
  df$site_id = i
  return(df)
}) %>% bind_rows() %>% pivot_wider(
  id_cols = site_id,
  names_from = name,
  values_from = c(max, min, sum, mean, sd, cv, sum_proportion),
  names_glue = "{name}_{.value}"
)

# Join raster stats
site_data = site_data %>% left_join(stats, by = "site_id")

# Visualize geospatial data ---------------------------------------------------------------------
message("Visualizing geospatial data")

# Order sites by development intensity
site_order = unique(site_data %>% arrange(desc(nlcd_developed_variable_intensity)) %>% pull(site_id))

nlcd_cols = grep("^nlcd_", names(site_data), value = TRUE)
nlcd_cols = setdiff(nlcd_cols, "nlcd_forest_and_wetlands")
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
    title = paste0("NLCD land cover composition (", buffer_radius, " m buffer)")
  ) +
  scale_fill_manual(values = c(
    "developed_variable_intensity" = "#eb0000",
    "developed_open_space"         = "#eb9999",
    "forest"                       = "#1c5f2c",
    "wetlands"                     = "#6c9fb8",
    "agriculture"                  = "#ead963",
    "barren_land"                  = "#b3ac9f",
    "grassland_herbaceous"         = "#dde9af",
    "open_water"                   = "#466b9f",
    "perennial_ice_snow"           = "#d1defa",
    "shrub_scrub"                  = "#af963c"
  )) +
  theme_minimal()

# Visualize land cover for all sites

sites_lc_df = bind_rows(
  lapply(names(geospatial_site_data), function(site_id) { # OR rast_nlcd_landcover (interaction zone only)
    if (site_id %in% c(sites_to_exclude)) return(data.frame())
    df = as.data.frame(geospatial_site_data[[site_id]]$rasters$rast_nlcd_landcover_site_buffer, xy = TRUE)
    value_col = names(df)[3]
    df = df %>% rename(value = all_of(value_col))  # rename for ggplot
    df$site <- site_id
    df$class = nlcd_metadata[match(df$value, nlcd_metadata$id), "class"]
    return(df)
  })
) %>% mutate(site = factor(site, levels = site_order))
present_levels = nlcd_metadata %>% filter(id %in% unique(sites_lc_df$value))

# Order sites by increasing amount of rast_nlcd_impervious_sum_proportion, nlcd_developed_variable_intensity, etc.
site_order = site_data %>% arrange((rast_nlcd_impervious_sum_proportion)) %>% pull(site_id)
sites_lc_df = sites_lc_df %>%
  mutate(site = factor(site, levels = site_order))

# Overwrite colors
# scales::show_col(present_levels$color)
# present_levels$color[1] = "#97B0BD"
# present_levels$color[10] = "purple"
# present_levels$color[11] = "purple"

p_sites = ggplot(sites_lc_df, aes(x = x, y = y, fill = factor(class))) + geom_raster() +
  scale_fill_manual(values = setNames(present_levels$color, present_levels$class), name = "Class") +
  facet_wrap(~site, scales = "free", ncol = 10, strip.position = "bottom") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD land cover configuration (", buffer_radius, " m buffer)"), x = "", y = "") +
  theme(
    strip.text = element_text(size = 6),
    panel.spacing = unit(0, "mm"),
    plot.margin = margin(0, 0, 0, 0)
  ); print(p_sites)

ggsave(paste0("site_landcover_", buffer_radius, "m.png"), p_sites, dpi = 300, width = 10, height = 5.3)

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
  labs(title = paste0("NLCD fractional impervious surface (", buffer_radius, " m buffer)"), x = "", y = "", fill = "Impervious %")

# Visualize canopy cover for all sites
sites_tcc_df = bind_rows(
  lapply(seq_along(geospatial_site_data), function(i) {
    raster_to_df(geospatial_site_data[[i]]$rasters$rast_usfs_canopycover, names(geospatial_site_data)[i])
  })
) %>% mutate(site = factor(site, levels = site_order))

ggplot(sites_tcc_df, aes(x = x, y = y, fill = class)) + geom_tile() +
  facet_wrap(~ site, scales = "free") +
  scale_fill_gradient(low = "lightgray", high = "darkgreen") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("USFS tree canopy cover (", buffer_radius, " m buffer)"), x = "", y = "", fill = "Cover %")

# Area of interaction zones
ggplot(site_data, aes(x = reorder(site_id, area_interaction_zone), y = area_interaction_zone)) + geom_col()
hist(site_data$area_interaction_zone, breaks = 20)

# Area of USFS riparian zones
ggplot(site_data, aes(x = reorder(site_id, area_riparian_usfs), y = area_riparian_usfs)) + geom_col()
hist(site_data$area_riparian_usfs, breaks = 20)

## Visualize relationships between raster summary stats

# Site attributes (e.g. BIBI)
ggplot(site_data, aes(x = rast_nlcd_impervious_sum_proportion, y = bibi)) +
  geom_point() + geom_smooth(method = "lm") + geom_text_repel(aes(label = site_id)) 

# # Comparing different rasters
# ggplot(site_data, aes(x = rast_landfire_treeheight_sum, y = rast_gedi_height_sum)) +
#   geom_point() + geom_abline(intercept = 0, slope = 1) +
#   geom_smooth(method = "lm") + geom_text_repel(aes(label = site_id))

# Cache site data --------------------------------------------------------------------

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
out_filepath = file.path(out_cache_dir, paste0(
  "site_data_", buffer_radius, "m.rds"))
saveRDS(site_data, out_filepath)
message(crayon::green("Cached", out_filepath))
