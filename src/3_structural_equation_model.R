# 3_explore_data.R ===================================================================
# Load ARU and PSSB site locations, calculate all site covariates from land cover
# data, and do some exploratory structural equation modeling
#
# Inputs:
overwrite_geospatial_site_data = TRUE
in_path_site_metadata   = "data/site_metadata.csv"
in_path_pssb_data       = "data/raw/pssb/ScoresByYear.csv"
in_path_nlcd_metadata   = "data/raw/nlcd_metadata.csv"
in_cache_detections     = "data/cache/1_preprocess_agg_pam_data/detections.rds"
in_cache_geospatial_dir = "data/cache/2_preprocess_geospatial_data"
in_path_species_list    = "data/pam/species_list.txt"
in_path_avonet_traits   = "data/traits/AVONET Supplementary dataset 1.xlsx"
in_path_elton_traits    = "data/traits/EltonTraits 1.0/BirdFuncDat.csv"
# Outputs:
out_cache_dir = "data/cache/3_structural_equation_model"

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
  "landscapemetrics", # landscape metrics
  "piecewiseSEM", # structural equation modeling
  "units"
)
sapply(pkgs, function(pkg) {
  if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
  as.character(packageVersion(pkg))
})

buffer_radius = set_units(150, m) # 5km site buffer (alternative: ~500 m insect emergence 90% flux range falloff)

theme_set(theme_minimal())

# Load ARU and PSSB site locations and define study area -----------------------------
message("Loading ARU and PSSB site locations")

site_metadata = read_csv(in_path_site_metadata, show_col_types = FALSE) %>% clean_names() %>% mutate(site_id = as.character(site_id))

crs_standard = "EPSG:32610" # shared coordinate reference system (metric)

# Create sf points for site locations (reference original crs 4326)
sites_aru = site_metadata %>%
  st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326)
sites_pssb = site_metadata %>%
  st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326)

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

# Load and process geospatial data ---------------------------------------------------
message("Loading geospatial data")

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

if (overwrite_geospatial_site_data) {

  geospatial_site_data = list()
  pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = nrow(sites_aru), clear = FALSE)
  for (s in 1:nrow(sites_aru)) {
  
    s_id = sites_aru$site_id[s]
    # print(s)
    # print(s_id)
    
    # Create a buffer around the ARU site # TODO: Center around ARU or PSSB?
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
    
    # 2. Buffer 50 m core emergence zone minimum from flowlines
    # TODO: Buffer 100m to include adjacent edge habitats, as below?
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
    
    # Only select nearest interaction zone
    interaction_zone = st_cast(interaction_zone, "POLYGON")
    nearest_interaction_zone = st_nearest_feature(st_transform(site, crs_standard), interaction_zone)
    interaction_zone = interaction_zone[nearest_interaction_zone, ]
    
    # 4. Clip interaction zone by 550m buffer from flowline (90% dispersal distance)
    max_dispersal_zone = st_union(st_buffer(flowlines, 550))
    interaction_zone = st_intersection(interaction_zone, st_transform(max_dispersal_zone, crs_standard))
    
    # 5. Clip interaction zone by site buffer
    interaction_zone = st_intersection(interaction_zone, st_transform(site_buffer, crs_standard))
    interaction_zone = st_cast(st_union(interaction_zone), "POLYGON")
    interaction_zone = interaction_zone[which.max(st_area(interaction_zone))] # clean up small fragments
    area_interaction_zone = st_area(interaction_zone)
    
    # Inspect interaction zone and component layers
    mapview(site, col.regions = "purple", layer.name = paste0("site_", s_id)) +
      mapview(site_pssb, col.regions = "orange") +
      mapview(site_buffer, alpha.regions = 0, lwd = 3, col.regions = "black", layer.name = "buffer") +
      # mapview(waterbodies, col.regions = "red") +
      # mapview(waterbodies_buffer, col.regions = "pink") +
      mapview(flowlines, zcol = "f_code") +
      mapview(max_dispersal_zone, col.regions = "gray", alpha = 0.2) +
      # mapview(riparian_zone, col.regions = "yellow") +
      mapview(emergence_zone, col.regions = "blue") +
      mapview(interaction_zone, col.regions = "green")
    
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
    # TODO: THIS ASSUMES A CIRCULAR AREA, NOT THE VARIABLE INTERACTION ZONE!
    density_roads_paved = as.numeric(sum(st_length(site_roads_paved)) / ((pi * buffer_radius^2) / set_units(1, km^2)))
    density_roads_major = as.numeric(sum(st_length(site_roads_major)) / ((pi * buffer_radius^2) / set_units(1, km^2)))
    units(density_roads_paved) <- "m/km^2"
    units(density_roads_major) <- "m/km^2"
    
    ## Calculate stats for continuous raster data
    
    # NLCD imperviousness
    site_nlcd_impervious  = transform_crop_and_mask(rast_data$rast_nlcd_impervious, site_buffer, interaction_zone)
    stats_nlcd_impervious = rast_stats(site_nlcd_impervious)
    # mapview(site) + mapview(site_buffer, alpha.regions = 0, lwd = 2) + mapview(site_nlcd_impervious)
    
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
    
    # Simple landcover classes for edge density
    r = site_nlcd_landcover
    # Reclassify specific values
    r[r == 11] <- 1 # Open water
    r[r %in% c(21, 22, 23, 24)] <- 2 # Developed 
    r[r %in% c(41, 42, 43, 90)] <- 3 # Forested
    # Set all other values to 0
    r[!r %in% c(1, 2, 3, NA)] <- 4 # Other
    
    # Landscape level edge density
    edge_density = lsm_l_ed(r)$value
    units(edge_density) = "m/ha"
    # contagion = lsm_l_contag(r)$value
    aggregation = lsm_l_ai(r)$value
    
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
      "riparian_zone"      = riparian_zone,
      "emergence_zone"     = emergence_zone,
      "max_dispersal_zone" = max_dispersal_zone,
      "flowlines"          = flowlines,
      "roads_paved"        = site_roads_paved,
      "roads_major"        = site_roads_major
    )
    
    geospatial_site_data[[as.character(site$site_id)]] = list(
      rasters = rasters,
      raster_stats = raster_stats,
      nlcd_cover = list(
        "nlcd_cover_detail" = cover_detail,
        "nlcd_cover_summary" = cover_summary
      ),
      vectors = vectors,
      edge_density = edge_density,
      aggregation = aggregation,
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
  
  # Cache results
  if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
  out_filepath = file.path(out_cache_dir, "geospatial_site_data.rds")
  saveRDS(geospatial_site_data, out_filepath)
  message(crayon::green("Cached", out_filepath))
  
} else {
  message("Loading geospatial site data from cache")
  in_filepath = file.path(out_cache_dir, "geospatial_site_data.rds")
  geospatial_site_data = readRDS(in_filepath)
  message(crayon::green("Loaded", in_filepath))
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
    edge_density          = geospatial_site_data[[i]]$edge_density,
    aggregation           = geospatial_site_data[[i]]$aggregation,
    tree_species_richness = geospatial_site_data[[i]]$tree_species_richness,
    # dist_road_major       = geospatial_site_data[[i]]$dist_road_major,
    # dist_road_paved       = geospatial_site_data[[i]]$dist_road_paved,
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
    df = as.data.frame(geospatial_site_data[[site_id]]$rasters$rast_nlcd_landcover, xy = TRUE)
    value_col = names(df)[3]
    df = df %>% rename(value = all_of(value_col))  # rename for ggplot
    df$site <- site_id
    df$class = nlcd_metadata[match(df$value, nlcd_metadata$id), "class"]
    return(df)
  })
) %>% mutate(site = factor(site, levels = site_order))
present_levels = nlcd_metadata %>% filter(id %in% unique(sites_lc_df$value))

# Order sites by decreasing amount of development
site_order = sites_lc_df %>% group_by(site) %>%
  summarise(developed_n = sum(grepl("^Developed", class)), .groups = "drop") %>%
  arrange(desc(developed_n)) %>% pull(site)
sites_lc_df = sites_lc_df %>%
  mutate(site = factor(site, levels = site_order))

ggplot(sites_lc_df, aes(x = x, y = y, fill = factor(class))) + geom_tile() +
  scale_fill_manual(values = setNames(present_levels$color, present_levels$class), name = "Class") +
  facet_wrap(~site, scales = "free") +
  theme_minimal() + theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = paste0("NLCD land cover configuration (", buffer_radius, " m buffer)"), x = "", y = "")

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

# TODO: assess adequacy of sampling effort to detect total species richness via sample size-based rarefaction and extrapolation sampling curves for species richness (iNEXT)? Quantify sample completeness to determine the estimated proportion of species detected from the predicted species pool by plotting sample coverage with respect to the number of sampling units.
if (FALSE) {
  library(iNEXT)
  
  zero_detections <- detections$long %>%
    group_by(common_name) %>%
    summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
    filter(total_detections == 0) %>%
    pull(common_name)
  
  detections_list <- detections$wide %>%
    group_split(site_id) %>%                                 # Split by site
    set_names(unique(detections$wide$site_id)) %>%            # Name list elements by site_id
    map(~ {
      df <- .x %>%
        select(-site_id) %>%                                  # Remove site_id column
        tibble::column_to_rownames("common_name")             # Set rownames as common_name
      return(df)
    })
  detections_incidence <- detections_list %>%
    map(~ as.data.frame((.x > 0) * 1))
  
  out.raw <- iNEXT(detections_incidence, q = 0, datatype="incidence_raw", endpoint=100)
  out.raw
  ggiNEXT(out.raw, type = 1) # sample size-based rarefaction/extrapolation curve
  ggiNEXT(out.raw, type = 2) # sample completeness curve
  ggiNEXT(out.raw, type = 3) # coverage-based rarefaction/extrapolation curve
  
}

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

# Load Stevan riparian dependent/obligate list
riparian_dep = read_csv("data/processed/Species_Habitat_List.csv") %>% clean_names() %>%
  rename(common_name = species) %>% mutate(common_name = tolower(common_name))
rip_deps = riparian_dep %>% filter(riparian_dependent_breeding == "Yes" | riparian_obligate_breeding == "Yes") %>% pull(common_name)

setdiff(rip_deps, species_names$common_name)

# Exclude certain sites from analysis ------------------------------------------------

# Inspect total detections per site
# total_detections_by_site = detections$long %>% group_by(site_id) %>%
#   summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
#   arrange(desc(total_detections))
# ggplot(total_detections_by_site, aes(x = reorder(site_id, total_detections), y = total_detections)) + geom_col()

# Exclude sites 257 and 259 that are dominated by agriculture
sites_to_exclude = c("257", "259")
message("Excluding agricultural site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data = site_data %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

# Exclude sites with incomplete surveys
sites_to_exclude = setdiff(unique(site_data$site_id), unique(detections$long$site_id))
sites_to_exclude = c(sites_to_exclude, "150") # ARU at site 150 destroyed by water damage
message("Excluding incomplete survey site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data = site_data %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

message("Retaining ", length(unique(site_data$site_id)), " sites")

# Visualize joint data ----------------------------------------------------------------
message("Visualizing joint data")

presence_absence = detections$long %>% group_by(site_id, common_name) %>%
  summarise(presence = if_else(sum(n_detections, na.rm = TRUE) > 0, 1, 0), .groups = "drop")

# Calculate richness of different groups per site
richness = presence_absence %>% group_by(site_id) %>% summarise(richness = sum(presence))
richness_invert_et = presence_absence %>% filter(common_name %in% invertivores_eltontraits) %>%
  group_by(site_id) %>% summarise(richness_invert_et = sum(presence))
richness_invert_avo = presence_absence %>% filter(common_name %in% invertivores_avonet) %>%
  group_by(site_id) %>% summarise(richness_invert_avo = sum(presence))
richness_ripdep = presence_absence %>% filter(common_name %in% rip_deps) %>%
  group_by(site_id) %>% summarise(richness_ripdep = sum(presence))

## NMDS visualization
library(vegan)
nmds = metaMDS(presence_absence %>%
          pivot_wider(names_from = common_name, values_from = presence, values_fill = 0) %>%
          column_to_rownames("site_id"), distance = "bray", k = 3, trymax = 100)
nmds$stress # stress > 0.2 suggests weak fit
nmds_scores = as.data.frame(scores(nmds, display = "sites"))
nmds_scores$site_id = rownames(nmds_scores)
nmds_plot_data = nmds_scores %>% left_join(site_metadata, by = "site_id")
species_scores = as.data.frame(scores(nmds, display = "species"))
species_scores$species = rownames(species_scores)

ggplot(nmds_plot_data %>% left_join(site_data, by = "site_id"), aes(NMDS1, NMDS2, color = bibi)) +
  geom_text(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), color = "gray", size = 3, vjust = -0.5, check_overlap = TRUE) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text(aes(label = site_id), vjust = -0.5, size = 3, color = "black") +
  scale_color_continuous(type = "viridis") +
  labs(color = "Impervious %")

ggplot(species_scores, aes(x = 0, y = 0)) +
  geom_segment(aes(xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "darkgray", alpha = 0.8) +
  geom_text(aes(x = NMDS1, y = NMDS2, label = species),
            color = "black", size = 3, vjust = -0.5, check_overlap = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
  scale_x_continuous(limits = c(-2, 2)) +
  scale_y_continuous(limits = c(-2, 2)) +
  labs(title = "Species vectors")

# Join with site data
site_data = left_join(site_data, richness, by = "site_id")
site_data = left_join(site_data, richness_invert_et, by = "site_id")
site_data = left_join(site_data, richness_invert_avo, by = "site_id")
site_data = left_join(site_data, richness_ripdep, by = "site_id")

# Richness across AVONET guilds per site
richness_by_trophic_niche = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, trophic_niche) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_primary_lifestyle = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, primary_lifestyle) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_ripdep = presence_absence %>% mutate(ripdep = common_name %in% c(rip_deps)) %>% group_by(site_id, ripdep) %>% summarize(count = sum(presence))

ggplot(richness_by_trophic_niche, aes(x = count, y = reorder(site_id, count), fill = trophic_niche)) + geom_col()
ggplot(richness_by_primary_lifestyle, aes(x = count, y = reorder(site_id, count), fill = primary_lifestyle)) + geom_col()
ggplot(richness_by_ripdep, aes(x = count, y = reorder(site_id, count), fill = ripdep)) + geom_col()

# Richness as a function of different predictors
ggplot(site_data, aes(x = rast_usfs_canopycover_sum_proportion, y = richness_invert_avo)) + geom_point() + geom_smooth() +
  labs(title = "Canopy cover") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data, aes(x = rast_nlcd_impervious_sum_proportion, y = richness_invert_avo)) + geom_point() + geom_smooth() +
  labs(title = "Local imperviousness") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data, aes(x = basin_impervious, y = richness_invert_avo)) + geom_point() + geom_smooth() +
  labs(title = "Basin imperviousness") + geom_text_repel(aes(label = site_id))

ggplot(site_data, aes(x = bibi, y = richness_invert_avo)) + geom_point() + geom_smooth() +
  labs(title = "BIBI") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data, aes(x = (density_roads_paved), y = richness_invert_avo)) + geom_point() + geom_smooth() +
  labs(title = "Density paved roads") + geom_text_repel(aes(label = site_id)) 

mapview(site_data %>% select(site_id, richness_invert_avo), zcol = "richness_invert_avo")

# Land cover visualization ordered by richness
site_order = site_data %>% arrange(desc(richness_invert_avo)) %>% pull(site_id)

sites_lc_df = bind_rows(
  lapply(names(geospatial_site_data), function(site_id) {
    # rast_nlcd_landcover (interaction zone only) OR rast_nlcd_landcover_site_buffer
    df = as.data.frame(geospatial_site_data[[site_id]]$rasters$rast_nlcd_landcover, xy = TRUE)
    value_col = names(df)[3]
    df = df %>% rename(value = all_of(value_col))  # rename for ggplot
    df$site <- site_id
    df$class = nlcd_metadata[match(df$value, nlcd_metadata$id), "class"]
    return(df)
  })
) %>% mutate(site = factor(site, levels = site_order))
present_levels = nlcd_metadata %>% filter(id %in% unique(sites_lc_df$value))

p = ggplot(sites_lc_df, aes(x = x, y = y, fill = factor(class))) + geom_tile() +
  scale_fill_manual(values = setNames(present_levels$color, present_levels$class), name = "Class") +
  facet_wrap(~site, scales = "free") +
  theme_minimal() + theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  labs(title = paste0("Land cover by invert richness (descending)"), x = "", y = ""); print(p)

sites_imp_df = bind_rows(
  lapply(seq_along(geospatial_site_data), function(i) {
    raster_to_df(geospatial_site_data[[i]]$rasters$rast_nlcd_impervious, names(geospatial_site_data)[i])
  })
) %>% mutate(site = factor(site, levels = site_order))

p = ggplot(sites_imp_df, aes(x = x, y = y, fill = class)) + geom_tile() +
  facet_wrap(~ site, scales = "free") +
  scale_fill_viridis_c(option = "inferno") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank()) +
  labs(title = paste0("Imperviousness by invert richness (descending)"), x = "", y = ""); print(p)
  
# Richness of guilds as a function of...
ggplot(left_join(richness_by_trophic_niche, site_data, by = "site_id"),
       aes(x = bibi, y = count, color = trophic_niche, fill = trophic_niche)) +
  geom_point() + geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data, by = "site_id"),
       aes(x = bibi, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_ripdep, site_data, by = "site_id"),
       aes(x = bibi, y = count, color = ripdep, fill = ripdep)) +
  geom_point() + geom_smooth(aes(group = ripdep), method = "lm", se = FALSE)

ggplot(left_join(richness_by_ripdep, site_data, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = ripdep, fill = ripdep)) +
  geom_point() + geom_smooth(aes(group = ripdep), method = "lm", se = FALSE)

# Wilson's Warbler presence/absence as a function of BIBI
ggplot(left_join(presence_absence %>% filter(common_name == "wilson's warbler"), site_data, by = "site_id"),
       aes(x = bibi, y = presence)) +
  geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

# Cache site data --------------------------------------------------------------------

if (!dir.exists(out_cache_dir)) dir.create(out_cache_dir, recursive = TRUE)
out_filepath = file.path(out_cache_dir, paste0(
  "site_data_", buffer_radius, "m_",
  detections$threshold_classifier_score, "t_", 
  detections$threshold_detected_days, "d.rds"))
saveRDS(site_data, out_filepath)
message(crayon::green("Cached", out_filepath))

stop("Finished!")

# Structural equation modeling -------------------------------------------------------

candidate_vars = c(
  ## Predictors
  "bibi"                 = "bibi",
  "site"                 = "site_id",
  "basin"                = "basin_name",
  # "elevation"            = "elevation",
  "imperviousness_basin" = "basin_impervious",
  "imperviousness_local" = "rast_nlcd_impervious_sum_proportion",
  "abund_dev_varint"     = "nlcd_developed_variable_intensity",
  "abund_dev_opensp"     = "nlcd_developed_open_space",
  "abund_forest"         = "nlcd_forest",
  "abund_wetland"        = "nlcd_wetlands",
  "abund_openwater"      = "nlcd_open_water",
  "canopy_usfs"          = "rast_usfs_canopycover_sum_proportion",
  "canopy_gedi"          = "rast_gedi_cover_sum",
  "height_landfire"      = "rast_landfire_treeheight_mean",
  "height_gedi"          = "rast_gedi_height_mean",
  "height_cv_landfire"   = "rast_landfire_treeheight_cv",
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
data = as.data.frame(site_data) %>% select(all_of(candidate_vars))

# Calculate pairwise collinearity among predictors
pairwise_collinearity = function(vars, threshold = 0.7) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

pairwise_collinearity(data %>% select(
  bibi,
  imperviousness_basin, imperviousness_local,
  canopy_usfs, abund_forest,
  height_gedi,
  fhd, height_cv_gedi, height_cv_gedi_mask,
  pavd20m, pavd5to10m,
  # elevation,
  edge_density, aggregation,
  stand_age_mean, stand_age_cv, stand_ba, stand_denall, stand_dencon, stand_denhw, stand_qmd, tree_richness
))

# TODO: VIF analyses for multicollinearity, e.g. sort(car::vif(m_bibi))

# Create and inspect structural equation model claims and fit
sem = psem(
  # Fit component regressions
  lm(bibi ~ imperviousness_local + canopy_usfs,
     data),
  glm(richness_invert_avo ~ bibi + imperviousness_local + canopy_usfs + height_cv_gedi_mask,
      data, family = poisson(link = "log"))
); plot(sem); print(summary(sem))
# - A significant independence claim from a test of directed separation suggests that the path is missing
# or misspecified.
# - A significant global Fisher’s C p-value (< 0.05) suggests that the modeled structure is statistically
# significantly different than the structure implied by the data, and that alternative pathways or causal
# links with missing variables warrant further exploration

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
# - type of stream (perennial, intermittent)?
# - year?

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
