#########################################################################################
## Clean and cache geospatial dependencies for the study area
in_data_geospatial = "/Volumes/gioj_work/riparian-avian-macroinvert-monitoring-kc/data/raw/geospatial"
out_dir = "data/cache/preprocess_geospatial_data"
out_filepath = "data/cache/preprocess_geospatial_data/geospatial_data.rds"
#########################################################################################

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

standard_crs_code = "EPSG:32610"

#########################################################################################
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
study_area = counties(state = "WA", year = 2024, cb = TRUE) %>% filter(NAME == "King") %>% st_transform(crs = standard_crs_code)

mapview(study_area, alpha.regions = 0, lwd = 2) +
  mapview(sites_aru, zcol = "dist_m", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

############################################################
# Load land cover and impervious surface data
message("Loading USGS NLCD land cover data")

lc_raw  = rast(paste0(in_data_geospatial, "/NLCD/Annual_NLCD_LndCov_2023_CU_C1V0.tif"))

template = project(vect(study_area), crs(lc_raw))
rast_landcover  = mask(crop(lc_raw, template), template)

# Factor land cover data
# https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
# rast_landcover = as.factor(rast_landcover)
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
# levels(rast_landcover) <- nlcd_levels

rast_landcover[] = as.integer(factor(rast_landcover[], levels = nlcd_levels$id))
levels(rast_landcover) = data.frame(id = 1:nrow(nlcd_levels), class = nlcd_levels$class, color = nlcd_levels$color)
levels(rast_landcover)[[1]] # verify

plot(rast_landcover)

############################################################
# Load impervious surface percentage
message("Loading USGS NLCD impervious surface data")

imp_raw = rast(paste0(in_data_geospatial, "/NLCD/Annual_NLCD_FctImp_2023_CU_C1V0.tif"))

template = project(vect(study_area), crs(imp_raw))
rast_impervious = mask(crop(imp_raw, template), template)

plot(rast_impervious)

############################################################
# Load tree canopy cover
message("Loading USFS tree canopy cover data")
tcc_raw = rast(paste0(in_data_geospatial, "/Forest Service Science TCC/science_tcc_conus_wgs84_v2023-5_20230101_20231231.tif"))

template = project(vect(study_area), crs(tcc_raw))
rast_canopycover = mask(crop(tcc_raw, template), template)

values(rast_canopycover) <- as.numeric(values(rast_canopycover))
rast_canopycover[rast_canopycover > 100] = NA
rast_canopycover[is.nan(rast_canopycover)] = NA

plot(rast_canopycover)

############################################################
# Load landfire geospatial data
message("Loading USDA/DOI Landfire geospatial data")

# Vegetation cover
# "Represents the vertically projected percent cover of the live canopy for a 30-m cell. rast_vegcover is produced separately for tree, shrub, and herbaceous lifeforms. Training data depicting percentages of canopy cover are obtained from plot-level ground-based visual assessments and lidar observations. These are combined with Landsat imagery (from multiple seasons), to inform models built independently for each lifeform. Tree, shrub, and herbaceous lifeforms each have a potential range from 10% to 100% (cover values less than 10% are binned into the 10% value). The three independent lifeform datasets are merged into a single product based on the dominant lifeform of each pixel. The rast_vegcover product is then reconciled through QA/QC measures to ensure lifeform is synchronized with Existing Vegetation Height (rast_vegheight)."
evc_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVC_250_CONUS/LC24_EVC_250.tif"))
template = project(vect(study_area), crs(evc_raw))
rast_vegcover = mask(crop(evc_raw, template), template)

evc_values = as.data.frame(values(rast_vegcover)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evc_table  = as.data.frame(cats(rast_vegcover)[[1]]) %>% select(VALUE, CLASSNAMES)
(evc_composition = evc_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evc_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_vegcover)

# Vegetation height
# "Represents the average height of the dominant vegetation for a 30-m cell. rast_vegheight is produced separately for tree, shrub, and herbaceous lifeforms using training data depicting the weighted average height by species cover and Existing Vegetation Type (rast_vegtype) lifeform. Decision tree models using field reference data, lidar, and Landsat are developed separately for each lifeform, then lifeform specific height class layers are merged along with land cover into a single rast_vegheight product based on the dominant lifeform of each pixel. rast_vegheight ranges are continuous for the herbaceous lifeform category ranging from 0.1 to 1 meter with decimeter increments, 0.1 to 3 meters for shrub lifeform, and 1 to 99 meters for tree lifeform. If the height values of each lifeform exceed the continuous value range, they are binned into the appropriate maximum height class. rast_vegheight is then reconciled through QA/QC measures to ensure lifeform is synchronized with Existing Vegetation Cover (rast_vegcover)."
evh_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVH_250_CONUS/LC24_EVH_250.tif"))
template = project(vect(study_area), crs(evh_raw))
rast_vegheight = mask(crop(evh_raw, template), template)

evh_table = as.data.frame(cats(rast_vegheight)[[1]])
(evh_table$CLASSNAMES[evh_table$VALUE %in% unique(values(rast_vegheight))])

evh_values = as.data.frame(values(rast_vegheight)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evh_table  = as.data.frame(cats(rast_vegheight)[[1]]) %>% select(VALUE, CLASSNAMES)
(evh_composition = evh_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evh_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_vegheight)

# Convert vegetation height data from categorical to continuous

tree_cats = cats(rast_vegheight)[[1]] %>% mutate(HEIGHT = NA)
tree_rows = grepl("^Tree Height", tree_cats$CLASSNAMES)
tree_cats$HEIGHT[tree_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", tree_cats$CLASSNAMES[tree_rows]))
rast_treeheight = classify(rast_vegheight, rcl = as.matrix(data.frame(from = tree_cats$VALUE, to = tree_cats$HEIGHT)), others = NA)
plot(rast_treeheight)

shrub_cats = cats(rast_vegheight)[[1]] %>% mutate(HEIGHT = NA)
shrub_rows = grepl("^Shrub Height", shrub_cats$CLASSNAMES)
shrub_cats$HEIGHT[shrub_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", shrub_cats$CLASSNAMES[shrub_rows]))
rast_shrubheight = classify(rast_vegheight, rcl = as.matrix(data.frame(from = shrub_cats$VALUE, to = shrub_cats$HEIGHT)), others = NA)
plot(rast_shrubheight)

herb_cats = cats(rast_vegheight)[[1]] %>% mutate(HEIGHT = NA)
herb_rows = grepl("^Herb Height", herb_cats$CLASSNAMES)
herb_cats$HEIGHT[herb_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", herb_cats$CLASSNAMES[herb_rows]))
rast_herbheight = classify(rast_vegheight, rcl = as.matrix(data.frame(from = herb_cats$VALUE, to = herb_cats$HEIGHT)), others = NA)
plot(rast_herbheight)


# Vegetation type
# "Represents the current distribution of the terrestrial ecological systems classification developed by NatureServe for the western hemisphere. In this context, a terrestrial ecological system is defined as a group of plant community types that tend to co-occur within landscapes with similar ecological processes, substrates, and/or environmental gradients. See the rast_vegtype product page (https://landfire.gov/vegetation/rast_vegtype) for more information about ecological systems and NVC classifications."
evt_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVT_250_CONUS/LC24_EVT_250.tif"))
template = project(vect(study_area), crs(evt_raw))
rast_vegtype = mask(crop(evt_raw, template), template)

evt_values = as.data.frame(values(rast_vegtype)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evt_table  = as.data.frame(cats(rast_vegtype)[[1]]) %>%  select(VALUE, EVT_NAME)
(evt_composition = evt_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evt_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_vegtype)

# Successional class
# "Categorizes current vegetation composition and structure into up to five successional classes, with successional classes defined in the appropriate Biophysical Settings (BpS) Model. There are two additional categories for uncharacteristic species (exotic or invasive vegetation), and uncharacteristic native vegetation cover, structure, or composition. Current successional classes and their historical reference conditions are compared to assess departure of vegetation characteristics. The classification schemes used to produce BpS and SClass may vary slightly between adjacent map zones, and reference conditions may be simulated independently in different map zones for the same BpS.
sc_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_SClass_250_CONUS/LC24_SCla_250.tif"))
template = project(vect(study_area), crs(sc_raw))
rast_sclass = mask(crop(sc_raw, template), template)

sc_table = as.data.frame(cats(rast_sclass)[[1]])
(sc_table$DESCRIPTIO[sc_table$VALUE %in% unique(values(rast_sclass))])

plot(rast_sclass)

############################################################
# Load landfire geospatial data
message("Loading King County DNRP hydrological data")

# basins = st_read(paste0(in_data_geospatial, "/King County DNRP/KC_basins.shp")

# areas = st_read(paste0(in_data_geospatial, "/NHD_H_1711_HU4_Shape/Shape/NHDArea.shp")

# flowlines = st_read(paste0(in_data_geospatial, "/NHD_H_1711_HU4_Shape/Shape/NHDFlowline.shp")

# "The riparian delineations we have are for both fixed widths as well as “functionally dynamic” buffers that may be more ecologically relevant."
# "As Abood and others (2012) describe, the variable-width, or “dynamic,” buffer is likely more ecologically relevant because it accounts for factors that affect how a stream interacts and is influenced by the riparian zone. A recent literature review of riparian buffers done by King County (2019b) highlights the range of widths needed to maintain various functions (e.g., erosion control, shade). The review did not evaluate this concept of dynamic buffers, but the findings illustrate that riparian functions are maintained at different widths depending on a range of factors including but not limited to the size of stream, soil composition, vegetation type and age, etc.
sf_ripfb = st_read(paste0(in_data_geospatial, "/King County DNRP/RiparianBuffer_basin.shp"), quiet = TRUE)
sf_ripfb = sf_ripfb %>% filter(lengths(st_intersects(geometry, st_as_sf(study_area) %>% st_transform(st_crs(sf_ripfb)))) > 0)

# "The pattern, quality, and connectivity of riparian areas can also be really interesting and I’ve wondered how that may affect birds and other wildlife. For instance, in Seattle, the streams can have pretty decent but very narrow riparian areas because they are in canyons. In other suburban watersheds, the riparian areas can be vegetated but lack the tall evergreens or complex structure - they may be wider and even “greener” but not as functional? Maybe birds with small territories can fare OK in urban riparian areas?"

############################################################
# Load NASA GEDI-Fusion forest structure data
#
# https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-gedi-fusion-structure-2236-1
message("Loading NASA GEDI-Fusion forest structure data")

# "Foliage height diversity (a unitless index)"
fhd_raw = rast(paste0(in_data_geospatial, "/NASA/gedifusion_fhd_2020.tif"))
template = project(vect(study_area), crs(fhd_raw))
rast_fhd = mask(crop(fhd_raw, template), template)

plot(rast_fhd)

# "Fractional canopy cover (proportion)"
cover_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_cover_2020.tif"))
template   = project(vect(study_area), crs(cover_raw))
rast_cover = mask(crop(cover_raw, template), template)

plot(rast_cover)

# "Corresponds to the height at which 98% of the waveform energy is captured - comparable to a canopy height measure"
height_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_rh98_2020.tif"))
template   = project(vect(study_area), crs(height_raw))
rast_height = mask(crop(height_raw, template), template)

plot(rast_height)

# "Plant area vegetation density (PAVD); the proportion of vegetation within the 5-10 m stratum above ground surface (PAVD 5-10 m)""
pavd55o10m_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_pavd5to10m_2020.tif"))
template   = project(vect(study_area), crs(pavd55o10m_raw))
rast_pavd55o10m = mask(crop(pavd55o10m_raw, template), template)

plot(rast_pavd55o10m)

# "The proportion of vegetation density (PAVD) greater than 20 m above ground surface, chosen to represent the presence of a mature upper canopy within different forest types"

pavd20m_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_pavd20m_2020.tif"))
template   = project(vect(study_area), crs(pavd20m_raw))
rast_pavd20m = mask(crop(pavd20m_raw, template), template)

plot(rast_pavd20m)

############################################################
# Project all data to same EPSG:32610 coordinate reference system
message("Projecting data to ", standard_crs_code)
study_area       = project(vect(study_area), standard_crs_code)

rast_data = list(
  rast_landcover   = project(rast_landcover,   standard_crs_code),
  rast_impervious  = project(rast_impervious,  standard_crs_code),
  rast_canopycover = project(rast_canopycover, standard_crs_code),
  rast_vegcover    = project(rast_vegcover,    standard_crs_code),
  rast_vegheight   = project(rast_vegheight,   standard_crs_code),
  rast_treeheight  = project(rast_treeheight,  standard_crs_code),
  rast_shrubheight = project(rast_shrubheight, standard_crs_code),
  rast_herbheight  = project(rast_herbheight,  standard_crs_code),
  rast_vegtype     = project(rast_vegtype,     standard_crs_code),
  rast_sclass      = project(rast_sclass,      standard_crs_code),
  rast_fhd         = project(rast_fhd,         standard_crs_code),
  rast_cover       = project(rast_cover,       standard_crs_code),
  rast_height      = project(rast_height,      standard_crs_code),
  rast_pavd55o10m  = project(rast_pavd55o10m,  standard_crs_code),
  rast_pavd20m     = project(rast_pavd20m,     standard_crs_code)
)

sf_ripfb = sf_ripfb %>% st_transform(crs = standard_crs_code)

# Cache data
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
for (n in names(rast_data)) {
  out_filepath = file.path(out_dir, paste0(n, ".tif"))
  writeRaster(rast_data[[n]], out_filepath, overwrite = TRUE)
  message("Cached ", out_filepath)
}

out_filepath = file.path(out_dir, "sf_ripfb.gpkg")
st_write(sf_ripfb, out_filepath, append = FALSE)
message("Cached ", out_filepath)

############################################################################################################
# Load ARU and PSSB site locations and calculate all site covariates from land cover data

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

primary_lifestyle_summary <- site_species_long %>%
  left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, primary_lifestyle) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = primary_lifestyle, values_from = count, values_fill = 0)

trophic_niche_summary <- site_species_long %>%
  left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, trophic_niche) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = trophic_niche, values_from = count, values_fill = 0)

trophic_niche_long <- trophic_niche_summary %>%
  pivot_longer(cols = -site_id, names_to = "trophic_niche", values_to = "count")

primary_lifestyle_long <- primary_lifestyle_summary %>%
  pivot_longer(cols = -site_id, names_to = "primary_lifestyle", values_to = "count")

# Plot trends in trophic niche richness
ggplot(left_join(trophic_niche_long, site_data %>% select(site_id, tcc_sum), by = "site_id"), aes(x = tcc_sum, y = count, color = trophic_niche, group = trophic_niche, fill = trophic_niche)) +
  geom_point() +
  geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

ggplot(left_join(trophic_niche_long, site_data %>% select(site_id, imp_sum), by = "site_id"), aes(x = imp_sum, y = count, color = trophic_niche, group = trophic_niche, fill = trophic_niche)) +
  geom_point() +
  geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

# Plot trends in trophic level richness
ggplot(left_join(primary_lifestyle_long, site_data %>% select(site_id, tcc_sum), by = "site_id"), aes(x = tcc_sum, y = count, color = primary_lifestyle, group = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() +
  geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(primary_lifestyle_long, site_data %>% select(site_id, imp_sum), by = "site_id"), aes(x = imp_sum, y = count, color = primary_lifestyle, group = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() +
  geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

# Calculate richness of different groups

richness = site_species_matrix %>% mutate(richness = rowSums(across(-site_id))) %>% select(site_id, richness)
richness_insectivore = site_insectivores_matrix %>% mutate(richness_insectivore = rowSums(across(-site_id))) %>% select(site_id, richness_insectivore)

# Join with other site data
site_data_bird = right_join(site_data, richness, by = "site_id")
site_data_bird = right_join(site_data_bird, richness_insectivore, by = "site_id")

ggplot(site_data_bird, aes(x = bibi, y = richness)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(site_data_bird, aes(x = bibi, y = richness_insectivore)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(site_data_bird, aes(x = tcc_sum, y = richness_insectivore)) +
  geom_point() + geom_smooth(method = "lm")

ggplot(site_data_bird, aes(x = treeheight_mean, y = richness_insectivore)) +
  geom_point() + geom_smooth(method = "lm")

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
                   "treeheight_mean", "treeheight_sd",
                   "nlcd_shrub_scrub", "nlcd_wetlands")

data = as.data.frame(site_data_bird) %>% janitor::clean_names() %>% select(
  all_of(candidate_vars)
)

# Explore pairwise collinearity
pairwise_collinearity = function(vars, threshold = 0.7) {
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
model_bibi = lm(bibi ~ imp_sum + tcc_sum,
                data)

model_alpha_total = glm(richness_insectivore ~ bibi + imp_sum + tcc_sum + treeheight_mean,
                        data, family = poisson(link = "log"))

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

# ALTERNATIVE
m1 = lm(bibi ~ imp_sum + tcc_sum, data)
m2 = lm(tcc_sum ~ treeheight_mean + imp_sum, data)
m3 = glm(richness_insectivore ~ bibi + imp_sum + tcc_sum + treeheight_mean, data, family = poisson(link = "log"))
sem_alt = psem(
  m1,
  m2,
  m3
)
plot(sem_alt)
print(summary(sem_alt))

