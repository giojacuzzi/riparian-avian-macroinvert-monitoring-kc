#########################################################################################
# Clean and cache geospatial dependencies for the study area
#
# Inputs:
in_data_geospatial = "/Volumes/gioj_work/riparian-avian-macroinvert-monitoring-kc/data/raw/geospatial"
in_path_nlcd_metadata = "data/raw/nlcd_metadata.csv"
# Outputs:
out_dir = "data/cache/2_preprocess_geospatial_data"
#########################################################################################

# Load required packages (automatically install any missing)
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

crs_standard = "EPSG:32610"

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
study_area = counties(state = "WA", year = 2024, cb = TRUE) %>% filter(NAME == "King") %>% st_transform(crs = crs_standard)

mapview(study_area, alpha.regions = 0, lwd = 2) +
  mapview(sites_aru, zcol = "dist_m", layer.name = "ARU") +
  mapview(sites_pssb, col.region = "blue", layer.name = "PSSB")

#########################################################################################
# Load land cover and impervious surface data
message("Loading USGS NLCD land cover data")

lc_raw  = rast(paste0(in_data_geospatial, "/NLCD/Annual_NLCD_LndCov_2023_CU_C1V0.tif"))

template = project(vect(study_area), crs(lc_raw))
rast_nlcd_landcover  = mask(crop(lc_raw, template), template)

# Factor land cover data
# https://www.mrlc.gov/data/legends/national-land-cover-database-class-legend-and-description
nlcd_metadata = read.csv(in_path_nlcd_metadata)

rast_nlcd_landcover[] = as.integer(factor(rast_nlcd_landcover[], levels = nlcd_metadata$id))
levels(rast_nlcd_landcover) = data.frame(id = 1:nrow(nlcd_metadata), class = nlcd_metadata$class, color = nlcd_metadata$color)
levels(rast_nlcd_landcover)[[1]] # verify

plot(rast_nlcd_landcover)

#########################################################################################
# Load impervious surface percentage
message("Loading USGS NLCD impervious surface data")

imp_raw = rast(paste0(in_data_geospatial, "/NLCD/Annual_NLCD_FctImp_2023_CU_C1V0.tif"))

template = project(vect(study_area), crs(imp_raw))
rast_nlcd_impervious = mask(crop(imp_raw, template), template)

plot(rast_nlcd_impervious)

#########################################################################################
# Load tree canopy cover
message("Loading USFS tree canopy cover data")
tcc_raw = rast(paste0(in_data_geospatial, "/Forest Service Science TCC/science_tcc_conus_wgs84_v2023-5_20230101_20231231.tif"))

template = project(vect(study_area), crs(tcc_raw))
rast_usfs_canopycover = mask(crop(tcc_raw, template), template)

values(rast_usfs_canopycover) <- as.numeric(values(rast_usfs_canopycover))
rast_usfs_canopycover[rast_usfs_canopycover > 100] = NA
rast_usfs_canopycover[is.nan(rast_usfs_canopycover)] = NA

plot(rast_usfs_canopycover)

#########################################################################################
# Load landfire geospatial data
message("Loading USDA/DOI Landfire geospatial data")

# Vegetation cover
# "Represents the vertically projected percent cover of the live canopy for a 30-m cell. rast_landfire_vegcover is produced separately for tree, shrub, and herbaceous lifeforms. Training data depicting percentages of canopy cover are obtained from plot-level ground-based visual assessments and lidar observations. These are combined with Landsat imagery (from multiple seasons), to inform models built independently for each lifeform. Tree, shrub, and herbaceous lifeforms each have a potential range from 10% to 100% (cover values less than 10% are binned into the 10% value). The three independent lifeform datasets are merged into a single product based on the dominant lifeform of each pixel. The rast_landfire_vegcover product is then reconciled through QA/QC measures to ensure lifeform is synchronized with Existing Vegetation Height (rast_landfire_vegheight)."
evc_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVC_250_CONUS/LC24_EVC_250.tif"))
template = project(vect(study_area), crs(evc_raw))
rast_landfire_vegcover = mask(crop(evc_raw, template), template)

evc_values = as.data.frame(values(rast_landfire_vegcover)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evc_table  = as.data.frame(cats(rast_landfire_vegcover)[[1]]) %>% select(VALUE, CLASSNAMES)
(evc_composition = evc_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evc_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_landfire_vegcover)

# Vegetation height
# "Represents the average height of the dominant vegetation for a 30-m cell. rast_landfire_vegheight is produced separately for tree, shrub, and herbaceous lifeforms using training data depicting the weighted average height by species cover and Existing Vegetation Type (rast_landfire_vegtype) lifeform. Decision tree models using field reference data, lidar, and Landsat are developed separately for each lifeform, then lifeform specific height class layers are merged along with land cover into a single rast_landfire_vegheight product based on the dominant lifeform of each pixel. rast_landfire_vegheight ranges are continuous for the herbaceous lifeform category ranging from 0.1 to 1 meter with decimeter increments, 0.1 to 3 meters for shrub lifeform, and 1 to 99 meters for tree lifeform. If the height values of each lifeform exceed the continuous value range, they are binned into the appropriate maximum height class. rast_landfire_vegheight is then reconciled through QA/QC measures to ensure lifeform is synchronized with Existing Vegetation Cover (rast_landfire_vegcover)."
evh_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVH_250_CONUS/LC24_EVH_250.tif"))
template = project(vect(study_area), crs(evh_raw))
rast_landfire_vegheight = mask(crop(evh_raw, template), template)

evh_table = as.data.frame(cats(rast_landfire_vegheight)[[1]])
(evh_table$CLASSNAMES[evh_table$VALUE %in% unique(values(rast_landfire_vegheight))])

evh_values = as.data.frame(values(rast_landfire_vegheight)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evh_table  = as.data.frame(cats(rast_landfire_vegheight)[[1]]) %>% select(VALUE, CLASSNAMES)
(evh_composition = evh_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evh_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_landfire_vegheight)

# Convert vegetation height data from categorical to continuous

tree_cats = cats(rast_landfire_vegheight)[[1]] %>% mutate(HEIGHT = NA)
tree_rows = grepl("^Tree Height", tree_cats$CLASSNAMES)
tree_cats$HEIGHT[tree_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", tree_cats$CLASSNAMES[tree_rows]))
rast_landfire_treeheight = classify(rast_landfire_vegheight, rcl = as.matrix(data.frame(from = tree_cats$VALUE, to = tree_cats$HEIGHT)), others = NA)
plot(rast_landfire_treeheight)

shrub_cats = cats(rast_landfire_vegheight)[[1]] %>% mutate(HEIGHT = NA)
shrub_rows = grepl("^Shrub Height", shrub_cats$CLASSNAMES)
shrub_cats$HEIGHT[shrub_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", shrub_cats$CLASSNAMES[shrub_rows]))
rast_landfire_shrubheight = classify(rast_landfire_vegheight, rcl = as.matrix(data.frame(from = shrub_cats$VALUE, to = shrub_cats$HEIGHT)), others = NA)
plot(rast_landfire_shrubheight)

herb_cats = cats(rast_landfire_vegheight)[[1]] %>% mutate(HEIGHT = NA)
herb_rows = grepl("^Herb Height", herb_cats$CLASSNAMES)
herb_cats$HEIGHT[herb_rows] = as.numeric(sub(".*= ([0-9\\.]+).*", "\\1", herb_cats$CLASSNAMES[herb_rows]))
rast_landfire_herbheight = classify(rast_landfire_vegheight, rcl = as.matrix(data.frame(from = herb_cats$VALUE, to = herb_cats$HEIGHT)), others = NA)
plot(rast_landfire_herbheight)

# Vegetation type
# "Represents the current distribution of the terrestrial ecological systems classification developed by NatureServe for the western hemisphere. In this context, a terrestrial ecological system is defined as a group of plant community types that tend to co-occur within landscapes with similar ecological processes, substrates, and/or environmental gradients. See the rast_landfire_vegtype product page (https://landfire.gov/vegetation/rast_landfire_vegtype) for more information about ecological systems and NVC classifications."
evt_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_EVT_250_CONUS/LC24_EVT_250.tif"))
template = project(vect(study_area), crs(evt_raw))
rast_landfire_vegtype = mask(crop(evt_raw, template), template)

evt_values = as.data.frame(values(rast_landfire_vegtype)) %>% rename(VALUE = 1) %>% filter(!is.na(VALUE))
evt_table  = as.data.frame(cats(rast_landfire_vegtype)[[1]]) %>%  select(VALUE, EVT_NAME)
(evt_composition = evt_values %>%
    group_by(VALUE) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(pcnt = (count / sum(count))) %>%
    left_join(evt_table, by = "VALUE") %>%
    arrange(desc(pcnt)))

plot(rast_landfire_vegtype)

# Successional class
# "Categorizes current vegetation composition and structure into up to five successional classes, with successional classes defined in the appropriate Biophysical Settings (BpS) Model. There are two additional categories for uncharacteristic species (exotic or invasive vegetation), and uncharacteristic native vegetation cover, structure, or composition. Current successional classes and their historical reference conditions are compared to assess departure of vegetation characteristics. The classification schemes used to produce BpS and SClass may vary slightly between adjacent map zones, and reference conditions may be simulated independently in different map zones for the same BpS.
sc_raw = rast(paste0(in_data_geospatial, "/Landfire/LF2024_SClass_250_CONUS/LC24_SCla_250.tif"))
template = project(vect(study_area), crs(sc_raw))
rast_landfire_sclass = mask(crop(sc_raw, template), template)

sc_table = as.data.frame(cats(rast_landfire_sclass)[[1]])
(sc_table$DESCRIPTIO[sc_table$VALUE %in% unique(values(rast_landfire_sclass))])

plot(rast_landfire_sclass)

#########################################################################################
# Load King County DNRP hydrological data
message("Loading King County DNRP hydrological data")

# basins = st_read(paste0(in_data_geospatial, "/King County DNRP/KC_basins.shp")
# areas = st_read(paste0(in_data_geospatial, "/NHD_H_1711_HU4_Shape/Shape/NHDArea.shp")
# flowlines = st_read(paste0(in_data_geospatial, "/NHD_H_1711_HU4_Shape/Shape/NHDFlowline.shp")

# "The riparian delineations we have are for both fixed widths as well as “functionally dynamic” buffers that may be more ecologically relevant."
# "As Abood and others (2012) describe, the variable-width, or “dynamic,” buffer is likely more ecologically relevant because it accounts for factors that affect how a stream interacts and is influenced by the riparian zone. A recent literature review of riparian buffers done by King County (2019b) highlights the range of widths needed to maintain various functions (e.g., erosion control, shade). The review did not evaluate this concept of dynamic buffers, but the findings illustrate that riparian functions are maintained at different widths depending on a range of factors including but not limited to the size of stream, soil composition, vegetation type and age, etc.
sf_ripfb = st_read(paste0(in_data_geospatial, "/King County DNRP/RiparianBuffer_basin.shp"), quiet = TRUE)
sf_ripfb = sf_ripfb %>% filter(lengths(st_intersects(geometry, st_as_sf(study_area) %>% st_transform(st_crs(sf_ripfb)))) > 0)

# "The pattern, quality, and connectivity of riparian areas can also be really interesting and I’ve wondered how that may affect birds and other wildlife. For instance, in Seattle, the streams can have pretty decent but very narrow riparian areas because they are in canyons. In other suburban watersheds, the riparian areas can be vegetated but lack the tall evergreens or complex structure - they may be wider and even “greener” but not as functional? Maybe birds with small territories can fare OK in urban riparian areas?"

#########################################################################################
# Load NASA GEDI-Fusion forest structure data
#
# https://www.earthdata.nasa.gov/data/catalog/ornl-cloud-gedi-fusion-structure-2236-1
message("Loading NASA GEDI-Fusion forest structure data")

# "Foliage height diversity (a unitless index)"
fhd_raw = rast(paste0(in_data_geospatial, "/NASA/gedifusion_fhd_2020.tif"))
template = project(vect(study_area), crs(fhd_raw))
rast_gedi_fhd = mask(crop(fhd_raw, template), template)

plot(rast_gedi_fhd)

# "Fractional canopy cover (proportion)"
cover_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_cover_2020.tif"))
template   = project(vect(study_area), crs(cover_raw))
rast_gedi_cover = mask(crop(cover_raw, template), template)

plot(rast_gedi_cover)

# "Corresponds to the height at which 98% of the waveform energy is captured - comparable to a canopy height measure"
height_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_rh98_2020.tif"))
template   = project(vect(study_area), crs(height_raw))
rast_gedi_height = mask(crop(height_raw, template), template)

plot(rast_gedi_height)

# "Plant area vegetation density (PAVD); the proportion of vegetation within the 5-10 m stratum above ground surface (PAVD 5-10 m)""
pavd55o10m_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_pavd5to10m_2020.tif"))
template   = project(vect(study_area), crs(pavd55o10m_raw))
rast_gedi_pavd5to10m = mask(crop(pavd55o10m_raw, template), template)

plot(rast_gedi_pavd5to10m)

# "The proportion of vegetation density (PAVD) greater than 20 m above ground surface, chosen to represent the presence of a mature upper canopy within different forest types"

pavd20m_raw  = rast(paste0(in_data_geospatial, "/NASA/gedifusion_pavd20m_2020.tif"))
template   = project(vect(study_area), crs(pavd20m_raw))
rast_gedi_pavd20m = mask(crop(pavd20m_raw, template), template)

plot(rast_gedi_pavd20m)

#########################################################################################
# Project all data to same EPSG:32610 coordinate reference system
message("Projecting data to ", crs_standard)
study_area       = project(vect(study_area), crs_standard)

names(rast_nlcd_landcover)       = "rast_nlcd_landcover"
names(rast_nlcd_impervious)      = "rast_nlcd_impervious"
names(rast_usfs_canopycover)     = "rast_usfs_canopycover"
names(rast_landfire_vegcover)    = "rast_landfire_vegcover"
names(rast_landfire_vegheight)   = "rast_landfire_vegheight"
names(rast_landfire_treeheight)  = "rast_landfire_treeheight"
names(rast_landfire_shrubheight) = "rast_landfire_shrubheight"
names(rast_landfire_herbheight)  = "rast_landfire_herbheight"
names(rast_landfire_vegtype)     = "rast_landfire_vegtype"
names(rast_landfire_sclass)      = "rast_landfire_sclass"
names(rast_gedi_fhd)             = "rast_gedi_fhd"
names(rast_gedi_cover)           = "rast_gedi_cover"
names(rast_gedi_height)          = "rast_gedi_height"
names(rast_gedi_pavd5to10m)      = "rast_gedi_pavd5to10m"
names(rast_gedi_pavd20m)         = "rast_gedi_pavd20m"

rast_data = list(
  rast_nlcd_landcover       = project(rast_nlcd_landcover,       crs_standard),
  rast_nlcd_impervious      = project(rast_nlcd_impervious,      crs_standard),
  rast_usfs_canopycover     = project(rast_usfs_canopycover,     crs_standard),
  rast_landfire_vegcover    = project(rast_landfire_vegcover,    crs_standard),
  rast_landfire_vegheight   = project(rast_landfire_vegheight,   crs_standard),
  rast_landfire_treeheight  = project(rast_landfire_treeheight,  crs_standard),
  rast_landfire_shrubheight = project(rast_landfire_shrubheight, crs_standard),
  rast_landfire_herbheight  = project(rast_landfire_herbheight,  crs_standard),
  rast_landfire_vegtype     = project(rast_landfire_vegtype,     crs_standard),
  rast_landfire_sclass      = project(rast_landfire_sclass,      crs_standard),
  rast_gedi_fhd             = project(rast_gedi_fhd,             crs_standard),
  rast_gedi_cover           = project(rast_gedi_cover,           crs_standard),
  rast_gedi_height          = project(rast_gedi_height,          crs_standard),
  rast_gedi_pavd5to10m      = project(rast_gedi_pavd5to10m,      crs_standard),
  rast_gedi_pavd20m         = project(rast_gedi_pavd20m,         crs_standard)
)

sf_ripfb = sf_ripfb %>% st_transform(crs = crs_standard)

#########################################################################################
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
