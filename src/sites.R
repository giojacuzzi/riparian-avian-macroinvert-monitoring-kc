# Install dependencies
packages = c(
  'dplyr',
  'mapview',
  'sf',
  'tigris',
  'gdalUtils',
  'raster'
)
install.packages(packages[!(packages %in% rownames(installed.packages()))])

# Load dependencies
library(dplyr)
library(mapview)
library(sf)
library(tigris)
library(gdalUtils)
library(raster)

# Load, clean, and process data ############################################################

# Load Puget Sound Stream Benthos data
data_path = "data/Puget Sound Stream Benthos - ScoresByYear.txt"
data_raw = read.table(data_path, header = TRUE, sep = "\t", fill = TRUE)
# View(data) # Run this line to inspect the data directly

# Clean data
data = data_raw
data[data == ''] = NA # replace missing values with NA
is_numeric <- function(x) { !is.na(as.numeric(as.character(x))) }
data = data[is_numeric(data$Latitude) & is_numeric(data$Longitude), ] # remove invalid coords
data = data[!is.na(data$Longitude) & !is.na(data$Latitude), ] # remove points with no coords

# Subset sites sampled at least once within past several years
years = c('X2023', 'X2022', 'X2021', 'X2020')
data = data[, c('Site.ID', 'Basin', 'Stream', 'Agency', 'Project', 'Latitude', 'Longitude', years)] # subset data of interest
data[, years] = lapply(data[, years], as.numeric) # convert to numeric
data = data[rowSums(!is.na(data[, years])) >= 1, ]

# Calculate mean B-IBI over the chosen years
data$mean_BIBI <- rowMeans(data[, years], na.rm = TRUE)

# Categorize B-IBI means as low/mid/high terciles
cutoffs = quantile(data$mean_BIBI, probs = c(1/3, 2/3), na.rm = TRUE)
cutoffs
data$tercile_BIBI = cut(data$mean_BIBI, breaks = c(-Inf, cutoffs, Inf), labels = c("Low", "Mid", "High"))

# Calculate B-IBI trend over the chosen years as linear regression slope
year_numbers <- as.numeric(gsub('X', '', years))
calculate_slope = function(row) { # function to calculate the slope for each row
  non_na_indices = which(!is.na(row))
  if (length(non_na_indices) < 2) return(NA)
  non_na_values = row[non_na_indices]
  non_na_years  = year_numbers[non_na_indices]
  fit = lm(non_na_values ~ non_na_years)
  return(coef(fit)[2])
}
data$trend_BIBI = apply(data[, years], 1, calculate_slope)

# Convert to sf points #####################################################################

# Convert to shapefile points
data_sf = st_as_sf(data, coords=c('Longitude', 'Latitude'), crs='NAD83', agr='constant')
data_sf$long = st_coordinates(data_sf$geometry)[,'X']
data_sf$lat  = st_coordinates(data_sf$geometry)[,'Y']

# Load King County shapefile
wa_counties = counties(state = 'WA', cb = TRUE)
king_county = wa_counties[wa_counties$NAME == 'King', ]

# Subset sites within King County
sites = st_intersection(data_sf, king_county)

# Print number of sites per tercile
summary(sites$tercile_BIBI)

# Plot sites on a map, colored by tercile
mapview(sites, zcol='tercile_BIBI')

# Plot other variables as separate map layers
mapview(sites, zcol='mean_BIBI') + mapview(sites, zcol='trend_BIBI')

# Load impervious surface data ###################################################

# Create virtual raster VRTs pointing to IMGs without any modification
imp_raster_imgfile = '/Volumes/gioj_t7_1/NLCD/NLCD_impervious_2021_release_all_files_20230630/nlcd_2021_impervious_l48_20230630.img'
imp_raster_file    = '_output/nlcd_2021_impervious_l48_20230630.vrt'

# Create virtual raster VRT pointing to IMG without any modification
gdalbuildvrt(
  gdalfile = imp_raster_imgfile,
  output.vrt = imp_raster_file
)

# Albers equal-area projection
# aea = '+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs'

# Use gdalwarp to extract the county area from the NLCD impervious percentage raster (already in Albers projection)
polygon_file = '_output/king_county.gpkg'
if (file.exists(polygon_file)) file.remove(polygon_file)
raster_file = '_output/impervious_surface_raster.tif'
if (file.exists(raster_file)) file.remove(raster_file)

st_write(st_union(king_county), dsn = polygon_file, driver = 'GPKG', append = F)
gdalwarp(
  srcfile = imp_raster_file, dstfile = raster_file,
  cutline = polygon_file, crop_to_cutline = T,
  tr = c(30, 30), dstnodata = 'None'
)

impervious_surface = raster(raster_file)
king_county_trs = st_transform(king_county, st_crs(impervious_surface))
king_county_mask   = rasterize(king_county_trs, impervious_surface)
impervious_surface = mask(impervious_surface, king_county_mask)
mapview(impervious_surface)

# Calculate local 100 meter buffer impervious coverage % ########################################
sites = st_transform(sites, crs = st_crs(impervious_surface))
raster_with_zeros <- reclassify(impervious_surface, cbind(NA, 0))

# Create buffer around each point and get raster values inside
local_imp_coverage <- sapply(1:nrow(sites), function(i) {
  print(paste('Calculating local impervious coverage for site', i))
  site <- sites[i, ]
  buffer <- st_buffer(site, 100) # 100 meter recording range radius
  local_mask <- rasterize(buffer, raster_with_zeros)
  masked_raster <- mask(raster_with_zeros, local_mask)
  mean_value <- mean(masked_raster[], na.rm = TRUE) * 0.01
  return(mean_value)
})

sites$local_imp_coverage = local_imp_coverage

# Plot mean_BIBI against local_imp_coverage
plot(sites$local_imp_coverage, sites$mean_BIBI, main = 'local % impervious (100m buffer) x mean B-IBI')
regression <- lm(mean_BIBI ~ local_imp_coverage, data = sites)
abline(regression, col = "red")
summary(regression)$r.squared

hist(sites$local_imp_coverage)

mapview(sites, zcol='local_imp_coverage')

# Calculate regional drainage impervious coverage % ####################################

# Load washershed boundaries (https://ecology.wa.gov/water-shorelines/water-supply/water-availability/watershed-look-up)
WBDHU12 = st_read('/Volumes/gioj_t7_1/WBDHU/WBDHU12/WBDHU12.shp')
WBDHU12 = st_transform(WBDHU12, crs=st_crs(king_county))
WBDHU12 = st_intersection(WBDHU12, king_county)
WBDHU12 = st_transform(WBDHU12, crs=st_crs(impervious_surface))

# Calculate the ratio of impervious surface area to total area for each watershed
# Initialize a vector to store the ratio for each watershed
total_area <- numeric(nrow(WBDHU12))
total_imp_area <- numeric(nrow(WBDHU12))
impervious_ratio <- numeric(nrow(WBDHU12))

# Iterate over each watershed
for (i in 1:nrow(WBDHU12)) {
  print(paste('Calculating regional drainage impervious coverage for drainage', i))
  
  # Extract the raster values within the current watershed polygon
  watershed_mask <- rasterize(WBDHU12[i, ], impervious_surface)
  masked_raster <- mask(impervious_surface, watershed_mask)
  imp_sum = sum(masked_raster[], na.rm = TRUE) * 0.01

  # Calculate the total impervious surface area in the watershed
  total_impervious_area <- imp_sum * 900 # 30m x 30m = 900 square meters
  # print(paste('total_impervious_area', total_impervious_area))
  
  # Calculate the total area of the watershed in square meters
  watershed_area <- length(na.omit(masked_raster[])) * 900 #st_area(WBDHU12[i, ])
  # print(paste('watershed_area       ', watershed_area))
  
  # Calculate the ratio of impervious surface area to total area
  ratio = total_impervious_area / watershed_area
  # print(paste('ratio                ', ratio))
  impervious_ratio[i] <- ratio
  total_area[i] <- watershed_area
  total_imp_area[i] <- total_impervious_area
}

# Add the calculated ratios as a new column in the WBDHU12 data frame
WBDHU12$drainage_imp_coverage <- impervious_ratio

# Intersect sites with watershed impervious coverage
WBDHU12 = st_transform(WBDHU12, crs = st_crs(sites))
sites = st_intersection(sites, WBDHU12)

# Plot mean_BIBI against drainage_imp_coverage
plot(sites$drainage_imp_coverage, sites$mean_BIBI, main = 'regional % impervious (drainage) x mean B-IBI')
regression <- lm(mean_BIBI ~ drainage_imp_coverage, data = sites)
abline(regression, col = "red")
summary(regression)$r.squared

hist(sites$drainage_imp_coverage)

mapview(WBDHU12, zcol='drainage_imp_coverage') + mapview(sites, zcol='tercile_BIBI')

# Plot and save results to file #############################################

par(mfrow = c(2, 2))

plot(sites$local_imp_coverage, sites$mean_BIBI, main = 'local % impervious (100m buffer) x mean B-IBI')
regression <- lm(mean_BIBI ~ local_imp_coverage, data = sites)
abline(regression, col = "red")

plot(sites$drainage_imp_coverage, sites$mean_BIBI, main = 'regional % impervious (drainage) x mean B-IBI')
regression <- lm(mean_BIBI ~ drainage_imp_coverage, data = sites)
abline(regression, col = "red")

hist(sites$local_imp_coverage)

hist(sites$drainage_imp_coverage)

par(mfrow = c(1, 1))
plot(sites$drainage_imp_coverage, sites$local_imp_coverage, main = 'regional % impervious (drainage) x local % impervious (100m buffer)')

results = sites %>% st_drop_geometry()
write.csv(results, "potential_sites.csv", row.names = FALSE)
