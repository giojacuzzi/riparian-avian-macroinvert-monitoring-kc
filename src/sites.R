# Install dependencies
packages = c(
  'dplyr',
  'mapview',
  'sf',
  'tigris'
)
install.packages(packages[!(packages %in% rownames(installed.packages()))])

# Load dependencies
library(dplyr)
library(mapview)
library(sf)
library(tigris)

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
