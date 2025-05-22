# Load dependencies
library(dplyr)
library(mapview)
library(sf)
library(tigris)
library(raster)
library(terra)
library(ggplot2)

# Load, clean, and process data ############################################################

# Load Puget Sound Stream Benthos data
# Permalink: https://pugetsoundstreambenthos.org/Download.aspx?page=Download%2FScoresByYear.ashx&TR=-1&minY=2020
# Raw .txt opened in Excel and saved as .csv to resolve tab-delimited formatting issues
data_path = "data/raw/benthos/ScoresByYear.csv"
data_raw = read.csv(data_path)
# View(data) # Run this line to inspect the data directly

# Clean data
data = data_raw %>% janitor::clean_names()

# Subset sites sampled at least once within past several years
years = c('x2024', 'x2023', 'x2022', 'x2021', 'x2020')
data = data[, c('site_id', 'basin', 'stream', 'agency', 'project', 'latitude', 'longitude', years)] # subset data of interest
data[, years] = lapply(data[, years], as.numeric) # convert to numeric
data = data[rowSums(!is.na(data[, years])) >= 1, ]

# Calculate mean B-IBI over the chosen years
data$mean_BIBI = rowMeans(data[, years], na.rm = TRUE)

# Categorize B-IBI means as low/mid/high terciles
tercile_cutoffs = quantile(data$mean_BIBI, probs = c(1/3, 2/3), na.rm = TRUE)
data$mean_tercile_BIBI = cut(data$mean_BIBI, breaks = c(-Inf, tercile_cutoffs, Inf), labels = c("Low", "Mid", "High"))

category_cutoffs <- c(0, 20, 40, 60, 80, 100)
category_labels <- c("Very Poor", "Poor", "Fair", "Good", "Excellent")
data$mean_category_BIBI = cut(data$mean_BIBI,
                            breaks = category_cutoffs,
                            labels = category_labels,
                            right = FALSE,  # interval includes the left, excludes the right
                            include.lowest = TRUE)

# Calculate B-IBI trend over the chosen years as linear regression slope
year_numbers <- as.numeric(gsub('x', '', years))
calculate_slope = function(row) { # function to calculate the slope for each row
  non_na_indices = which(!is.na(row))
  if (length(non_na_indices) < 2) return(NA)
  non_na_values = row[non_na_indices]
  non_na_years  = year_numbers[non_na_indices]
  fit = lm(non_na_values ~ non_na_years)
  return(coef(fit)[2])
}
data$trend_slope_BIBI = apply(data[, years], 1, calculate_slope)

# Convert to sf points #####################################################################

# Convert to shapefile points
data_sf = st_as_sf(data, coords=c('longitude', 'latitude'), crs='NAD83', agr='constant')
data_sf$long = st_coordinates(data_sf$geometry)[,'X']
data_sf$lat  = st_coordinates(data_sf$geometry)[,'Y']

# Load King County shapefile
wa_counties = counties(state = 'WA', cb = TRUE)
king_county = wa_counties[wa_counties$NAME == 'King', ]

# Subset sites within King County
sites = st_intersection(data_sf, king_county)

# Print number of sites per tercile and category
summary(sites$mean_tercile_BIBI)
summary(sites$mean_category_BIBI)

# Plot sites on a map, colored by tercile
mapview(sites, zcol='mean_tercile_BIBI')
mapview(sites, zcol='mean_category_BIBI')

# Plot other variables as separate map layers
mapview(sites, zcol='mean_BIBI') + mapview(sites, zcol='trend_slope_BIBI')

# Load impervious surface data ###################################################

# Create virtual raster VRTs pointing to IMGs without any modification
# https://www.mrlc.gov/downloads/sciweb1/shared/mrlc/metadata/Annual_NLCD_FctImp_2023_CU_C1V0.xml
imp_raster_imgfile = 'data/raw/environment/Annual_NLCD_FctImp_2023_CU_C1V0.tif'

# Mask the raster to King County extent and consolidate coordinate systems
imp_raster  <- rast(imp_raster_imgfile)
king_county <- st_transform(king_county, crs(imp_raster))
imp_raster  <- crop(imp_raster, vect(king_county))
imp_raster  <- mask(imp_raster, vect(king_county))
sites <- st_transform(sites, crs(imp_raster))
# mapview(aggregate(imp_raster, fact = 4, fun = mean)) + mapview(sites)

# Calculate mean impervious surface coverage for site buffers ###################################################

# 100 meter buffer (~3.14 ha) roughly corresponds to the home range size of most riparian-associated species under consideration
buffer_sizes = c(100, 250, 500, 1000) # meters
sites_buffered = sites
for (buffer_size in buffer_sizes) {
  # Create circular buffers around site points
  buffers <- st_buffer(sites, dist = buffer_size)
  imp_buffers <- mask(imp_raster, vect(buffers))

  # Calculate mean impervious surface coverage within buffers
  imp_mean <- terra::extract(imp_raster, vect(buffers), fun = mean, na.rm = TRUE)
  col_name = paste0("imp_mean_", buffer_size, "m")
  sites_buffered[[col_name]] <- imp_mean[,2]
}
mapview(imp_buffers) + mapview(sites_buffered, zcol = col_name)

# Plot mean_BIBI against local_imp_coverage
for (buffer_size in buffer_sizes) {
  temp = sites_buffered
  temp = temp %>% mutate(imp_mean = temp[[paste0("imp_mean_", buffer_size,"m")]])
  p <- ggplot(sites_buffered %>% mutate(imp_mean = temp[[paste0("imp_mean_", buffer_size,"m")]]),
              aes(x = imp_mean, y = mean_BIBI)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "blue") +
    labs(title = paste0("Mean impervious surface coverage x mean B-IBI (", buffer_size," m)")) +
    xlim(0,100) + ylim(0,100)
  model <- lm(mean_BIBI ~ imp_mean, data = temp)
  rsq <- summary(model)$r.squared
  p = p + annotate("text", x = Inf, y = -Inf, label = paste0("R² = ", round(rsq, 3)), hjust = 1.1, vjust = -1.1, size = 5, color = "blue")
  print(p + theme_minimal())
}

# Note that B-IBI is most strongly affected by urbanization at the watershed, not local level

# Get previously surveyed sites
chosen_site_metadata = read.csv("data/site_metadata.csv") %>% janitor::clean_names()
aru_sites = st_as_sf(chosen_site_metadata, coords=c('aru_long', 'aru_lat'), crs='NAD83', agr='constant')
aru_sites$long = st_coordinates(aru_sites$geometry)[,'X']
aru_sites$lat  = st_coordinates(aru_sites$geometry)[,'Y']

sites_buffered = sites_buffered %>% mutate(previously_surveyed = site_id %in% aru_sites$benthos_site)
mapview(imp_buffers) + mapview(sites_buffered, zcol = "previously_surveyed") + mapview(aru_sites)

# Make failed surveys available again
sites_lost = c(171, 167, 152)
sites_buffered[sites_buffered$site_id %in% sites_lost, "previously_surveyed"] = FALSE

# Exclude sites with abs(trend_slope_BIBI) > 10
sites_available = sites_buffered %>% filter(abs(trend_slope_BIBI) < 2 | previously_surveyed == TRUE | site_id %in% sites_lost)

buffer_size = 100
p = ggplot(data = sites_available %>%
             mutate(imp_mean = sites_available[[paste0("imp_mean_", buffer_size,"m")]],
                    bibi_data_for_2024 = !is.na(sites_available$x2024)
                    ),
           aes(x = imp_mean, y = mean_BIBI)) +
  geom_point(aes(color = previously_surveyed, shape = bibi_data_for_2024), size = 3) +
  geom_text(aes(label = site_id), vjust = -1, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  labs(
    title = "Old and candidate sites",
    x = paste0("Mean impervious surface coverage (", buffer_size, " m)"),
    y = "Mean B-IBI",
    color = "Previously surveyed by ARU",
    shape = "B-IBI data collected in 2024"
  ) +
  theme_minimal(); print(p)

# View selected sites
site_list = c(171, 167, 152)

buffer_size = 100
sites_selected = sites_available %>% filter(site_id %in% site_list | previously_surveyed == TRUE)
p = ggplot(data = sites_selected %>%
             mutate(imp_mean = sites_selected[[paste0("imp_mean_", buffer_size,"m")]]),
           aes(x = imp_mean, y = mean_BIBI)) +
  geom_point(aes(color = !previously_surveyed), size = 3) +
  geom_text(aes(label = site_id), vjust = -1, size = 3) +
  scale_color_manual(values = c("FALSE" = "gray", "TRUE" = "blue")) +
  labs(
    title = "New and old sites",
    x = paste0("Mean impervious surface coverage (", buffer_size, " m)"),
    y = "Mean B-IBI",
    color = "New site selected to survey"
  ) +
  theme_minimal(); print(p)

mapview(sites_available %>% filter(site_id %in% site_list | previously_surveyed == TRUE), zcol = "previously_surveyed") 


# ## DISTRIBUTION OF LOCAL % IMPERVIOUS ############################################
# for (buffer_size in buffer_sizes) {
#   in_col_name = paste0("imp_mean_", buffer_size, "m")
#   out_col_name = paste0("imp_category_", buffer_size, "m")
#   breaks_imp_lmh = seq(min(sites_buffered[[in_col_name]]), max(sites_buffered[[in_col_name]]), length.out=4)
#   labels_imp_lmh = c("Low", "Mid", "High")
#   sites_buffered[[out_col_name]] = cut(sites_buffered[[in_col_name]], breaks = breaks_imp_lmh, labels = labels_imp_lmh, right = FALSE, include.lowest = TRUE)
#   summary(sites_buffered[[out_col_name]])
# }

# # COUNTS #####################################################################
# 
# # Number of potential sites per IMP/BIBI LMH combination
# buffer_size = 100
# sites_test %>% count(paste0("imp_category_", buffer_size, "m"), mean_tercile_BIBI)
# 
# # Print each row corresponding to the combinations
# counted_df <- sites %>%
#   group_by(imp_lmh, bibi_lmh) %>%
#   summarize(n = n(), .groups = "drop")
# for (i in 1:nrow(counted_df)) {
#   combination <- counted_df[i, c("imp_lmh", "bibi_lmh")]
#   subset_df <- sites[sites$imp_lmh == combination$imp_lmh & sites$bibi_lmh == combination$bibi_lmh, ]
#   print(subset_df)
# }
