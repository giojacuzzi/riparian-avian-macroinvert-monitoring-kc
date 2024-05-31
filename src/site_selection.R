data_path = "_output/potential_sites.csv"

library(sf)
library(mapview)
library(dplyr)

sites = read.csv(data_path)

columns = c(
  "Site.ID",
  "long",
  "lat",
  "Name",
  "Basin",
  "Stream",
  "Agency",
  "Project",
  "X2023",
  "X2022",
  "X2021",
  "X2020",
  "mean_BIBI",
  "tercile_BIBI",
  "trend_BIBI",
  "local100_imp_coverage",
  "local200_imp_coverage",
  "drainage_imp_coverage"
)

sites = sites[, columns]

# Filter for DNRP sites
sites = sites[sites$Agency == 'King County - DNRP',]

## DISTRIBUTION OF MEAN B-IBI ############################################

# https://pugetsoundstreambenthos.org/About-BIBI.aspx
#
# CATEGORY RANKING - Modified from Karr et al. (1986) by Morley (2000).
# [80,100] - Excellent
# [60,80)  - Good
# [40,60)  - Fair
# [20,40)  - Poor
# [0,20)   - Very poor
rank_breaks_bibi_category = c(0,20,40,60,80,100)
rank_labels_bibi_category = c("Very poor", "Poor", "Fair", "Good", "Excellent")
hist(sites$mean_BIBI, breaks = rank_breaks_bibi_category)
sites$rank_bibi_category = cut(sites$mean_BIBI, breaks = rank_breaks_bibi_category, labels = rank_labels_bibi_category, right = FALSE, include.lowest = TRUE)
summary(sites$rank_bibi_category)

# LOW-MID-HIGH RANKING
# [0, 33)  - Low (Very Poor-Poor)
# [33,66)  - Mid (Poor-Good)
# [66,100] - High (Good-Excellent)
rank_breaks_bibi_lmh = c(0,33,66,100)
rank_labels_bibi_lmh = c("Low", "Mid", "High")
hist(sites$mean_BIBI, breaks = rank_breaks_bibi_lmh, freq = TRUE)
sites$rank_bibi_lmh = cut(sites$mean_BIBI, breaks = rank_breaks_bibi_lmh, labels = rank_labels_bibi_lmh, right = FALSE, include.lowest = TRUE)
summary(sites$rank_bibi_lmh)

# TODO: Select sites such that there are 10 per LOW-MID-HIGH ranking and 6 per CATEGORY ranking

## DISTRIBUTION OF LOCAL % IMPERVIOUS ############################################

hist(sites$local100_imp_coverage)

rank_breaks_imp_lmh = c(0.0,0.2,0.4,1.0)
rank_labels_imp_lmh = c("Low", "Mid", "High")
sites$rank_imp_lmh = cut(sites$local100_imp_coverage, breaks = rank_breaks_imp_lmh, labels = rank_labels_imp_lmh, right = FALSE, include.lowest = TRUE)
summary(sites$rank_imp_lmh)

plot(sites$local100_imp_coverage, sites$mean_BIBI, main = 'local % impervious (100m buffer) x mean B-IBI')
regression <- lm(mean_BIBI ~ local100_imp_coverage, data = sites)
abline(regression, col = "red")

# COUNTS #####################################################################
sites %>% count(rank_imp_lmh, rank_bibi_lmh)

counted_df <- sites %>%
  group_by(rank_imp_lmh, rank_bibi_lmh) %>%
  summarize(n = n(), .groups = "drop")

# View the result
print(counted_df)

# Print each row corresponding to the combinations
for (i in 1:nrow(counted_df)) {
  combination <- counted_df[i, c("rank_imp_lmh", "rank_bibi_lmh")]
  subset_df <- sites[sites$rank_imp_lmh == combination$rank_imp_lmh & sites$rank_bibi_lmh == combination$rank_bibi_lmh, ]
  print(subset_df)
}

# MAPPING ####################################################################

sites_sf = st_as_sf(sites, coords=c('long', 'lat'), crs='NAD83', agr='constant')
sites_sf$long = st_coordinates(sites_sf$geometry)[,'X']
sites_sf$lat  = st_coordinates(sites_sf$geometry)[,'Y']

mapview(sites_sf, zcol='rank_imp_lmh') + mapview(sites_sf, zcol='rank_bibi_lmh')
