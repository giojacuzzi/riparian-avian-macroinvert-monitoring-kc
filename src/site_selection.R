# Path to .csv output of potential_sites.R
data_path = "_output/potential_sites.csv"

# Install dependencies
packages = c(
  'dplyr',
  'mapview',
  'sf'
)
install.packages(packages[!(packages %in% rownames(installed.packages()))])

# Load dependencies
library(sf)
library(mapview)
library(dplyr)

# Read potential sites from file
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

# Filter for King County DNRP sites only
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
breaks_bibi_category = c(0,20,40,60,80,100)
labels_bibi_category = c("Very poor", "Poor", "Fair", "Good", "Excellent")
hist(sites$mean_BIBI, breaks = breaks_bibi_category)
sites$bibi_category = cut(sites$mean_BIBI, breaks = breaks_bibi_category, labels = labels_bibi_category, right = FALSE, include.lowest = TRUE)
summary(sites$bibi_category)

# LOW-MID-HIGH RANKING
# [0, 33)  - Low (Very Poor-Poor)
# [33,66)  - Mid (Poor-Good)
# [66,100] - High (Good-Excellent)
breaks_bibi_lmh = c(0,33,66,100)
labels_bibi_lmh = c("Low", "Mid", "High")
hist(sites$mean_BIBI, breaks = breaks_bibi_lmh, freq = TRUE)
sites$bibi_lmh = cut(sites$mean_BIBI, breaks = breaks_bibi_lmh, labels = labels_bibi_lmh, right = FALSE, include.lowest = TRUE)
summary(sites$bibi_lmh)

## DISTRIBUTION OF LOCAL % IMPERVIOUS ############################################

hist(sites$local100_imp_coverage)

breaks_imp_lmh = c(0.0,0.2,0.4,1.0)
labels_imp_lmh = c("Low", "Mid", "High")
sites$imp_lmh = cut(sites$local100_imp_coverage, breaks = breaks_imp_lmh, labels = labels_imp_lmh, right = FALSE, include.lowest = TRUE)
summary(sites$imp_lmh)
sites$impDrainage_lmh = cut(sites$drainage_imp_coverage, breaks = breaks_imp_lmh, labels = labels_imp_lmh, right = FALSE, include.lowest = TRUE)
summary(sites$impDrainage_lmh)

plot(sites$local100_imp_coverage, sites$mean_BIBI, main = 'local % impervious (100m buffer) x mean B-IBI')
regression <- lm(mean_BIBI ~ local100_imp_coverage, data = sites)
abline(regression, col = "red")

# COUNTS #####################################################################

# Number of potential sites per IMP/BIBI LMH combination
sites %>% count(imp_lmh, bibi_lmh)

# Print each row corresponding to the combinations
counted_df <- sites %>%
  group_by(imp_lmh, bibi_lmh) %>%
  summarize(n = n(), .groups = "drop")
for (i in 1:nrow(counted_df)) {
  combination <- counted_df[i, c("imp_lmh", "bibi_lmh")]
  subset_df <- sites[sites$imp_lmh == combination$imp_lmh & sites$bibi_lmh == combination$bibi_lmh, ]
  print(subset_df)
}

# MAPPING ####################################################################

# Map all potential sites
sites_sf = st_as_sf(sites, coords=c('long', 'lat'), crs='NAD83', agr='constant')
sites_sf$long = st_coordinates(sites_sf$geometry)[,'X']
sites_sf$lat  = st_coordinates(sites_sf$geometry)[,'Y']
mapview(sites_sf, zcol='imp_lmh') + mapview(sites_sf, zcol='bibi_lmh')

# mapview(sites_sf[sites_sf$imp_lmh == 'Low' & sites_sf$bibi_lmh == 'Low', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'Low' & sites_sf$bibi_lmh == 'Mid', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'Low' & sites_sf$bibi_lmh == 'High', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'Mid' & sites_sf$bibi_lmh == 'Low', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'Mid' & sites_sf$bibi_lmh == 'Mid', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'Mid' & sites_sf$bibi_lmh == 'High', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'High' & sites_sf$bibi_lmh == 'Low', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'High' & sites_sf$bibi_lmh == 'Mid', ])
# mapview(sites_sf[sites_sf$imp_lmh == 'High' & sites_sf$bibi_lmh == 'High', ])

# Selected sites
# NOTE: Aim for 3 sites per IMP/BIBI combination, 10 per LMH gradient each.
# Avoid sites with erratic BIBI trends, and try to get ~6 sites per bibi_category as well.
selected_site_ids = c(
  122, 210, 211, # LL (pass: 190, 1628)
  121, 136, 193, # LM (pass: multiple)
  125, 155, 158, 160, 2217, # LH (pass: multiple)
  209, 213, 189, # ML (pass: 162, 208, 214, 1635)
  159, 167, 251, 252, 264, # MM (pass: 2366)
  150, 153, 262,      # MH
  161, 188, 191, 216, # HL
  171, 174, 270,      # HM
  152                 # HH
)
selected_sites = sites %>% filter(Site.ID %in% selected_site_ids)

# Selected sites count per gradient category
selected_sites %>% count(imp_lmh, bibi_lmh)
selected_sites %>% count(imp_lmh)
selected_sites %>% count(bibi_lmh)
selected_sites %>% count(bibi_category)

# MAP SELECTED SITES
selected_sites_sf = st_as_sf(selected_sites, coords=c('long', 'lat'), crs='NAD83', agr='constant')
selected_sites_sf$long = st_coordinates(selected_sites_sf$geometry)[,'X']
selected_sites_sf$lat  = st_coordinates(selected_sites_sf$geometry)[,'Y']
mapview(selected_sites_sf, zcol='imp_lmh') + mapview(selected_sites_sf, zcol='bibi_lmh')

# Save selected sites to file with specific columns
final_sites = selected_sites[, c('Site.ID', 'lat', 'long', 'Stream', 'Basin', 'mean_BIBI', 'local100_imp_coverage', 'bibi_category', 'bibi_lmh', 'imp_lmh')] %>% arrange(Basin, Site.ID)
write.csv(final_sites, '/Users/giojacuzzi/Downloads/final_sites.csv', row.names = FALSE)
