source("src/global.R")

site_data_reach = readRDS("data/cache/3_calculate_vars/NEW_site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/NEW_site_data_5000m.rds")

in_path_nlcd_metadata   = "data/raw/nlcd_metadata.csv"
in_cache_geospatial_dir = "data/cache/2_preprocess_geospatial_data"

# Retrieve study area
study_area = st_read(paste0(in_cache_geospatial_dir, "/sf_studyarea.gpkg"), quiet = TRUE)

counties = counties(state = "WA", cb = TRUE) %>% st_transform(crs = crs_standard)
county_king = counties %>% filter(NAME == "King")

library(rnaturalearth)
states <- ne_states(country = "United States of America", returnclass = "sf") %>% filter(!(name_en %in% c("Alaska", "Hawaii")))
states$washington = FALSE
washington = states %>% filter(name == "Washington")

countries <- ne_countries(scale = "large", returnclass = "sf") %>% filter(name_en %in% c("Canada", "Mexico", "United States of America"))

bbox <- st_bbox(st_transform(site_data_reach, st_crs(washington)))

ggplot() +
  geom_sf(data = washington, fill = "lightgray", color = "black") +
  geom_sf(data = site_data_reach %>% st_transform(st_crs(washington)), aes(color = bibi), size = 2) +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), ylim = c(bbox["ymin"], bbox["ymax"]))

states_sf <- states(cb = TRUE, resolution = "5m", year = 2024) %>% filter(NAME == "Washington")

sf_flowline = st_read(paste0(in_cache_geospatial_dir, "/sf_flowline.gpkg"), quiet = TRUE) %>% clean_names() %>%
  filter(f_code %in% c(
    46000, 46003, 46006, 46007, # include only stream/river flowlines
    55800 # NOTE: Sammamish River is designated as arificial (55800) but is covered by USFS riparian area
  ))
sf_flowline$f_code = factor(sf_flowline$f_code)

# Transform flowlines to the CRS of your Washington polygon
sf_flowline_trans <- st_transform(sf_flowline, st_crs(washington))

# Optional: filter to just the study area to reduce plotting load
sf_flowline_sub <- st_intersection(sf_flowline_trans, st_transform(study_area, st_crs(sf_flowline_trans)))

message("- NHD waterbodies")
sf_waterbody = st_read(paste0(in_cache_geospatial_dir, "/sf_waterbody.gpkg"), quiet = TRUE) %>% clean_names()

sf_waterbody <- st_transform(sf_waterbody, st_crs(washington))

counties = counties %>% filter(NAME %in% c("King", "Pierce", "Snohomish"))

bbox_poly <- st_as_sfc(bbox)          # returns sfc_POLYGON
bbox_poly <- st_sf(geometry = bbox_poly)  # convert to sf for ggplot

large_waterbodies = sf_waterbody %>% filter(gnis_name %in% c("Puget Sound", "Lake Washington Ship Canal", "Lake Union", "Lake Washington", "Lake Sammamish"))


# Load cached basin sf objects and retain only those sampled
message("- HUC 12 basins")
sf_basins12d = st_read(paste0(in_cache_geospatial_dir, "/sf_basins12d.gpkg"), quiet = TRUE) %>%
  clean_names() %>% select(huc12, name, area_sq_km) %>% mutate(basin_name = name, basin_area = area_sq_km)
sf_basins12d = sf_basins12d %>% filter(lengths(st_intersects(., st_transform(sites_aru, st_crs(sf_basins12d)))) > 0)

# Load cached raster data
message("- Rasters")
rast_filepaths = list.files(in_cache_geospatial_dir, pattern = "^rast_.*\\.tif$", full.names = TRUE)
rast_data = lapply(rast_filepaths, rast)
names(rast_data) = gsub("\\.tif$", "", basename(rast_filepaths))

library(ggspatial)

pt_regions <- st_as_sf(data.frame(
  name = c("Washington"),
  lon  = c(-119.695328),
  lat  = c(47.497497)
), coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(washington))

p_region = ggplot() +
  geom_sf(data = countries, color = "grey10", fill = "gray80") +
  geom_sf(data = states, aes(fill = (iso_3166_2 == "US-WA")), color = "grey10") +
  geom_sf(data = bbox_poly, fill = "transparent", color = "red", linewidth = 0.75) +
  scale_fill_manual(values = c("gray80", "white")) +
  coord_sf(expand = FALSE, xlim = c(-126, -116), ylim = c(45, 50)) +
  geom_sf_text(data = pt_regions, aes(label = name), size = 1, fontface = "bold") +
  labs(x = "", y = "") +
  theme_sleek() +
  theme(legend.position = "none", axis.ticks = element_blank(), axis.text  = element_blank(),
        plot.background  = element_blank(),
        panel.border = element_rect(fill = NA, colour = "grey50", linewidth = 1)); print(p_region)

pt_cities <- st_as_sf(data.frame(
  name = c("Seattle", "Bellevue", "Bothell", "Renton", "Kent", "Enumclaw", "Snoqualmie", "Duvall", "Issaquah"),
  lon  = c(-122.316966, -122.210994, -122.197858, -122.186078, -122.244956, -121.998776, -121.874129, -121.969865, -122.062781),
  lat  = c(47.619934, 47.627918, 47.785999, 47.479029, 47.356098, 47.216641, 47.528785, 47.737381, 47.537253)
), coords = c("lon", "lat"), crs = 4326) %>% st_transform(st_crs(washington))

ggplot() +
  geom_sf(data = counties %>% st_transform(st_crs(washington)), fill = "white", color = "transparent") +
  geom_sf(data = sf_flowline_sub, color = "lightskyblue2", linewidth = 0.2) +
  geom_sf(data = sf_waterbody, color = "transparent", fill = "lightskyblue2", linewidth = 0) +
  geom_sf(data = sf_basins12d, color = "grey60", fill = "transparent", linewidth = 0.35) +
  geom_sf(data = counties %>% st_union() %>% st_transform(st_crs(washington)), fill = "transparent", color = "grey30") +
  geom_sf(data = large_waterbodies, color = "grey30", fill = "lightskyblue2") +
  geom_sf(data = site_data_reach %>% st_transform(st_crs(washington)), aes(fill = bibi), color = "grey10", shape = 21, size = 3) +
  coord_sf(expand = FALSE, xlim = c(-122.44068, -121.56035), ylim = c(47.19100, 47.79312)) +
  scale_fill_viridis_c(option = "viridis") +
  geom_sf_text(data = pt_cities, aes(label = name), size = 1.5) +
  labs(x = "", y = "", fill = "B-IBI") +
  annotation_scale(location = "bl", width_hint = 0.25) +
  annotation_north_arrow(location = "br", which_north = "true", style = north_arrow_minimal,
    height = unit(0.75, "cm"), width = unit(0.75, "cm")
  ) +
  theme_sleek() +
  inset_element(
    p_region,
    left = 0.65, right = 1.0,
    bottom = 0.75, top = 1.0,
    align_to = "panel"
  )






# site_id_focal = "262"
# site_focal = site_data_reach %>% filter(site_id == site_id_focal)
# site_focal_buffer = st_buffer(site_focal, 550)
# site_focal = site_focal %>% st_transform(st_crs(washington))
# 
# r_utm <- project(rast_data$rast_riparian, crs_standard)
# df_riparian <- as.data.frame(r_utm, xy = TRUE)
# df_riparian$rast_riparian = round(df_riparian$rast_riparian)
# 
# r = project(rast_data$rast_riparian, crs(sf_flowline_sub))
# df_r = as.data.frame(r, xy = TRUE)
# 
# pnt_aru = sites_aru %>% filter(site_id == site_id_focal)
# buffer_aru = st_buffer(pnt_aru %>% st_transform(crs_standard), 550)
# pnt_pssb = sites_pssb %>% filter(site_id == site_id_focal)
# bbox_buffer = st_bbox(buffer_aru %>% st_transform(st_crs(sf_flowline_sub)))
# 
# r <- rast_data$rast_riparian
# buffer_rast_crs <- st_transform(buffer_aru, crs(r))
# buffer_vect <- vect(buffer_rast_crs)
# r_crop <- crop(r, buffer_vect)
# r_mask <- mask(r_crop, buffer_vect)
# r_df <- as.data.frame(r_mask %>% project(crs(sf_flowline_sub)), xy = TRUE)
# 
# sf_flowline_crop = st_intersection(sf_flowline_sub %>% st_transform(st_crs(buffer_aru)), buffer_aru)
# 
# sf_waterbody_crop = st_intersection(sf_waterbody %>% st_transform(st_crs(buffer_aru)), buffer_aru)
# 
# emergence_zone = st_make_valid(st_union(st_buffer(sf_flowline_crop, 100)))
# emergence_zone = st_intersection(emergence_zone, buffer_aru %>% st_transform(st_crs(emergence_zone)))
# 
# ggplot() +
#   geom_sf(data = buffer_aru %>% st_transform(st_crs(sf_flowline_sub)), color = "black", fill = "white", linewidth = 0.5) +
#   geom_raster(data = r_df, aes(x = x, y = y), fill = "forestgreen") +
#   geom_sf(data = emergence_zone %>% st_transform(st_crs(sf_flowline_sub)), color = "transparent", fill = "lightskyblue2", alpha = 0.5, linewidth = 0.5) +
#   geom_sf(data = sf_flowline_crop %>% st_transform(st_crs(sf_flowline_sub)), color = "lightskyblue2", linewidth = 1) +
#   geom_sf(data = sf_waterbody_crop %>% st_transform(st_crs(sf_flowline_sub)), color = "lightskyblue2", fill = "lightskyblue2", linewidth = 0) +
#   geom_sf(data = pnt_aru, size = 3, color = "red") +
#   geom_sf(data = pnt_pssb, size = 3, color = "blue") +
#   coord_sf(xlim = c(bbox_buffer["xmin"], bbox_buffer["xmax"]),
#            ylim = c(bbox_buffer["ymin"], bbox_buffer["ymax"]),
#            datum = st_crs(buffer_aru),
#            expand = FALSE)
# 
# ggplot() +
#   geom_raster(data = df_riparian, aes(x = x, y = y), fill = "forestgreen") +
#   coord_sf(datum = st_crs(washington), expand = FALSE)
# 
# ggplot() +
#   geom_sf(data = test %>% st_transform(st_crs(washington)), fill = "white", color = "transparent") +
#   # geom_sf(data = sf_flowline_sub, color = "lightskyblue2", linewidth = 0.2) +
#   # geom_sf(data = sf_waterbody, color = "lightskyblue2", fill = "lightskyblue2") +
#   geom_sf(data = sf_basins12d, color = "darkgray", fill = "transparent", linewidth = 0.35) +
#   geom_sf(data = test %>% st_transform(st_crs(washington)), fill = "transparent", color = "black") +
#   geom_sf(data = large_waterbodies, color = "black", fill = "lightskyblue2") +
#   geom_sf(data = site_data_reach %>% st_transform(st_crs(washington)), 
#           aes(color = bibi), size = 3) +
#   coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]),
#            ylim = c(bbox["ymin"], bbox["ymax"]))
