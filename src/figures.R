source("src/global.R")

site_data_reach = readRDS("data/cache/3_calculate_vars/site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/site_data_5000m.rds")

in_path_nlcd_metadata   = "data/raw/nlcd_metadata.csv"
in_cache_geospatial_dir = "data/cache/2_preprocess_geospatial_data"

out_dir = "data/cache/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

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

bbox_poly = st_as_sfc(bbox)
bbox_poly = st_sf(geometry = bbox_poly)

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

fig_2A = ggplot() +
  geom_sf(data = counties %>% st_transform(st_crs(washington)), fill = "white", color = "transparent") +
  geom_sf(data = sf_flowline_sub, color = "lightskyblue2", linewidth = 0.2) +
  geom_sf(data = sf_waterbody, color = "transparent", fill = "lightskyblue2", linewidth = 0) +
  geom_sf(data = sf_basins12d, color = "grey70", fill = "transparent", linewidth = 0.35) +
  geom_sf(data = counties %>% st_union() %>% st_transform(st_crs(washington)), fill = "transparent", color = "grey60") +
  geom_sf(data = large_waterbodies, color = "grey60", fill = "lightskyblue2") +
  geom_sf(data = site_data_basin %>% filter(!site_id %in% c(sites_to_exclude)) %>%
            st_transform(st_crs(washington)), aes(size = bibi, fill = rast_nlcd_impervious_mean), color = "grey10", shape = 21, alpha = 0.75) +
  scale_fill_gradientn(
    colors = c("forestgreen", "lightcoral", "red", "darkred"),
    values = scales::rescale(c(0, 35, 65, 90)),
    na.value = "white"
  ) +
  coord_sf(expand = FALSE, xlim = c(-122.45, -121.6), ylim = c(47.25, 47.78)) +
  geom_sf_text(data = pt_cities, aes(label = name), size = 1.5) +
  labs(x = "", y = "", fill = "ISC", size = "B-IBI") +
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
  ) + theme(legend.position = "none"); print(fig_2A)

ggsave(paste0(out_dir, "/fig_2A.pdf"), fig_2A, width = 6, height = 6)

# Inset
world = ne_countries(scale = "medium", returnclass = "sf")
usa = rnaturalearth::ne_states(country = "United States of America", returnclass = "sf")
wa = usa[usa$name == "Washington", ]
outline = st_sfc(st_point(c(0,0)), crs = "+proj=ortho +lat_0=47 +lon_0=-120") |> st_buffer(dist = 6378137)
fig_2A_inset = ggplot() +
  geom_sf(data = world, fill = "gray90", color = "gray40", linewidth = 0.2) +
  geom_sf(data = wa, fill = "red", color = "gray20") +
  geom_sf(data = outline, fill = NA, color = "gray20", linewidth = 0.2) +
  coord_sf(crs = "+proj=ortho +lat_0=47 +lon_0=-120"); print(fig_2A_inset)

ggsave(paste0(out_dir, "/fig_2A_inset.pdf"), fig_2A_inset, width = 2, height = 2)

