# 1_preprocess_agg_pam_data.R ==========================================================
# Derive and cache putative species detection histories and diversity metrics from aggregated PAM data
#
# Inputs:
# Aggregated predictions from classifier and metadata
path_prediction_data = "data/cache/0_aggregate_raw_pam_data/prediction_data.feather"
path_species_list = "data/pam/species_list.txt"
path_avonet_traits = "data/traits/AVONET Supplementary dataset 1.xlsx"
# Naive thresholds
threshold_classifier_score = 0.95 # Naive classifier minimum confidence score threshold to assume binary presence/absence
threshold_detected_days = 5 # Minimum number of unique days detected to retain species detections at a site
# Outputs:
out_detect_hist_data = "data/cache/1_preprocess_agg_pam_data/detections.rds"

# Load required packages (automatically install any missing)
pkgs = c(
  "tidyverse",  # data manipulation
  "arrow"       # cache data compression
)
sapply(pkgs, function(pkg) {
  if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
  library(pkg, character.only = TRUE)
  as.character(packageVersion(pkg))
})

# Load classifier prediction data
message("Loading classifier prediction data")
prediction_data = read_feather(path_prediction_data) %>% rename(site_id = site)

# Clean data --------------------------------------------------------------------------------

# Obtain putative detections with naive threshold
message("Obtaining putative detections with naive threshold")
detections = prediction_data %>% filter(confidence >= threshold_classifier_score)

# NOTE: Shared survey start/end dates were previously determined in 0_agg_raw_pam_data.R
start_date_2024 = min(detections$date[detections$season == "2024"])
end_date_2024 = max(detections$date[detections$season == "2024"])
start_date_2025 = min(detections$date[detections$season == "2025"])
end_date_2025 = max(detections$date[detections$season == "2025"])

# Only retain detections for species detected a minimum number of days at a site
message("Discarding  detections for species that were detected during less than ", threshold_detected_days, " surveys (days) at a site")
detections = detections %>%
  group_by(site_id, common_name) %>%
  summarise(n_surveys = n_distinct(date), .groups = "drop") %>%
  filter(n_surveys >= threshold_detected_days) %>%
  inner_join(detections, by = c("site_id", "common_name"))

# Determine detected species and detections per site ----------------------------------------

species_names = read_lines(path_species_list) %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>%
  filter(common_name %in% sort(unique(detections$common_name)))

message(nrow(species_names), " species detected:")
species_summary = detections %>% group_by(common_name) %>%
  summarise(
    n_sites = n_distinct(site_id),
    n_detections = n()
  ) %>% arrange(desc(n_sites))
print(species_summary, n = nrow(species_summary))

message("Detections per site:")
site_summary = detections %>% group_by(site_id) %>%
  summarise(
    n_species = n_distinct(common_name),
    n_detections = n()
  ) %>% arrange(desc(n_species))
print(site_summary, n = nrow(site_summary))

## Format detection data into species x site x survey ----------------------------------------

# Determine survey numbers by year
clean_detections <- detections %>%
  mutate(
    survey_num = case_when(
      year(date) == 2024 ~ as.integer(difftime(date, start_date_2024, units = "days")) + 1,
      year(date) == 2025 ~ as.integer(difftime(date, start_date_2025, units = "days")) + 1,
      TRUE ~ NA
    )
  ) %>% arrange(site_id, survey_num)

# Make a complete grid of site × species x surveys
template_all_site_species_surveys <- expand.grid(
  site_id = unique(clean_detections$site_id),
  common_name = unique(clean_detections$common_name),
  survey_num = 1:max(clean_detections$survey_num), # include all survey numbers
  stringsAsFactors = FALSE
)

# Summarize detections (long)
species_site_survey_long = clean_detections %>%
  group_by(site_id, common_name, survey_num) %>%
  summarise(n_detections = n(), .groups = "drop") %>%
  right_join(template_all_site_species_surveys, by = c("site_id", "common_name", "survey_num")) %>%
  replace_na(list(n_detections = 0)) %>% arrange(site_id, common_name, survey_num)

# Visualize detections per species across all sites as a function of survey number
ggplot(species_site_survey_long %>% group_by(common_name, survey_num) %>%
       summarise(total_detections = sum(n_detections), .groups = "drop"),
       aes(x = survey_num, y = total_detections)) +
  geom_line() + facet_wrap(~ common_name, scales = "free_y") +
  labs(
    title = "Species detections across all sites and surveys",
    x = "Survey number",
    y = "Total detections across sites"
  )

# Summarize detections (wide)
species_site_survey_wide = species_site_survey_long %>%
  pivot_wider(
    names_from = survey_num,
    values_from = n_detections,
    values_fill = 0
  )
survey_cols = setdiff(names(species_site_survey_wide), c("site_id", "common_name"))
species_site_survey_wide = species_site_survey_wide %>%
  relocate(all_of(survey_cols[order(as.numeric(survey_cols))]), .after = common_name)

# Visualize assemblage composition across sites ----------------------------------------

# Load species trait metadata
avonet = readxl::read_xlsx(path_avonet_traits, sheet = "AVONET2_eBird") %>%
janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  filter(scientific_name %in% species_names$scientific_name) %>%
  select(scientific_name, family, order, mass, habitat, habitat_density, migration, trophic_level, trophic_niche, primary_lifestyle)
species_metadata = left_join(species_names, avonet, by = "scientific_name")

# Summarize number of species per site and trophic niche
trophic_niche_per_site = left_join(species_site_survey_long, species_metadata, by = "common_name") %>% filter(n_detections > 0) %>%
  group_by(site_id, common_name, trophic_niche) %>% summarise(total_detections = sum(n_detections), .groups = "drop") %>%
  group_by(site_id, trophic_niche) %>% summarise(species_count = n_distinct(common_name), .groups = "drop")
print(ggplot(trophic_niche_per_site %>% left_join(site_summary, by = "site_id"),
       aes(x = species_count, y = reorder(site_id, n_species), fill = trophic_niche)) +
  geom_bar(stat = "identity") +
  theme_minimal())

# Summarize number of species per site and primary lifestyle
primary_lifestyle_per_site = left_join(species_site_survey_long, species_metadata, by = "common_name") %>% filter(n_detections > 0) %>%
  group_by(site_id, common_name, primary_lifestyle) %>% summarise(total_detections = sum(n_detections), .groups = "drop") %>%
  group_by(site_id, primary_lifestyle) %>% summarise(species_count = n_distinct(common_name), .groups = "drop")
print(ggplot(primary_lifestyle_per_site %>% left_join(site_summary, by = "site_id"),
       aes(x = species_count, y = reorder(site_id, n_species), fill = primary_lifestyle)) +
  geom_bar(stat = "identity") +
  theme_minimal())

# Cache putative species detection history and diversity data ------------------------------

detect_hist_data = list(
  threshold_classifier_score = threshold_classifier_score,
  threshold_detected_days = threshold_detected_days,
  long = species_site_survey_long,
  wide = species_site_survey_wide
)

if (!dir.exists(dirname(out_detect_hist_data))) dir.create(dirname(out_detect_hist_data), recursive = TRUE)
saveRDS(detect_hist_data, out_detect_hist_data)
message(crayon::green("Cached ", out_detect_hist_data))
