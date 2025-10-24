# 0_agg_raw_pam_data.R ===============================================================
# Aggregate raw prediction (.csv) data and associated metadata
#
# Directory structure is ".../season/deployment/serialno/date/serialno_date_time.csv"
# For example: "predictions/2020/Deployment1/S4A04271_20200412_Data/S4A04271_20200411/S4A04271_20200411_235938.BirdNET.results.csv"
#
# Defines a "season" (primary period) by year, and a "survey" (secondary period) as a
# continuous 24h recording period (from 00:00:00 to 23:59:59).Assumes that the entire
# active survey period is represented with individual .csv files (i.e. a recording with
# no detections should still be represented by an empty .csv file)
#
# Inputs:
# Root directory containing raw prediction data with metadata in structure defined above
root_dir_in = "data/pam/predictions"
pred_filetype = ".BirdNET.results.csv"
#
# Outputs:
# Cached dataframe of prediction file counts (i.e. recordings) per unit-survey
# Cached dataframe of all predictions
path_out_dir = "data/cache/0_aggregate_raw_pam_data"

path_out_survey_file_counts = paste0(path_out_dir, '/survey_file_counts.feather')
path_out_prediction_data = paste0(path_out_dir, '/prediction_data.feather')

library(arrow)
library(tidyverse)
library(janitor)
library(progress)

# Parse metadata from survey .csv path
parse_metadata = function(file_path) {
  
  # Get path relative to root
  rel_path = sub(paste0("^", root_dir_in, "/?"), "", file_path)
  
  # Split the file path into metadata components
  metadata_components = unlist(strsplit(rel_path, .Platform$file.sep))
  
  # Extract metadata
  season = metadata_components[1]
  site   = metadata_components[2]
  date_time = metadata_components[3]
  date_time = strsplit(date_time, pred_filetype)[[1]]
  date = strsplit(date_time, '_')[[1]][1]
  time = strsplit(date_time, '_')[[1]][2]
  datetime_start = as.POSIXct(paste(date,time), format = "%Y%m%d %H%M%S", tz = "UTC")
  yday_start = yday(datetime_start)
  
  return(list(
    season = season,
    site = site,
    datetime_start = datetime_start,
    yday_start = yday_start
  ))
}


message("Locating all ", pred_filetype, " prediction files")
survey_files = list.files(path = root_dir_in, pattern = paste0("\\", pred_filetype, "$"), recursive = TRUE, full.names = TRUE)
message("Located ", length(survey_files), " prediction files")

message("Parsing metadata for each prediction file")
metadata_list = lapply(seq_along(survey_files), function(f) {
  c(parse_metadata(survey_files[[f]]), file_path = survey_files[[f]])
})
predictions_metadata = data.table::rbindlist(metadata_list, fill = TRUE)

n_sites = length(unique(predictions_metadata$site))
n_seasons = length(unique(predictions_metadata$season))

message("Located predictions for ", n_sites, " sites across ", n_seasons, " seasons")

message("Visualizing survey effort")
survey_summary = predictions_metadata %>%
  group_by(site) %>%
  summarise(
    n_unique_dates = n_distinct(as.Date(datetime_start)),  # number of unique calendar dates
    min_yday = min(yday_start),
    min_date = min(as.Date(datetime_start)),
    max_yday = max(yday_start),
    max_date = max(as.Date(datetime_start)),
    .groups = "drop"
  )
ggplot(survey_summary) +
  geom_segment(aes(x = min_date, xend = max_date, 
                   y = reorder(site, min_date), yend = site), linewidth = 2) +
  labs(x = "Date", y = "Site", title = "Survey effort") +
  theme_minimal() +
  theme(legend.position = "none")

# Discard site 3097
sites_to_discard = c("3097")
message("Discarding ", length(sites_to_discard), " site(s) with insufficient survey effort:")
message(sites_to_discard)
predictions_metadata = predictions_metadata %>% filter(!(site %in% sites_to_discard))
survey_summary = predictions_metadata %>%
  group_by(site) %>%
  summarise(
    n_unique_dates = n_distinct(as.Date(datetime_start)),  # number of unique calendar dates
    min_yday = min(yday_start),
    min_date = min(as.Date(datetime_start)),
    max_yday = max(yday_start),
    max_date = max(as.Date(datetime_start)),
    .groups = "drop"
  )

# Filter for 2024 only
metadata_2024 <- predictions_metadata %>%
  filter(season == "2024") %>%
  mutate(date = as.Date(datetime_start))

# Count how many sites have data per date
common_dates <- metadata_2024 %>%
  group_by(date) %>%
  summarise(n_sites_present = n_distinct(site), .groups = "drop")

# Earliest date where all sites surveyed
earliest_common_date_2024 <- common_dates %>%
  filter(n_sites_present == length(unique(metadata_2024$site))) %>%
  summarise(min_date = min(date)) %>% pull(min_date)

latest_common_date_2024 <- common_dates %>%
  filter(n_sites_present == length(unique(metadata_2024$site))) %>%
  summarise(max_date = max(date)) %>% pull(max_date)

# Filter for 2025 only
metadata_2025 <- predictions_metadata %>%
  filter(season == "2025") %>%
  mutate(date = as.Date(datetime_start))

# Count how many sites have data per date
common_dates <- metadata_2025 %>%
  group_by(date) %>%
  summarise(n_sites_present = n_distinct(site), .groups = "drop")

# Earliest date where all sites surveyed
earliest_common_date_2025 <- common_dates %>%
  filter(n_sites_present == length(unique(metadata_2025$site))) %>%
  summarise(min_date = min(date)) %>% pull(min_date)

latest_common_date_2025 <- common_dates %>%
  filter(n_sites_present == length(unique(metadata_2025$site))) %>%
  summarise(max_date = max(date)) %>% pull(max_date)

start_date_2024 = earliest_common_date_2024 + 1
end_date_2024   = latest_common_date_2024 - 2

message("2024: start date ", start_date_2024, ", end date ", end_date_2024)
message("Days: ", end_date_2024 - start_date_2024)
ggplot(survey_summary %>% filter(lubridate::year(min_date) == 2024)) +
  geom_segment(aes(x = min_date, xend = max_date, 
                   y = reorder(site, min_date), yend = site), linewidth = 2) +
  scale_x_date(limits = as.Date(c("2024-06-01", "2024-09-14"))) +
  geom_vline(xintercept = start_date_2024, color = "red") +
  geom_vline(xintercept = end_date_2024, color = "red") +
  labs(x = "Date", y = "Site", title = "2024 survey effort") +
  theme_minimal() +
  theme(legend.position = "none")

start_date_2025 = earliest_common_date_2025 + 1
end_date_2025   = latest_common_date_2025 - 1

message("2025: start date ", start_date_2025, ", end date ", end_date_2025)
message("Days: ", end_date_2025 - start_date_2025)
ggplot(survey_summary %>% filter(lubridate::year(min_date) == 2025)) +
  geom_segment(aes(x = min_date, xend = max_date, 
                   y = reorder(site, min_date), yend = site), linewidth = 2) +
  scale_x_date(limits = as.Date(c("2025-06-01", "2025-09-14"))) +
  geom_vline(xintercept = start_date_2025, color = "red") +
  geom_vline(xintercept = end_date_2025, color = "red") +
  labs(x = "Date", y = "Site", title = "2025 survey effort") +
  theme_minimal() +
  theme(legend.position = "none")

# Filter for predictions within analysis window
message("Filtering for predictions within analysis window per season")
predictions_filtered <- predictions_metadata %>%
  mutate(date = as.Date(datetime_start)) %>%
  filter(
    (season == "2024" & date >= start_date_2024 & date <= end_date_2024) |
      (season == "2025" & date >= start_date_2025 & date <= end_date_2025)
  )

# # Count the number of prediction files (recordings) per unit-survey pairing
# message("Caching prediction file (i.e. recording) counts per unit survey to ", path_out_survey_file_counts)
# survey_file_counts = predictions_filtered %>%
#   group_by(season, site, date) %>%
#   summarise(n_prediction_files = n_distinct(file_path), .groups = "drop") %>%
#   mutate(season = factor(season), site = factor(site)) %>%
#   arrange(season, site, date)
# if (!dir.exists(path_out_dir)) dir.create(path_out_dir, recursive = TRUE)
# arrow::write_feather(survey_file_counts, path_out_survey_file_counts)

# Read and aggregate all predictions
message('Aggregating prediction data (this may take some time)')
n = nrow(predictions_filtered)
all_data_list = vector("list", n)  # preallocate list
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = n, clear = FALSE)
for (i in seq_len(n)) {
  file_path = predictions_filtered[[i, 'file_path']]
  rec_time_start = predictions_filtered[[i, 'datetime_start']]
  # Read the prediction history data
  prediction_data = read_csv(file_path, show_col_types = FALSE) %>% mutate(across(contains("Confidence"), ~ as.numeric(.)))
  if (nrow(prediction_data) > 0) {
    # Compute prediction time, add metadata, and drop irrelevant columns
    prediction_data = prediction_data %>% mutate(
      time     = as.POSIXct(rec_time_start + `Start (s)`),
      season   = predictions_filtered[[i, 'season']],
      yday_start   = predictions_filtered[[i, 'yday_start']],
      date = predictions_filtered[[i, 'date']],
      site     = predictions_filtered[[i, 'site']]
    )
  }
  prediction_data = prediction_data %>% select(-`Start (s)`, -`End (s)`, -`Scientific name`)
  all_data_list[[i]] = prediction_data
  pb$tick()
}
prediction_data = bind_rows(all_data_list)
prediction_data = prediction_data %>% clean_names()

message(nrow(prediction_data), " predictions aggregated")

# Write results to cache
dir.create(dirname(path_out_prediction_data), recursive = TRUE, showWarnings = FALSE)
write_feather(prediction_data, path_out_prediction_data)
message(crayon::green("Cached prediction data to", path_out_prediction_data))
