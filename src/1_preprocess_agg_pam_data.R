# 1_preprocess_agg_pam_data.R ==========================================================
# Derive and cache putative species detection histories and diversity metrics from aggregated PAM data

# INPUTS:
# Aggregated raw predictions from classifier and metadata
path_prediction_data = "data/cache/0_aggregate_raw_pam_data/prediction_data.feather"
path_species_list = "data/pam/species_list.txt"
path_avonet_traits = "data/traits/AVONET Supplementary dataset 1.xlsx"
# Naive thresholds
threshold_min_classifier_score = 0.5 # Naive classifier minimum confidence score threshold to assume binary presence/absence. # "For most false-positive models in our study, using a mid-range threshold of 0.50 or above generally yielded stable estimates." (Katsis et al. 2025)
threshold_min_detected_days = 2 # Minimum number of unique days detected to retain species detections at a site
# Classifier calibration (Platt scaling)
path_validation_data = "/Users/giojacuzzi/Library/CloudStorage/GoogleDrive-giojacuzzi@gmail.com/My Drive/Research/Projects/C5 - Riparian avian macroinvert monitoring/validations" # annotated Raven Pro selection tables organized by [common_name]/[file].selections.txt
use_platt_scaling = TRUE
display_plots = TRUE
tp_min_prob = 0.95 # Desired minimum probability of true positive

# OUTPUTS:
calibration_tag = ifelse(use_platt_scaling, "calibrated", "naive")
out_detect_hist_data = paste0("data/cache/1_preprocess_agg_pam_data/detections_", calibration_tag, "_", threshold_min_classifier_score, ".rds")

source("src/global.R")

# Load classifier prediction data
message("Loading classifier prediction data")
prediction_data = read_feather(path_prediction_data) %>% rename(site_id = site) %>% mutate(common_name = tolower(common_name))

# Helper functions --------------------------------------------------------------------------------

conf_to_logit = function(c) {
  c = min(max(c, 0.00001), 1.0 - 0.00001) # prevent undefined logit for extreme scores due to rounding error
  return(log(c / (1 - c)))
}

logit_to_conf = function(l) {
  return(1 / (1 + exp(-l)))
}

# Calculate species-specific probabilistic thresholds via Platt scaling -----------------------------

message("Loading validation data")

# Raven Pro selection tables. All corresponding audio segments must contain at least one label (unlabeled segments are ignored).
# Selection tables must have "Begin File" column, with files prefixed with the confidence score, e.g. "0.973_N_original_file_timestart_timeend.wav".
validation_files = list.files(path_validation_data, "\\.selections\\.txt$", recursive = TRUE, full.names = TRUE)
validation_files = data.frame(label_predicted = tolower(basename(dirname(validation_files))), path_validation = validation_files)

validation_data = validation_files %>%
  # Combine all validation selection tables
  rowwise() %>% mutate(data = list(
    read_tsv(path_validation, col_types = cols()) %>%
      janitor::clean_names() %>%
      select(begin_file, label) %>%
      rename(label_truth = label)
  )) %>%
  unnest(data) %>% ungroup() %>% select(-path_validation) %>%
  # Extract confidence score from audio segment file name
  mutate(confidence = as.numeric(str_extract(begin_file, "^[0-9.]+")))

if (use_platt_scaling) {
  calibration_results = data.frame()
  for (label in validation_files$label) {
    message("Calculating species-specific threshold for \"", label, "\"")
    
    # Obtain validation data for this label
    label_data = validation_data[validation_data$label_predicted == label, ]
    
    # Exclude "unknown" labels
    data = label_data %>% filter(label_truth != "unknown")
    
    # Create binary indicator for each unique audio file for true presence (aggregating all validations for each file)
    data = data %>% group_by(begin_file) %>%
      summarise(
        label_truth = if_else(any(label_truth == label), 1, 0),
        confidence = unique(confidence)
      )
    
    # Convert confidence scores to logit scale
    data = data %>% mutate(score = sapply(confidence, conf_to_logit))
    
    stopifnot(all(unique(data$label_truth) %in% c(0, 1)))
    
    # Calculate PR AUC and ROC AUC of classifier (from raw scores)
    auc_pr = pr.curve(scores.class0 = subset(data, label_truth == 1) %>% pull(score),
                      scores.class1 = subset(data, label_truth == 0) %>% pull(score))$auc.integral
    
    auc_roc = roc.curve(scores.class0 = subset(data, label_truth == 1) %>% pull(score),
                        scores.class1 = subset(data, label_truth == 0) %>% pull(score))$auc
    
    message("  PR-AUC ", round(auc_pr,3), ", ROC-AUC ", round(auc_roc,3))
    
    # Perform Platt scaling to determine threshold for desired probability of true positive
    threshold = Inf # default threshold is infinite (i.e. do not retain detections unless a species is validated)
    precision_threshold <- recall_threshold <- precision_tmin <- recall_tmin <- NA
    model_warning = FALSE
    tryCatch(
      {
        regression = glm(label_truth ~ score, data, family = binomial(link = "logit"))
        intercept   = as.numeric(coef(regression)[1])
        coefficient = as.numeric(coef(regression)[2])
        threshold_logit = (log(tp_min_prob / (1 - tp_min_prob)) - intercept) / coefficient # logit scale
        threshold       = logit_to_conf(threshold_logit) # confidence scale
        
        message("  ", round(threshold, 3), " threshold to achieve Pr(TP)>=", tp_min_prob)
        
        # Calculate estimated precision and recall at this threshold
        # Predicted positive/negative based on threshold
        calc_precision_recall = function(d, t_logit) {
          predicted = ifelse(d$score >= t_logit, 1, 0)
          TP = sum(predicted == 1 & d$label_truth == 1)
          FP = sum(predicted == 1 & d$label_truth == 0)
          FN = sum(predicted == 0 & d$label_truth == 1)
          return(data.frame(
            precision = TP / (TP + FP),
            recall    = TP / (TP + FN)
          ))
        }
        perf_t = calc_precision_recall(data, threshold_logit)
        precision_threshold  = perf_t$precision
        recall_threshold     = perf_t$recall
        message("  Performance at threshold ", round(threshold,3), ":\n  Precision ", round(precision_threshold,3), "\n  Recall ", round(recall_threshold,3))
        
        # Calculate at minimum threshold
        perf_tmin = calc_precision_recall(data, conf_to_logit(threshold_min_classifier_score))
        precision_tmin = perf_tmin$precision
        recall_tmin    = perf_tmin$recall
        message("  Performance at minimum threshold ", threshold_min_classifier_score, ":\n  Precision ", round(precision_tmin,3), "\n  Recall ", round(recall_tmin,3))
        
        if (display_plots) {
          
          ## Visualize precision and recall performance as a function of score
          data_sorted = data %>% mutate(score = logit_to_conf(score)) %>% arrange(desc(score))
          n_pos = sum(data_sorted$label_truth == 1)
          data_sorted = data_sorted %>% mutate(
            tp = cumsum(label_truth == 1),
            fp = cumsum(label_truth == 0),
            recall = tp / n_pos,
            precision = ifelse(tp + fp == 0, 1, tp / (tp + fp))  # handle division by zero
          )
          plt = ggplot(data_sorted, aes(x = score)) +
            geom_line(aes(y = precision, color = "Precision"), linewidth = 1) +
            geom_line(aes(y = recall, color = "Recall"), linewidth = 1) +
            labs(title = paste0(label), x = "Score", y = "Performance", color = "Metric"); print(plt)
          
          ## Visualize the logistic regression and data
          x_range = seq(min(data$score), max(data$score), length.out = 100)
          
          # Predict fitted values and calculate confidence intervals on probability scale
          pred = predict(regression, newdata = data.frame(score = x_range), type = "link", se.fit = TRUE)
          
          # Calculate confidence intervals on probability scale
          z = qnorm(0.975)  # 95% CI
          upper  = logit_to_conf(pred$fit + z * pred$se.fit)
          lower  = logit_to_conf(pred$fit - z * pred$se.fit)
          y_pred = logit_to_conf(pred$fit)
          
          regression_df = data.frame(score = x_range, prob = y_pred, lower = lower, upper = upper, warning = model_warning)
          plt = ggplot() +
            geom_vline(xintercept = 0, color = "gray", linetype = "solid", linewidth = 0.5) +
            geom_point(data = data, aes(x = score, y = label_truth), shape = 1, alpha = 0.5) +
            geom_hline(yintercept = tp_min_prob, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_vline(xintercept = threshold_logit, color = "black", linetype = "dashed", linewidth = 0.5) +
            geom_ribbon(data = regression_df, aes(x = score, ymin = lower, ymax = upper, fill = warning), alpha = 0.2) +
            geom_line(data = regression_df, aes(x = score, y = prob, color = warning), linewidth = 1) +
            scale_color_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            scale_fill_manual(values = c("FALSE" = "blue", "TRUE" = "red")) +
            labs(
              x = "Score (logit)", y = "True positive probability",
              title = paste0(label), subtitle = paste0("Threshold p(TP) ≥ ", tp_min_prob, " = ", round(threshold, 3))
            ); print(plt)
        }
      },
      warning = function(w) {
        model_warning <<- TRUE
        message(crayon::yellow("WARNING", conditionMessage(w)))
      },
      error = function(e) {
        stop(crayon::red("ERROR", conditionMessage(e)))
      }
    )
    
    # Enforce minimum threshold for classes with no negatives in validation data
    if (label %in% c(
      "swainson's thrush",
      "black-throated gray warbler",
      "hammond's flycatcher",
      "hutton's vireo",
      "northern flicker",
      "pacific wren",
      "pacific-slope flycatcher",
      "violet-green swallow"
    )) {
      threshold = threshold_min_classifier_score
    }
    
    calibration_results = rbind(calibration_results, data.frame(
      common_name   = label,
      n_pos         = sum(data$label_truth == 1),
      n_neg         = sum(data$label_truth == 0),
      auc_pr        = auc_pr,
      auc_roc       = auc_roc,
      warning       = model_warning,
      tp_min_prob   = tp_min_prob,
      threshold     = threshold,
      threshold_min = threshold_min_classifier_score,
      precision_threshold  = precision_threshold,
      recall_threshold     = recall_threshold,
      precision_tmin = precision_tmin,
      recall_tmin    = recall_tmin
    ))
  }
  
  # Clean results
  calibration_results = calibration_results %>%
    mutate(across(everything(), ~ ifelse(is.nan(.) | . == "NaN", NA, .)))
  message("Calibration results:")
  print(calibration_results)
}

# Apply thresholds ----------------------------------------------------------------------------------

species_names = read_lines(path_species_list) %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>%
  mutate(scientific_name = tolower(scientific_name), common_name = tolower(common_name))

# Manually exclude any species determined not present at any site
# TODO
detections = prediction_data %>% filter(!common_name %in% c(
  "fox sparrow",
  "lincoln's sparrow",
  "barn swallow",
  "cliff swallow",
  "house wren",
  "marsh wren",
  "northern rough-winged swallow",
  "purple martin",
  "townsend's warbler",
  "tree swallow",
  "western meadowlark",
  "brewer's blackbird",
  "marbled murrelet",
  "pied-billed grebe",
  "ruddy duck",
  "bufflehead",
  "american pipit",
  "ruby-crowned kinglet"
))

species_names = species_names %>% filter(common_name %in% detections$common_name)

if (use_platt_scaling) {
  
  # Obtain putative detections with calibrated thresholds
  message("Obtaining putative detections with calibrated / minimum thresholds")
  
  species_thresholds = species_names %>% mutate(threshold = threshold_min_classifier_score)
  
  # Calibrated threshold for common_name in calibration_results
  species_thresholds = species_thresholds %>% left_join(
      calibration_results %>% select(common_name, threshold_calibrated = threshold),
      by = "common_name"
    ) %>% mutate(threshold = coalesce(threshold_calibrated, threshold)) %>% select(-threshold_calibrated)
  
  # Manually exclude detections of specific species for which completely manual validation is necessary via Inf threshold
  species_thresholds[species_thresholds$common_name %in% c(
    "american dipper",
    "common yellowthroat",
    "macgillivray's warbler",
    "orange-crowned warbler",
    "red-eyed vireo",
    "cassin's vireo",
    "red-breasted sapsucker",
    "vaux's swift",
    "merlin",
    "spotted sandpiper",
    "northern pygmy-owl"
  ), "threshold"] = Inf
  
  # Enforce the minimum threshold
  species_thresholds = species_thresholds %>% mutate(threshold = pmax(threshold, threshold_min_classifier_score))
  message("Final species thresholds:")
  print(species_thresholds, n = nrow(species_names))
  detections = prediction_data %>% left_join(species_thresholds %>% select(common_name, threshold), by = "common_name") %>%
    filter(confidence >= threshold) %>% select(-threshold)

} else {
  
  # Obtain putative detections with naive minimum threshold
  message("Obtaining putative detections with naive minimum threshold ", threshold_min_classifier_score)
  detections = prediction_data %>% filter(confidence >= threshold_min_classifier_score)
  
}

# Remove putative detections for which we have ground-truth validation data
validation_data = validation_data %>% separate(begin_file,
                                               into = c("confidence_str", "x2", "site_id", "date_str", "time_start", "time_start_s", "time_end_s"),
                                               sep = "_",
                                               remove = FALSE) %>%
  mutate(
    confidence = as.numeric(confidence_str),
    site_id = as.character(site_id),
    date = ymd(date_str),
    season = as.character(year(date)),
    time_start_dt = ymd_hms(paste(date_str, time_start)),
    time_start_s = as.numeric(str_remove(time_start_s, "s")),
    time = time_start_dt + seconds(time_start_s)
  ) %>%
  select(label_truth, label_predicted, confidence, begin_file, time, season, date, site_id) %>%
  rename(common_name = label_predicted, file = begin_file) %>% arrange(desc(confidence))

detections = detections %>% anti_join(validation_data, by = c("site_id", "common_name", "time"))

# Add presence detections for which we have ground-truth validation data
v = validation_data %>% filter(!label_truth %in% c("unknown", "not_focal", "not_target"), label_truth == common_name)
v$confidence = 2.0 # 2.0 confidence value indicates manual validation

detections = detections %>%
  bind_rows(v %>%
      mutate(
        yday_start = NA_real_
      ) %>%
      select(common_name, confidence, file, time, season, yday_start, date, site_id)
  )

# Only retain detections for species detected a minimum number of days at a site
message("Discarding detections for species that were detected during less than ", threshold_min_detected_days, " surveys (days) at a site")
detections = detections %>%
  group_by(site_id, common_name) %>%
  summarise(n_surveys = n_distinct(date), .groups = "drop") %>%
  filter(n_surveys >= threshold_min_detected_days) %>%
  inner_join(detections, by = c("site_id", "common_name"))

# Determine detected species and detections per site ----------------------------------------

species_names = species_names %>%
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

# NOTE: Shared survey start/end dates were previously determined in 0_agg_raw_pam_data.R
start_date_2024 = min(detections$date[detections$season == "2024"])
end_date_2024   = max(detections$date[detections$season == "2024"])
start_date_2025 = min(detections$date[detections$season == "2025"])
end_date_2025   = max(detections$date[detections$season == "2025"])

# Determine survey numbers by year
clean_detections = detections %>%
  mutate(
    survey_num = case_when(
      year(date) == 2024 ~ as.integer(difftime(date, start_date_2024, units = "days")) + 1,
      year(date) == 2025 ~ as.integer(difftime(date, start_date_2025, units = "days")) + 1,
      TRUE ~ NA
    )
  ) %>% arrange(site_id, survey_num)

# Make a complete grid of site × species x surveys
template_all_site_species_surveys = expand.grid(
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
  mutate(scientific_name = tolower(scientific_name)) %>%
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
  threshold_min_classifier_score = threshold_min_classifier_score,
  threshold_min_detected_days = threshold_min_detected_days,
  long = species_site_survey_long,
  wide = species_site_survey_wide
)

if (!dir.exists(dirname(out_detect_hist_data))) dir.create(dirname(out_detect_hist_data), recursive = TRUE)
saveRDS(detect_hist_data, out_detect_hist_data)
message(crayon::green("Cached", out_detect_hist_data))
