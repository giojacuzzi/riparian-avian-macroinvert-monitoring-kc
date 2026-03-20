# 6_predictions.R ======================================================================
# Use ICLUS impervious surface projections to predict future changes in B-IBI and predator
# richness across sites
#
# Input
path_sem = "data/cache/4_sem/sem_diet.rds"
path_projections_reach = "data/cache/3_calculate_vars/imp_projections_550m.rds"
path_projections_watershed = "data/cache/3_calculate_vars/imp_projections_5000m.rds"
use_2023_observed_as_baseline = FALSE
# Deltas between 2023 observed response variables (canopy, bibi, richness) and 2100 predictions
# may be positive OR negative because it includes residual error from component model fit
# Canopy increasing reflects the presence of negative residuals in the 2023 data, not that the model actually predicts a gain in cover.
#
# Output
out_dir = "data/cache/6_prediction"

# Load projections
projections_reach     = readRDS(path_projections_reach)
projections_watershed = readRDS(path_projections_watershed)

source("src/global.R")

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Site-specific % change in ICLUS impervious surface coverage from 2020 was calculated for each decade and
# multiplied with "ground truth" 2023 NLCD measurements to estimate future coverage at each site.
projections_reach %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))
projections_watershed %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))

# Load sem for prediction and extract component models
sem = readRDS(path_sem)
m_canopy = sem[[1]]
m_bibi   = sem[[2]]
m_rich   = sem[[3]]
data     = sem[[4]]

# Only retain sites with projections
data = data %>% filter(site_id %in% projections_reach$site_id) %>% arrange(site_id)

# For each scenario
scenario_data_pred = list()
for (s in c("a1", "a2", "b1", "b2")) {

  message("Projecting conditions for scenario ", s)
  
  # Get projected imp_reach (local) and imp_basin (landscape)
  
  data_pred = projections_reach %>% filter(scenario == s) %>% mutate(imp_reach = sum_proportion)
  if (use_2023_observed_as_baseline) {
    data_pred = data_pred %>% filter(year != "2023")
  }
  
  # Note that percent_change_sum_prop is actually a multiplicative factor of change relative to the 2020 baseline
  # Reach-level ISC is projected to increase fourfold at the most impacted stream (mean increase by a factor of 1.63 across all streams)
  summary(projections_reach %>% filter(scenario == s, year == "2100") %>% pull(percent_change_sum_prop))
  
  # Relative to baseline conditions, both reach- and catchment-level ISC are expected to increase on average by roughly an additional 10% of surface area cover by the year 2100.
  summary(projections_reach %>% filter(scenario == s, year == "2100") %>%
            mutate(
              baseline_sum_proportion = sum_proportion / percent_change_sum_prop,
              delta_sum_proportion = sum_proportion - baseline_sum_proportion
            ) %>% pull(delta_sum_proportion))
  
  p_proj_imp_reach = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = imp_reach, group = site_id)) +
    # geom_line(color = "black", alpha = 0.2) +
    geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("gray20", 0.2)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(expand = c(0,0)) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", linewidth = 1.2) +
    # geom_text(
    #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
    #   aes(label = site_id), hjust = -0.1, size = 2, color = "gray60"
    # ) +
    labs(x = "Year", y = "", title = "Impervious % (reach)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  data_pred = data_pred %>% left_join(projections_watershed %>%
                                        rename(imp_basin = sum_proportion) %>%
                                        select(-percent_change_sum_prop),
                                      by = c("scenario", "site_id", "year"))
  
  # Note that percent_change_sum_prop is actually a multiplicative factor of change relative to the 2020 baseline
  # ...while catchment-level ISC will increase twofold (mean 1.58).
  summary(projections_watershed %>% filter(scenario == s, year == "2100") %>% pull(percent_change_sum_prop))
  
  summary(projections_watershed %>% filter(scenario == s, year == "2100") %>%
            mutate(
              baseline_sum_proportion = sum_proportion / percent_change_sum_prop,
              delta_sum_proportion = sum_proportion - baseline_sum_proportion
            ) %>% pull(delta_sum_proportion))
  
  p_proj_imp_basin = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = imp_basin, group = site_id)) +
    # geom_line(color = "black", alpha = 0.2) +
    geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("gray20", 0.2)) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_discrete(expand = c(0,0)) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", linewidth = 1.2) +
    # geom_text(
    #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
    #   aes(label = site_id), hjust = -0.1, size = 2, color = "gray60"
    # ) +
    labs(x = "Year", y = "", title = "Impervious % (watershed)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Predict canopy_reach from imp_reach
  data_pred = data_pred %>% mutate(
    tcc_reach = predict(m_canopy, newdata = data.frame(imp_reach = imp_reach))
  )
  
  p_proj_canopy_reach = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = plogis(tcc_reach), group = site_id)) +
    # geom_line(color = "black", alpha = 0.2) +
    geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("gray20", 0.2)) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_continuous(limits = c(0,1)) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "forestgreen", linewidth = 1.2) +
    # geom_text(
    #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
    #   aes(label = site_id), hjust = -0.1, size = 2, color = "gray60"
    # ) +
    labs(x = "Year", y = "", title = "Canopy % (reach)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Plot time series environmental projections
  fig_S3A = p_proj_imp_basin + p_proj_imp_reach + p_proj_canopy_reach
  print(fig_S3A)
  
  ggsave(paste0(out_dir, "/fig_S3A_", s, ".pdf"), fig_S3A, width = 7, height = 3)
  
  # TCC change relative to baseline
  baseline_imp_reach = projections_reach %>% filter(scenario == s, year == "2100") %>%
    mutate(
      baseline_sum_proportion = sum_proportion / percent_change_sum_prop
    )
  compare_tcc_reach = data.frame(
    site_id = baseline_imp_reach$site_id,
    baseline_tcc_reach = plogis(predict(m_canopy, newdata = data.frame(imp_reach = baseline_imp_reach$baseline_sum_proportion)))
  )
  compare_tcc_reach = compare_tcc_reach %>%
    left_join(data_pred %>% filter(year == "2100") %>% select(site_id, tcc_reach) %>% mutate(tcc_reach_2100 = plogis(tcc_reach))) %>%
    mutate(delta_tcc_reach = round(tcc_reach_2100 - baseline_tcc_reach, 3),
           change_factor = tcc_reach_2100 / baseline_tcc_reach)
  # Consequently, the total amount of TCC within the reach was predicted to decrease by up to 22%
  summary(compare_tcc_reach$delta_tcc_reach)
  summary(compare_tcc_reach$change_factor)
  
  # Predict bibi from predicted canopy_reach and projected imp_basin
  data_pred = data_pred %>% mutate(
    bibi = predict(m_bibi, newdata = data.frame(tcc_reach = tcc_reach, imp_basin = imp_basin))
  )
  
  p_pred_bibi = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = plogis(bibi) * 100, group = site_id)) +
    # geom_line(color = "black", alpha = 0.2) +
    geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "royalblue", linewidth = 1.2) +
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100), expand = c(0, 0)) +
    # geom_text(
    #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
    #   aes(label = site_id), hjust = -0.1, size = 3, color = "gray60"
    # ) +
    scale_x_discrete(expand = c(0, 0)) + 
    labs(x = "", y = "Projected B-IBI") +
    theme(axis.text.x = element_blank(), #element_text(angle = 45, hjust = 1)
          axis.ticks.x = element_blank(),
          plot.margin = margin(0,0,0,0)
          )
  
  # Predict rich_predator from predicted bibi, predicted canopy_reach, and projected imp_reach
  data_pred = data_pred %>% mutate(rich_predator = predict(m_rich, newdata = data.frame(
    bibi = bibi,
    tcc_reach = tcc_reach,
    imp_reach = imp_reach
  ), type = "response"))
  
  p_pred_rich_predator = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = rich_predator, group = site_id)) +
    # geom_line(color = "black", alpha = 0.2) +
    geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
    stat_summary(aes(group = 1), fun = mean, geom = "line", color = "orange", size = 1.2) +
    # geom_text(
    #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
    #   aes(label = site_id), hjust = -0.1, size = 3, color = "gray60"
    # ) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ceiling(max(data_pred$rich_predator))), breaks = seq(0,12,2)) + 
    # scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 0, x, ""), expand = c(0, 0)) + 
    scale_x_discrete(expand = c(0,0)) +
    labs(x = "Year", y = "Projected predator richness") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(0,0,0,0))
  
  # Plot time series predictions for responses
  # print(p_pred_bibi + p_pred_rich_predator)
  
  # Calculate deltas from 2023 to 2100
  calc_delta = function(d, var, year_start = "2023", year_end = "2100") {
    delta_name = paste0("delta_", var)
    delta_pct_name = paste0("deltapcnt_", var)
    
    d %>%
      filter(year %in% c(year_start, year_end)) %>%
      select(site_id, year, all_of(var)) %>%
      pivot_wider(
        names_from = year,
        values_from = all_of(var),
        names_prefix = paste0(var, "_")
      ) %>%
      mutate(
        !!delta_name := .[[paste0(var, "_", year_end)]] - .[[paste0(var, "_", year_start)]],
        !!delta_pct_name := 100 * (.[[paste0(var, "_", year_end)]] - .[[paste0(var, "_", year_start)]]) /
          .[[paste0(var, "_", year_start)]]
      )
  }
  
  if (use_2023_observed_as_baseline) {
    # Add observed 2023 data to data_pred
    observed_2023 = data %>% mutate(
      year = "2023",       # add year as character to match data_pred
      scenario = s         # same scenario
    ) %>%
      select(site_id, scenario, year, imp_basin, imp_reach, canopy_reach, bibi, rich_predator)
    observed_2023 = observed_2023 %>% filter(site_id %in% unique(data_pred$site_id))
    
    data_pred = bind_rows(data_pred, observed_2023)
  }
  
  delta_vars = c("imp_basin", "imp_reach", "tcc_reach", "bibi", "rich_predator") %>%
    map(~ calc_delta(data_pred, .x)) %>% reduce(full_join, by = "site_id")
  
  summary(delta_vars$delta_imp_basin)
  summary(delta_vars$delta_imp_reach)
  summary(delta_vars$delta_tcc_reach)
  summary(delta_vars$delta_bibi)
  summary(delta_vars$delta_rich_predator)
  summary(delta_vars$deltapcnt_rich_predator)
  
  p_bibi_delta = ggplot(
    delta_vars %>%
      mutate(site_id = factor(site_id, levels = delta_vars$site_id[order(delta_vars$bibi_2023, decreasing = TRUE)])) %>%
      pivot_longer(c(bibi_2023, bibi_2100), names_to = "year", values_to = "bibi"),
    aes(x = site_id, y = plogis(bibi) * 100, group = site_id)) +
    geom_line(color = "gray10") +
    geom_point(aes(color = year)) +
    scale_color_manual(values = c("bibi_2023" = "gray10", "bibi_2100" = "royalblue")) +
    scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100), expand = c(0,0)) +
    labs(x = "", y = "", color = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0,0,0,0),
          panel.grid.major.y = element_line(color = "grey90", linewidth = 0.4),
          legend.position = "none")
  
  p_rich_predator_delta = ggplot(delta_vars %>%
           mutate(site_id = factor(site_id, levels = delta_vars$site_id[order(delta_vars$rich_predator_2023, decreasing = TRUE)])) %>%
           pivot_longer(c(rich_predator_2023, rich_predator_2100), names_to = "year", values_to = "rich_predator"),
    aes(x = site_id, y = rich_predator, group = site_id)) +
    geom_line() +
    geom_point(aes(color = year)) +
    scale_color_manual(values = c("rich_predator_2023" = "gray10", "rich_predator_2100" = "orange")) +
    scale_y_continuous(breaks = seq(0, 14, by = 1), limits = c(0, ceiling(max(data_pred$rich_predator))), expand = c(0,0)) +
    labs(x = "", y = "", color = "") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = margin(0,0,0,0),
          legend.position = "bottom"
          )
  
  # Plot time series and delta projections for responses
  p_proj_responses = (p_pred_bibi + p_bibi_delta) / (p_pred_rich_predator + p_rich_predator_delta)
  fig_5 = p_proj_responses + plot_annotation(tag_levels = "A")
  print(fig_5)
  ggsave(paste0(out_dir, "/fig_5_", s, ".pdf"), fig_5, width = 8.2, height = 7.5)

  scenario_data_pred[[s]] = data_pred
}

# Local projections across scenarios
scenario_isc_local = ggplot(projections_reach %>% filter(year != "2023"), aes(x = year, y = sum_proportion * 100)) +
  stat_summary(
    aes(group = str_to_title(scenario), color = str_to_title(scenario)),
    fun = mean, geom = "line", linewidth = 1,
  ) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "", y = "Projected mean ISC local", color = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); print(scenario_isc_local)

# Landscape projections across scenarios
scenario_isc_land = ggplot(projections_watershed %>% filter(year != "2023"), aes(x = year, y = sum_proportion * 100)) +
  stat_summary(
    aes(group = str_to_title(scenario), color = str_to_title(scenario)),
    fun = mean, geom = "line", linewidth = 1,
  ) +
  scale_color_viridis_d() +
  scale_y_continuous(limits = c(0,100)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x = "", y = "Projected mean ISC landscape", color = "Scenario") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); print(scenario_isc_land)


all_scenarios = imap_dfr(
  scenario_data_pred,
  ~ .x %>%
    mutate(scenario = .y)   # .y is the list name (e.g. "a1")
)

scenario_canopy = ggplot(all_scenarios %>% filter(year != "2023"),
                       aes(x = year, y = plogis(tcc_reach) * 100, color = str_to_title(scenario), group = str_to_title(scenario))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "", y = "Projected mean TCC local", color = "Scenario") +
  scale_color_viridis_d(option = "viridis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); print(scenario_canopy)

scenario_bibi = ggplot(all_scenarios %>% filter(year != "2023"),
  aes(x = year, y = plogis(bibi) * 100, color = str_to_title(scenario), group = str_to_title(scenario))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "", y = "Projected mean B-IBI", color = "Scenario") +
  scale_color_viridis_d(option = "viridis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); print(scenario_bibi)

scenario_rich = ggplot(all_scenarios %>% filter(year != "2023"),
                       aes(x = year, y = rich_predator, color = str_to_title(scenario), group = str_to_title(scenario))) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  scale_y_continuous(breaks = seq(0,12,2), limits = c(0, ceiling(max(data_pred$rich_predator))), expand = c(0, 0)) +
  scale_x_discrete(expand = c(0, 0)) +
  labs(x = "", y = "Projected mean predator richness", color = "Scenario") +
  scale_color_viridis_d(option = "viridis") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); print(scenario_rich)

fig_S3B = (
  scenario_isc_local + theme(legend.position = "none") +
  scenario_isc_land + theme(legend.position = "none") +
  scenario_canopy + theme(legend.position = "none")
) / (
  scenario_bibi + theme(legend.position = "none") +
  scenario_rich
); print(fig_S3B)

ggsave(paste0(out_dir, "/fig_S3B.pdf"), fig_S3B, width = 7, height = 6)

# Landscape watershed-only restoration (e.g. integrated watershed management plans) -----------------------------
# Assuming imp_reach and canopy_reach values cannot be changed (e.g. local-scale restoration is
# not feasible due to private land ownership), estimate per site what value of imp_basin
# is needed to maintain the same amount of richness as in 2023

data_pred = scenario_data_pred[["b1"]]

# Filter to the 2023 richness and projected imp_basin in 2100
sites_2100_basin = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, scenario, year, imp_reach, tcc_reach, imp_basin, bibi, rich_predator) %>%
  pivot_wider(names_from = year, values_from = c(imp_reach, tcc_reach, imp_basin, bibi, rich_predator), names_sep = "_")

# Function that returns richness given imp_reach and site-specific values
richness_for_site_basin = function(imp_basin_try, imp_reach_2100, m_canopy, m_bibi, m_rich) {
  # predict canopy
  tcc_try = predict(m_canopy, newdata = data.frame(
    imp_reach = imp_reach_2100
  ))
  # predict bibi
  bibi_try = predict(m_bibi, newdata = data.frame(
    tcc_reach = tcc_try,
    imp_basin = imp_basin_try
  ))
  # predict richness
  predict(m_rich, newdata = data.frame(
    bibi = bibi_try,
    tcc_reach = tcc_try,
    imp_reach = imp_reach_2100
  ), type = "response")
}

# Solve for required imp_basin target
sites_2100_basin = sites_2100_basin %>% rowwise() %>%
  mutate(
    imp_basin_target = tryCatch({
      uniroot(
        function(x)
          richness_for_site_basin(imp_basin_try = x, imp_reach_2100 = imp_reach_2100, m_canopy = m_canopy, m_bibi = m_bibi, m_rich = m_rich) - rich_predator_2023,
        lower = 0, upper = 1
      )$root
    },
    error = function(e) NA)
  ) %>% ungroup()

# Compute deltas
sites_2100_basin = sites_2100_basin %>% mutate(imp_basin_target_delta = imp_basin_2023 - imp_basin_target)

# Summary statistics
summary(sites_2100_basin %>% pull(imp_basin_target_delta))

# NOTE: NA's suggest that no amount of basin level restoration can maintain richness levels
# in the face of how degraded the stream reach becomes by 2100
sites_2100_basin = sites_2100_basin %>% mutate(restoration_possible = ifelse(is.na(imp_basin_target), "No", "Yes"))
site_id_no_watershed_restoration_possible = sites_2100_basin %>% filter(restoration_possible == "No") %>% pull(site_id)

# Plot
sites_2100_basin = sites_2100_basin %>% arrange(imp_basin_2023) %>% mutate(site_id = factor(site_id, levels = site_id))
sites_2100_basin_long = sites_2100_basin %>%
  pivot_longer(
    cols = c(imp_basin_2023, imp_basin_target, imp_basin_2100),
    names_to = "type",
    values_to = "imp_basin"
  ) %>%
  mutate(type = factor(
    type,
    levels = c("imp_basin_2023", "imp_basin_target", "imp_basin_2100"),
    labels = c("Current (2023)", "Restoration to maintain richness", "No action (2100)")
  ))

p_restore_basin_imp = ggplot(sites_2100_basin_long, aes(y = site_id)) +
  geom_segment(data = sites_2100_basin, aes(x = imp_basin_2023, xend = imp_basin_2100, y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100_basin, aes(x = imp_basin_2023, xend = imp_basin_target, y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = imp_basin, color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Watershed impervious %", y = "", color = "", shape = "") +
  theme(legend.position = "none")

# Local reach-only restoration (e.g. public lands, private land easements, land trusts) ---------------------
# Assuming imp_basin values are inevitable in 2100, estimate per site what value of imp_reach
# is needed to maintain the same amount of richness as in 2023

# Filter to the 2023 richness and projected imp_basin in 2100
sites_2100_reach = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, scenario, year, imp_reach, tcc_reach, imp_basin, bibi, rich_predator) %>%
  pivot_wider(names_from = year, values_from = c(imp_reach, tcc_reach, imp_basin, bibi, rich_predator), names_sep = "_")

# Function that returns richness given imp_reach and site-specific values
richness_for_site_reach = function(imp_reach_try, imp_basin_2100, m_canopy, m_bibi, m_rich) {
  # predict canopy
  tcc_try = predict(m_canopy, newdata = data.frame(
    imp_reach = imp_reach_try
  ))
  # predict bibi
  bibi_try = predict(m_bibi, newdata = data.frame(
    tcc_reach = tcc_try,
    imp_basin = imp_basin_2100
  ))
  # predict richness
  predict(m_rich, newdata = data.frame(
    bibi = bibi_try,
    tcc_reach = tcc_try,
    imp_reach = imp_reach_try
  ), type = "response")
}

# Solve for required imp_reach target
sites_2100_reach = sites_2100_reach %>% rowwise() %>%
  mutate(
    imp_reach_target = tryCatch({ 
      uniroot(
        function(x) richness_for_site_reach(x, imp_basin_2100, m_canopy, m_bibi, m_rich) - rich_predator_2023,
        lower = 0, upper = 1
      )$root
    }, error = function(e) NA)
  ) %>% ungroup()

sites_2100_reach = sites_2100_reach %>% mutate(imp_reach_target_delta = imp_reach_2023 - imp_reach_target)
sites_2100_reach = sites_2100_reach %>% mutate(restoration_possible = ifelse(is.na(imp_reach_target), "No", "Yes"))
site_id_no_reach_restoration_possible = sites_2100_reach %>% filter(restoration_possible == "No") %>% pull(site_id)
sites_2100_reach = sites_2100_reach %>% arrange(imp_reach_2023) %>%
  mutate(site_id = factor(site_id, levels = site_id))

sites_2100_long = sites_2100_reach %>%
  pivot_longer(cols = c(imp_reach_2023, imp_reach_target, imp_reach_2100), names_to = "type", values_to = "imp_reach") %>%
  mutate(type = factor(type,
                       levels = c("imp_reach_2023", "imp_reach_target", "imp_reach_2100"),
                       labels = c("Current (2023)", "Restoration to maintain richness", "No action (2100)")))

p_restore_reach_imp = ggplot(sites_2100_long, aes(y = site_id)) +
  geom_segment(data = sites_2100_reach, aes(x = imp_reach_2023, xend = imp_reach_2100, y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100_reach, aes(x = imp_reach_2023, xend = imp_reach_target, y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = imp_reach, color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Reach impervious %", y = "Site", color = "", shape = "") +
  theme(legend.position = "bottom")

summary(sites_2100_reach %>% pull(imp_reach_target_delta))

# Align plots to reach sites
site_order_reach = levels(sites_2100_reach$site_id)

# Apply to basin data (both wide and long)
sites_2100_basin_aligned = sites_2100_basin
sites_2100_basin_aligned$site_id <- factor(sites_2100_basin_aligned$site_id, levels = site_order_reach)
sites_2100_basin_aligned_long = sites_2100_basin_long
sites_2100_basin_aligned_long$site_id <- factor(sites_2100_basin_aligned_long$site_id, levels = site_order_reach)

p_restore_basin_aligned_imp = ggplot(sites_2100_basin_aligned_long, aes(y = site_id)) +
  geom_segment(data = sites_2100_basin_aligned,
               aes(x = imp_basin_2023, xend = imp_basin_2100,
                   y = site_id, yend = site_id,
                   color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100_basin_aligned,
               aes(x = imp_basin_2023, xend = imp_basin_target,
                   y = site_id, yend = site_id),
               color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = imp_basin, color = type, shape = type)) +
  scale_color_manual(values = c(
    "Current (2023)" = "black",
    "Restoration to maintain richness" = "forestgreen",
    "No action (2100)" = "grey40"
  )) +
  scale_shape_manual(values = c(
    "Current (2023)" = 19,
    "Restoration to maintain richness" = 19,
    "No action (2100)" = 19
  )) +
  labs(x = "Watershed impervious %", y = "", color = "", shape = "") +
  theme(legend.position = "none")

print(p_restore_reach_imp + p_restore_basin_aligned_imp)

# Align plots to basin sites
site_order_basin <- levels(sites_2100_basin$site_id)
sites_2100_reach_aligned <- sites_2100_reach
sites_2100_reach_aligned$site_id <- factor(
  sites_2100_reach_aligned$site_id,
  levels = site_order_basin
)
sites_2100_long_aligned <- sites_2100_long
sites_2100_long_aligned$site_id <- factor(
  sites_2100_long_aligned$site_id,
  levels = site_order_basin
)
p_restore_reach_aligned_imp =
  ggplot(sites_2100_long_aligned, aes(y = site_id)) +
  
  geom_segment(
    data = sites_2100_reach_aligned,
    aes(x = imp_reach_2023,
        xend = imp_reach_2100,
        y = site_id,
        yend = site_id,
        color = restoration_possible)
  ) +
  
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  
  geom_segment(
    data = sites_2100_reach_aligned,
    aes(x = imp_reach_2023,
        xend = imp_reach_target,
        y = site_id,
        yend = site_id),
    color = "forestgreen",
    linetype = "dotted"
  ) +
  
  geom_point(aes(x = imp_reach, color = type, shape = type)) +
  
  scale_color_manual(values = c(
    "Current (2023)" = "black",
    "Restoration to maintain richness" = "forestgreen",
    "No action (2100)" = "grey40"
  )) +
  
  scale_shape_manual(values = c(
    "Current (2023)" = 19,
    "Restoration to maintain richness" = 19,
    "No action (2100)" = 19
  )) +
  
  labs(x = "Reach impervious %", y = "Site", color = "", shape = "") +
  theme(legend.position = "bottom")

print(p_restore_reach_aligned_imp + p_restore_basin_imp)

# Find canopy reach target values
sites_2100_reach = sites_2100_reach %>% rowwise() %>%
  mutate(
    tcc_reach_target = predict(m_canopy, newdata = data.frame(
      imp_reach = imp_reach_target
    ))
  ) %>% ungroup()
sites_2100_reach = sites_2100_reach %>% mutate(tcc_reach_target_delta = tcc_reach_2023 - tcc_reach_target)
sites_2100_reach = sites_2100_reach %>% arrange(imp_reach_2023) %>%
  mutate(site_id = factor(site_id, levels = site_id))

sites_2100_long = sites_2100_reach %>%
  pivot_longer(cols = c(tcc_reach_2023, tcc_reach_target, tcc_reach_2100), names_to = "type", values_to = "tcc_reach") %>%
  mutate(type = factor(type,
                       levels = c("tcc_reach_2023", "tcc_reach_target", "tcc_reach_2100"),
                       labels = c("Current (2023)", "Restoration to maintain richness", "No action (2100)")))

p_restore_reach_tcc = ggplot(sites_2100_long, aes(y = site_id)) +
  geom_segment(data = sites_2100_reach, aes(x = plogis(tcc_reach_2023), xend = plogis(tcc_reach_2100), y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100_reach, aes(x = plogis(tcc_reach_2023), xend = plogis(tcc_reach_target), y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = plogis(tcc_reach), color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Reach canopy %", y = "Site", color = "", shape = "") +
  theme(legend.position = "bottom")

# Plot both
p_restore_imp = p_restore_reach_imp + p_restore_basin_imp
print(p_restore_imp + plot_annotation(tag_levels = "A"))

# Both reach (local) and watershed (landscape) restoration -----------------------------------------------------------------

site_restorations = data.frame(
  site_id = data$site_id,
  cannot_restore_reach_only = data$site_id %in% site_id_no_reach_restoration_possible,
  cannot_restore_basin_only = data$site_id %in% site_id_no_watershed_restoration_possible
)
site_restorations$cannot_restore_either_only = site_restorations$cannot_restore_reach_only & site_restorations$cannot_restore_basin_only

site_restorations$scales_required = ifelse(site_restorations$cannot_restore_either_only, "Multi-scale", "Single scale")

site_restorations = site_restorations %>%
  mutate(
    needed_restoration = factor(
      case_when(
        cannot_restore_reach_only & cannot_restore_basin_only ~ "Both",
        cannot_restore_reach_only ~ "Catchment",
        cannot_restore_basin_only ~ "Reach",
        TRUE ~ "Either"
      ),
      levels = c("Reach", "Catchment", "Either", "Both")
    )
  )

site_restorations = site_restorations %>% filter(site_id %in% delta_vars$site_id)
site_restorations = site_restorations %>% left_join(delta_vars %>% select(site_id, imp_basin_2023, imp_basin_2100, imp_reach_2023, imp_reach_2100), by = "site_id")

site_restorations = site_restorations %>%
  mutate(
    delta_imp_basin = imp_basin_2023 - imp_basin_2100,
    delta_imp_reach = imp_reach_2023 - imp_reach_2100,
    delta_imp_max = pmax(abs(delta_imp_basin), abs(delta_imp_reach))
  )

## Sensitivity analysis richness ~ imp
pred_grid = expand.grid(
  imp_reach = seq(0, 1.0, length.out = 200),
  imp_basin = seq(0, 1.0, length.out = 200)
)
#  Predict canopy_reach
pred_grid = pred_grid %>%
  mutate(
    tcc_reach = predict(m_canopy, newdata = data.frame(imp_reach = imp_reach)),
  )
# Predict bibi
pred_grid = pred_grid %>%
  mutate(
    bibi = predict(m_bibi, newdata = data.frame(tcc_reach = tcc_reach, imp_basin = imp_basin))
  )
# Predict predator richness
pred_grid = pred_grid %>%
  mutate(
    rich_predator = predict(m_rich,
                            newdata = data.frame(
                              bibi = bibi,
                              tcc_reach = tcc_reach,
                              imp_reach = imp_reach
                            ), type = "response")
  )

fig_6C = ggplot(pred_grid, aes(x = imp_reach, y = imp_basin, z = rich_predator)) +
  # geom_raster(aes(fill = rich_predator), interpolate = TRUE) +
  # scale_fill_gradient(low = alpha("red", 0.7), high = alpha("forestgreen", 0.7)) +
  geom_contour(aes(linetype = "Predator richness"), color = "gray80", breaks = seq(floor(min(pred_grid$rich_predator)), ceiling(max(pred_grid$rich_predator)), by = 1)) +
  geom_text_contour(aes(z = rich_predator),
                    breaks = seq(floor(min(pred_grid$rich_predator)), ceiling(max(pred_grid$rich_predator)), by = 1),
                    stroke = 0.0, color = "gray80", check_overlap = FALSE, skip = 0) +
  geom_segment(data = site_restorations,
    aes(x = imp_reach_2023,  y = imp_basin_2023, xend = imp_reach_2100, yend = imp_basin_2100,
        color = needed_restoration),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  # geom_point(data = site_restorations, inherit.aes = FALSE, aes(x = imp_reach_2023, y = imp_basin_2023, color = restore_need)) +
  geom_point(data = site_restorations, inherit.aes = FALSE, aes(x = imp_reach_2100, y = imp_basin_2100, color = needed_restoration, shape = needed_restoration), size = 2, alpha = 0.9) +
  # scale_linetype_manual(name = NULL, values = "solid") +
  # scale_color_manual(values = c("gray40", "magenta3", "lightseagreen")) +
  scale_x_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_color_manual(values = c("firebrick", "firebrick", "gray40", "purple3")) +
  scale_shape_manual(values = c(17, 15, 16, 18)) +
  labs(x = "Reach impervious %", y = "Catchment impervious %", color = "Required restoration", shape = "Required restoration", linetype = "Contour") +
  coord_fixed(ratio = 1); print(fig_6C)

fig_6A = p_restore_reach_imp
fig_6B = p_restore_basin_aligned_imp

fig_6AB = (fig_6A | fig_6B) + plot_annotation(tag_levels = "A")
ggsave(paste0(out_dir, "/fig_6AB.pdf"), fig_6AB, width = 6, height = 6)
ggsave(paste0(out_dir, "/fig_6C.pdf"), fig_6C +
         theme(plot.margin = margin(0, 0, 0, 0, unit = "pt")), width = 6, height = 6)
