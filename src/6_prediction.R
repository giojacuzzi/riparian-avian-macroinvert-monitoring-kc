# 6_predictions.R ---------------------------------------------------------------------------------------------
# Use ICLUS impervious surface projections to predict future changes in B-IBI and predator richness across sites

source("src/global.R")

# NOTE: 2023 values in plots are currently model predictions, not ground-truth
# Deltas between 2023 observed response variables (canopy, bibi, richness) and 2100 predictions
# may be positive OR negative because it includes residual error from component model fit
# Canopy increasing reflects the presence of negative residuals in the 2023 data, not that the model actually predicts a gain in cover.
use_2023_observed_as_baseline = FALSE

# Load projections
projections_reach     = readRDS("data/cache/3_calculate_vars/imp_projections_550m.rds")
projections_watershed = readRDS("data/cache/3_calculate_vars/imp_projections_5000m.rds")

# Site-specific % change in ICLUS impervious surface coverage from 2020 was calculated for each decade and
# multiplied with "ground truth" 2023 NLCD measurements to estimate future coverage at each site.
projections_reach %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))
projections_watershed %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))

# Load sem for prediction and extract component models
sem = readRDS("data/cache/4_sem/sem_predator.rds")
m_canopy = sem[[1]]
m_bibi   = sem[[2]]
m_rich   = sem[[3]]
data     = sem[[4]]

# Only retain sites with projections
data = data %>% filter(site_id %in% projections_reach$site_id) %>% arrange(site_id)

# Select a scenario
s = "a2"

# Get projected imp_reach and imp_basin

data_pred = projections_reach %>% filter(scenario == s) %>%
  mutate(imp_reach = sum_proportion)
if (use_2023_observed_as_baseline) {
  data_pred = data_pred %>% filter(year != "2023")
}

p_proj_imp_reach = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = imp_reach, group = site_id)) +
  # geom_line(color = "black", alpha = 0.2) +
  geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0)) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", size = 1.2) +
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

p_proj_imp_basin = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = imp_basin, group = site_id)) +
  # geom_line(color = "black", alpha = 0.2) +
  geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_discrete(expand = c(0,0)) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", size = 1.2) +
  # geom_text(
  #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
  #   aes(label = site_id), hjust = -0.1, size = 2, color = "gray60"
  # ) +
  labs(x = "Year", y = "", title = "Impervious % (watershed)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Predict canopy_reach from imp_reach
data_pred = data_pred %>% mutate(
  canopy_reach = predict(m_canopy, newdata = data.frame(imp_reach = imp_reach))
)

p_proj_canopy_reach = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = plogis(canopy_reach), group = site_id)) +
  # geom_line(color = "black", alpha = 0.2) +
  geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1)) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "forestgreen", size = 1.2) +
  # geom_text(
  #   data = data_pred %>% group_by(site_id) %>% filter(year == max(year)) %>% ungroup(),
  #   aes(label = site_id), hjust = -0.1, size = 2, color = "gray60"
  # ) +
  labs(x = "Year", y = "", title = "Canopy % (reach)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot time series environmental projections
print(p_proj_imp_basin + p_proj_imp_reach + p_proj_canopy_reach)

# Predict bibi from predicted canopy_reach and projected imp_basin
data_pred = data_pred %>% mutate(
  bibi = predict(m_bibi, newdata = data.frame(canopy_reach = canopy_reach, imp_basin = imp_basin))
)

p_pred_bibi = ggplot(data_pred %>% filter(year != "2023"), aes(x = year, y = plogis(bibi) * 100, group = site_id)) +
  # geom_line(color = "black", alpha = 0.2) +
  geom_smooth(aes(group = site_id), method = "loess", se = FALSE, size = 0.5, color = alpha("black", 0.2)) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "royalblue", size = 1.2) +
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
  canopy_reach = canopy_reach,
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
  scale_y_continuous(expand = c(0, 0), limits = c(0, 12), breaks = seq(0,12,2)) + 
  scale_x_discrete(labels = function(x) ifelse(seq_along(x) %% 2 == 0, x, ""), expand = c(0, 0)) + 
  labs(x = "Year", y = "Projected predator richness") +
  theme(plot.margin = margin(0,0,0,0))

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

delta_vars = c("imp_basin", "imp_reach", "canopy_reach", "bibi", "rich_predator") %>%
  map(~ calc_delta(data_pred, .x)) %>% reduce(full_join, by = "site_id")

summary(delta_vars$delta_imp_basin)
summary(delta_vars$delta_imp_reach)
summary(delta_vars$delta_canopy_reach)
summary(delta_vars$delta_bibi)
summary(delta_vars$delta_rich_predator)
summary(delta_vars$deltapcnt_rich_predator)

p_bibi_delta = ggplot(
  delta_vars %>%
    mutate(site_id = factor(site_id, levels = delta_vars$site_id[order(delta_vars$bibi_2023, decreasing = TRUE)])) %>%
    pivot_longer(c(bibi_2023, bibi_2100), names_to = "year", values_to = "bibi"),
  aes(x = site_id, y = plogis(bibi) * 100, group = site_id)) +
  # geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "#FFCCCC", alpha = 0.1) +
  # geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "#FFE5CC", alpha = 0.1) +
  # geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "#FFFFCC", alpha = 0.1) +
  # geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "#CCFFCC", alpha = 0.1) +
  # geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "#CCE5FF", alpha = 0.1) +
  # geom_hline(yintercept = 80, color = "#CCFFCC", alpha = 1) +
  # geom_hline(yintercept = 60, color = "#FFFFCC", alpha = 1) +
  # geom_hline(yintercept = 40, color = "#FFE5CC", alpha = 1) +
  # geom_hline(yintercept = 20, color = "#FFCCCC", alpha = 1) +
  geom_line(color = "gray10") +
  geom_point(aes(color = year)) +
  scale_color_manual(values = c("bibi_2023" = "gray10", "bibi_2100" = "royalblue")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100), limits = c(0, 100), expand = c(0,0)) +
  labs(x = "", y = "", color = "") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
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
  scale_y_continuous(breaks = seq(0, 14, by = 1), limits = c(0, 12), expand = c(0,0)) +
  labs(x = "Site", y = "", color = "") +
  theme(axis.text.x = element_text(size = 8, angle = 90, hjust = 1),
        plot.margin = margin(0,0,0,0),
        legend.position = "bottom"
        )

# Plot time series and delta projections for responses
p_proj_responses = (p_pred_bibi + p_bibi_delta) / (p_pred_rich_predator + p_rich_predator_delta)
print(p_proj_responses + plot_annotation(tag_levels = "A"))

# Watershed-only restoration (e.g. integrated watershed management plans) -----------------------------
# Assuming imp_reach and canopy_reach values cannot be changed (e.g. local-scale restoration is
# not feasible due to private land ownership), estimate per site what value of imp_basin
# is needed to maintain the same amount of richness as in 2023

# Filter to the 2023 richness and projected imp_basin in 2100
sites_2100 = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, scenario, year, imp_reach, canopy_reach, imp_basin, bibi, rich_predator) %>%
  pivot_wider(names_from = year, values_from = c(imp_reach, canopy_reach, imp_basin, bibi, rich_predator), names_sep = "_")

# Function that returns richness given imp_reach and site-specific values
richness_for_site_basin = function(imp_basin_try, imp_reach_2100, m_canopy, m_bibi, m_rich) {
  # predict canopy
  canopy_try = predict(m_canopy, newdata = data.frame(
    imp_reach = imp_reach_2100
  ))
  # predict bibi
  bibi_try = predict(m_bibi, newdata = data.frame(
    canopy_reach = canopy_try,
    imp_basin = imp_basin_try
  ))
  # predict richness
  predict(m_rich, newdata = data.frame(
    bibi = bibi_try,
    canopy_reach = canopy_try,
    imp_reach = imp_reach_2100
  ), type = "response")
}

# Solve for required imp_basin target
sites_2100 = sites_2100 %>% rowwise() %>%
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
sites_2100 = sites_2100 %>% mutate(imp_basin_target_delta = imp_basin_2023 - imp_basin_target)

# Summary statistics
summary(sites_2100 %>% pull(imp_basin_target_delta))

# NOTE: NA's suggest that no amount of basin level restoration can maintain richness levels
# in the face of how degraded the stream reach becomes by 2100
sites_2100 = sites_2100 %>% mutate(restoration_possible = ifelse(is.na(imp_basin_target), "No", "Yes"))
site_id_no_watershed_restoration_possible = sites_2100 %>% filter(restoration_possible == "No") %>% pull(site_id)

# Plot
sites_2100 = sites_2100 %>% arrange(imp_basin_2023) %>% mutate(site_id = factor(site_id, levels = site_id))
sites_2100_basin_long = sites_2100 %>%
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
  geom_segment(data = sites_2100, aes(x = imp_basin_2023, xend = imp_basin_2100, y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100, aes(x = imp_basin_2023, xend = imp_basin_target, y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = imp_basin, color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Watershed impervious %", y = "", color = "", shape = "") +
  theme(legend.position = "none")

# Reach-only restoration (e.g. public lands, private land easements, land trusts) ---------------------
# Assuming imp_basin values are inevitable in 2100, estimate per site what value of imp_reach
# is needed to maintain the same amount of richness as in 2023

# Filter to the 2023 richness and projected imp_basin in 2100
sites_2100 = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, scenario, year, imp_reach, canopy_reach, imp_basin, bibi, rich_predator) %>%
  pivot_wider(names_from = year, values_from = c(imp_reach, canopy_reach, imp_basin, bibi, rich_predator), names_sep = "_")

# Function that returns richness given imp_reach and site-specific values
richness_for_site_reach = function(imp_reach_try, imp_basin_2100, m_canopy, m_bibi, m_rich) {
  # predict canopy
  canopy_try = predict(m_canopy, newdata = data.frame(
    imp_reach = imp_reach_try
  ))
  # predict bibi
  bibi_try = predict(m_bibi, newdata = data.frame(
    canopy_reach = canopy_try,
    imp_basin = imp_basin_2100
  ))
  # predict richness
  predict(m_rich, newdata = data.frame(
    bibi = bibi_try,
    canopy_reach = canopy_try,
    imp_reach = imp_reach_try
  ), type = "response")
}

# Solve for required imp_reach target
sites_2100 = sites_2100 %>% rowwise() %>%
  mutate(
    imp_reach_target = tryCatch({ # TODO: Some sites returning NA for imp_reach_target
      # Is this because, given unimpeded basin urbanization, no possible reach-scale imperviousness reduction can restore 2023 richness levels?
      uniroot(
        function(x) richness_for_site_reach(x, imp_basin_2100, m_canopy, m_bibi, m_rich) - rich_predator_2023,
        lower = 0, upper = 1
      )$root
    }, error = function(e) NA)
  ) %>% ungroup()

sites_2100 = sites_2100 %>% mutate(imp_reach_target_delta = imp_reach_2023 - imp_reach_target)
sites_2100 = sites_2100 %>% mutate(restoration_possible = ifelse(is.na(imp_reach_target), "No", "Yes"))
site_id_no_reach_restoration_possible = sites_2100 %>% filter(restoration_possible == "No") %>% pull(site_id)
sites_2100 = sites_2100 %>% arrange(imp_reach_2023) %>%
  mutate(site_id = factor(site_id, levels = site_id))

sites_2100_long = sites_2100 %>%
  pivot_longer(cols = c(imp_reach_2023, imp_reach_target, imp_reach_2100), names_to = "type", values_to = "imp_reach") %>%
  mutate(type = factor(type,
                       levels = c("imp_reach_2023", "imp_reach_target", "imp_reach_2100"),
                       labels = c("Current (2023)", "Restoration to maintain richness", "No action (2100)")))

p_restore_reach_imp = ggplot(sites_2100_long, aes(y = site_id)) +
  geom_segment(data = sites_2100, aes(x = imp_reach_2023, xend = imp_reach_2100, y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100, aes(x = imp_reach_2023, xend = imp_reach_target, y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = imp_reach, color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Reach impervious %", y = "Site", color = "", shape = "") +
  theme(legend.position = "bottom")

summary(sites_2100 %>% pull(imp_reach_target_delta))

# Find canopy reach target values
sites_2100 = sites_2100 %>% rowwise() %>%
  mutate(
    canopy_reach_target = predict(m_canopy, newdata = data.frame(
      imp_reach = imp_reach_target
    ))
  ) %>% ungroup()
sites_2100 = sites_2100 %>% mutate(canopy_reach_target_delta = canopy_reach_2023 - canopy_reach_target)
sites_2100 = sites_2100 %>% arrange(imp_reach_2023) %>%
  mutate(site_id = factor(site_id, levels = site_id))

sites_2100_long = sites_2100 %>%
  pivot_longer(cols = c(canopy_reach_2023, canopy_reach_target, canopy_reach_2100), names_to = "type", values_to = "canopy_reach") %>%
  mutate(type = factor(type,
                       levels = c("canopy_reach_2023", "canopy_reach_target", "canopy_reach_2100"),
                       labels = c("Current (2023)", "Restoration to maintain richness", "No action (2100)")))

p_restore_reach_canopy = ggplot(sites_2100_long, aes(y = site_id)) +
  geom_segment(data = sites_2100, aes(x = plogis(canopy_reach_2023), xend = plogis(canopy_reach_2100), y = site_id, yend = site_id, color = restoration_possible)) +
  scale_color_manual(values = c("Yes" = "grey40", "No" = "firebrick")) +
  new_scale_color() +
  geom_segment(data = sites_2100, aes(x = plogis(canopy_reach_2023), xend = plogis(canopy_reach_target), y = site_id, yend = site_id), color = "forestgreen", linetype = "dotted") +
  geom_point(aes(x = plogis(canopy_reach), color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration to maintain richness" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration to maintain richness" = 19, "No action (2100)" = 19)) +
  labs(x = "Reach canopy %", y = "Site", color = "", shape = "") +
  theme(legend.position = "bottom")

# Plot both
#
# Interesting points / findings:
# - Reach-scale restoration includes illustrated reductions in impervious % as well as commensurate increases in canopy %.
# - The magnitude of restoration required is not necessarily contingent on the magnitude of future impairment.
# - For many sites, even dramatic restoration efforts at the watershed scale cannot
# mitigate the overwhelming impacts of local-scale degradation within the reach.
p_restore_imp = p_restore_reach_imp + p_restore_basin_imp
print(p_restore_imp + plot_annotation(tag_levels = "A"))

p_total = (p_proj_responses | p_restore_imp) + plot_annotation(tag_levels = "A")
print(p_total)
ggsave("p_total.pdf", p_total, width = 13, height = 7)

# Both reach and watershed restoration -----------------------------------------------------------------

site_restorations = data.frame(
  site_id = data$site_id,
  restore_need_watershed = data$site_id %in% site_id_no_reach_restoration_possible,
  restore_need_reach     = data$site_id %in% site_id_no_watershed_restoration_possible
)
site_restorations = site_restorations %>% filter(site_id %in% delta_vars$site_id)
site_restorations$restorable = !site_restorations$restore_need_watershed & !site_restorations$restore_need_reach
site_restorations = site_restorations %>% left_join(delta_vars %>% select(site_id, imp_basin_2023, imp_basin_2100, imp_reach_2023, imp_reach_2100), by = "site_id")

site_restorations$restore_need = with(site_restorations,
                                     ifelse(restore_need_watershed, "restore_need_watershed",
                                            ifelse(restore_need_reach, "restore_need_reach",
                                                   "restorable single scale"))
)

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
    canopy_reach = predict(m_canopy, newdata = data.frame(imp_reach = imp_reach)),
  )
# Predict bibi
pred_grid = pred_grid %>%
  mutate(
    bibi = predict(m_bibi, newdata = data.frame(canopy_reach = canopy_reach, imp_basin = imp_basin))
  )
# Predict predator richness
pred_grid = pred_grid %>%
  mutate(
    rich_predator = predict(m_rich,
                            newdata = data.frame(
                              bibi = bibi,
                              canopy_reach = canopy_reach,
                              imp_reach = imp_reach
                            ), type = "response")
  )

# Richness responds to imperviousness at both the reach and basin scales.
# Tradeoff: As imp_basin increases, imp_reach must decrease to maintain the same richness, and vice-versa.
# Responses are not symmetric between scales: The slope of contour lines is < -1 (to maintain the same amount of richness, you can have more
# imp_basin than imp_reach), indicating that richness is more sensitive to imp_reach than imp_basin.
# A larger change in imp_basin is needed to compensate for a smaller change in imp_reach.
# Also, the slope gets steeper as impervious levels get smaller (bottom left corner), suggesting that
# this sensitivity is more pronounced at lower levels of imperviousness (relatively pristine sites).
# Therefore, management interventions reducing imp_reach have the largest potential benefit in less disturbed sites.
# The distance between contours increases as imperviousness increases across both scales --> 1) declines in richness slow as a stream becomes
# heavily degraded, and 2) greater amounts of restoration are needed to conserve richness at more imperiled sites (conversely, pristine or
# moderately impacted streams benefit from relatively small interventions).
# Streams that require multi-scale restoration are 1) projected to experience dramatic changes in
# imperviousness along at one or both scales (their lines are very long), and/or 2) already have relatively
# low levels of imperviousness at one or both scales, and so are limited in their "wiggle room / floor" for restoration adjustment along
# one or both scales (axes)
p_imp_rich = ggplot(pred_grid, aes(x = imp_reach, y = imp_basin, z = rich_predator)) +
  # geom_raster(aes(fill = rich_predator), interpolate = TRUE) +
  # scale_fill_gradient(low = alpha("red", 0.7), high = alpha("forestgreen", 0.7)) +
  geom_contour(aes(linetype = "Predator richness"), color = "gray70", breaks = seq(floor(min(pred_grid$rich_predator)), ceiling(max(pred_grid$rich_predator)), by = 1)) +
  geom_text_contour(aes(z = rich_predator),
                    breaks = seq(floor(min(pred_grid$rich_predator)), ceiling(max(pred_grid$rich_predator)), by = 1),
                    stroke = 0.0, color = "gray70", check_overlap = FALSE, skip = 0) +
  geom_segment(data = site_restorations,
    aes(x = imp_reach_2023,  y = imp_basin_2023, xend = imp_reach_2100, yend = imp_basin_2100, color = restore_need),
    inherit.aes = FALSE, alpha = 0.5
  ) +
  # geom_point(data = site_restorations, inherit.aes = FALSE, aes(x = imp_reach_2023, y = imp_basin_2023, color = restore_need)) +
  geom_point(data = site_restorations, inherit.aes = FALSE, aes(x = imp_reach_2100, y = imp_basin_2100, color = restore_need), size = 2, alpha = 0.9) +
  scale_linetype_manual(name = NULL, values = "solid") +
  scale_color_manual(values = c("gray40", "magenta3", "lightseagreen")) +
  labs(x = "Reach impervious %", y = "Watershed impervious %", color = "Restoration scale needed") +
  coord_fixed(ratio = 1); print(p_imp_rich)

