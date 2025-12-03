# 6_predictions.R ---------------------------------------------------------------------------------------------
# Use ICLUS impervious surface projections to predict future changes in B-IBI and predator richness across sites

source("src/global.R")

# Load projections
projections_reach     = readRDS("data/cache/3_calculate_vars/imp_projections_550m.rds")
projections_watershed = readRDS("data/cache/3_calculate_vars/imp_projections_5000m.rds")

# TODO: Compare projected 2020 imp to current 2024 imp at each site to find bias offset

# Site-specific % change in ICLUS impervious surface coverage from 2020 was calculated for each decade and
# multiplied with "ground truth" 2023 NLCD measurements to estimate future coverage at each site.
projections_reach %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))
projections_watershed %>% group_by(scenario, year) %>% summarise(sum_sum_proportion = mean(sum_proportion, na.rm = TRUE))

# Load sem for prediction and extract component models
sem = readRDS("data/cache/4_sem/sem_predator.rds")
m_bibi = sem[[1]]
m_rich_predator = sem[[2]]
data = sem[[3]]

# Select a scenario
s = "a2"

# Predict canopy_reach from imp_reach using arcsinsqrt for proportions
ggplot(data, aes(x = imp_reach, y = canopy_reach)) +
  geom_point() + geom_smooth(method = "lm")

data$imp_reach_tr    = asin(sqrt(data$imp_reach))
data$canopy_reach_tr = asin(sqrt(data$canopy_reach))
m_canopy_reach = lm(canopy_reach_tr ~ imp_reach_tr, data)
summary(m_canopy_reach)

data_pred = projections_reach %>% filter(scenario == s) %>%
  mutate(
    imp_reach = sum_proportion,
    imp_reach_tr  = asin(sqrt(imp_reach)),
    canopy_tr_pred = predict(m_canopy_reach, newdata = data.frame(imp_reach_tr = imp_reach_tr)),
    canopy_reach = (sin(canopy_tr_pred))^2          # back-transform to 0–1
  ) %>% select(site_id, scenario, year, imp_reach, canopy_reach)
summary(data_pred$canopy_reach)

ggplot(data_pred, aes(x = year, y = imp_reach, group = site_id)) +
  geom_line(color = "black", alpha = 0.2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", size = 1.2) +
  labs(x = "Year", y = "Predicted reach impervious %")

# Get projected imp_basin
data_pred = data_pred %>% left_join(projections_watershed %>%
                                    rename(imp_basin = sum_proportion) %>%
                                    select(-percent_change_sum_prop),
                                    by = c("scenario", "site_id", "year"))

ggplot(data_pred, aes(x = year, y = imp_basin, group = site_id)) +
  geom_line(color = "black", alpha = 0.2) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "firebrick", size = 1.2) +
  labs(x = "Year", y = "Projected watershed impervious %")

# Predict bibi from predicted canopy_reach and projected imp_basin
data_pred = data_pred %>% mutate(bibi = predict(m_bibi, newdata = data.frame(
                        canopy_reach = canopy_reach,
                        imp_basin = imp_basin
                      ))) # TODO: arcsinsqrt transformation?

ggplot(data_pred, aes(x = year, y = bibi, group = site_id)) +
  geom_line(aes(color = site_id), alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Predicted B-IBI")

# Predict rich_predator from predicted bibi, predicted canopy_reach, and projected imp_reach
data_pred = data_pred %>% mutate(rich_predator = predict(m_rich_predator, newdata = data.frame(
  bibi = bibi,
  canopy_reach = canopy_reach,
  imp_reach = imp_reach
), type = "response"))

ggplot(data_pred, aes(x = year, y = rich_predator, group = site_id)) +
  geom_line(aes(color = site_id), alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "gray10", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Predicted predator richness")

# Calculate deltas from 2023 to 2100
calc_delta = function(d, var, year_start = "2023", year_end = "2100") {
  delta_name = paste0("delta_", var)
  d %>% filter(year %in% c(year_start, year_end)) %>% select(site_id, year, all_of(var)) %>%
    pivot_wider(names_from = year, values_from = all_of(var), names_prefix = paste0(var, "_")) %>%
    mutate(!!delta_name := .[[paste0(var, "_", year_end)]] - .[[paste0(var, "_", year_start)]])
}
calc_delta <- function(d, var, year_start = "2023", year_end = "2100") {
  delta_name <- paste0("delta_", var)
  delta_pct_name <- paste0("deltapcnt_", var)
  
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

delta_vars = c("imp_basin", "imp_reach", "canopy_reach", "bibi", "rich_predator") %>%
  map(~ calc_delta(data_pred, .x)) %>% reduce(full_join, by = "site_id")

summary(delta_vars$delta_imp_basin)
summary(delta_vars$delta_imp_reach)
summary(delta_vars$delta_canopy_reach)
summary(delta_vars$delta_bibi)
summary(delta_vars$delta_rich_predator)
summary(delta_vars$deltapcnt_rich_predator)

ggplot(
  delta_vars %>%
    mutate(site_id = factor(site_id, levels = delta_vars$site_id[order(delta_vars$bibi_2023)])) %>%
    pivot_longer(c(bibi_2023, bibi_2100), names_to = "year", values_to = "bibi"),
  aes(x = site_id, y = bibi, group = site_id)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "#FFCCCC", alpha = 0.1) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "#FFE5CC", alpha = 0.1) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "#FFFFCC", alpha = 0.1) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "#CCFFCC", alpha = 0.1) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "#CCE5FF", alpha = 0.1) +
  geom_line(color = "gray10") +
  geom_point(aes(color = year)) +
  scale_color_manual(values = c("bibi_2023" = "gray10", "bibi_2100" = "firebrick")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot(delta_vars %>%
         mutate(site_id = factor(site_id, levels = delta_vars$site_id[order(delta_vars$rich_predator_2023)])) %>%
         pivot_longer(c(rich_predator_2023, rich_predator_2100), names_to = "year", values_to = "rich_predator"),
  aes(x = site_id, y = rich_predator, group = site_id)) +
  geom_line() +
  geom_point(aes(color = year)) +
  scale_color_manual(values = c("rich_predator_2023" = "gray10", "rich_predator_2100" = "firebrick")) +
  scale_y_continuous(breaks = seq(0, 11, by = 2), limits = c(0, 11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Estimated % change in bird richness per site in 2050 as a function of % change in imp basin (A), imp local (B)
ggplot(delta_vars, aes(x = delta_imp_basin, y = deltapcnt_rich_predator)) +
  geom_point()
ggplot(delta_vars, aes(x = delta_imp_reach, y = deltapcnt_rich_predator)) +
  geom_point()
ggplot(delta_vars, aes(x = delta_imp_reach, y = delta_imp_basin)) +
  geom_point(aes(color = deltapcnt_rich_predator)) +
  scale_color_gradient(low = "red", high = "black")

## Sensitivity analysis richness ~ imp
pred_grid = expand.grid(
  imp_reach = seq(0, 1.0, length.out = 200),
  imp_basin = seq(0, 1.0, length.out = 200)
)
#  Predict canopy_reach
pred_grid = pred_grid %>%
  mutate(
    imp_reach_tr   = asin(sqrt(imp_reach)),
    canopy_tr_pred = predict(m_canopy_reach, newdata = data.frame(imp_reach_tr = imp_reach_tr)),
    canopy_reach   = sin(canopy_tr_pred)^2
  )
# Predict bibi
pred_grid = pred_grid %>%
  mutate(
    bibi = predict(m_bibi, newdata = data.frame(canopy_reach = canopy_reach, imp_basin = imp_basin))
  )
# Predict predator richness
pred_grid = pred_grid %>%
  mutate(
    rich_predator = predict(m_rich_predator, 
                            newdata = data.frame(
                              bibi = bibi,
                              canopy_reach = canopy_reach,
                              imp_reach = imp_reach
                            ), type = "response")
  )

# Richness depends on imperviousness at both the reach and basin scales.
# As imp_basin increases, imp_reach must decrease to maintain the same richness, and vice-versa.
# The slope of contour lines is < -1 (to maintain the same amount of richness, you can have more
# imp_basin than imp_reach), indicating that richness is more sensitive to imp_reach than imp_basin.
# A larger change in imp_basin is needed to compensate for a smaller change in imp_reach.
# Also, the slope gets steeper as impervious levels get smaller (bottom left corner), suggesting that
# this sensitivity is more pronounced at lower levels of imperviousness (relatively pristine sites).
# Therefore, management interventions reducing imp_reach have the largest potential benefit in less disturbed sites.
ggplot(pred_grid, aes(x = imp_reach, y = imp_basin, z = rich_predator)) +
  geom_raster(aes(fill = rich_predator), interpolate = TRUE) +
  scale_fill_gradient(low = alpha("red", 0.75), high = alpha("white", 0.75)) +
  geom_contour(color = "gray20", breaks = seq(floor(min(pred_grid$rich_predator)),
                                             ceiling(max(pred_grid$rich_predator)), by = 1)) +
  metR::geom_text_contour(aes(z = rich_predator), 
                          breaks = seq(floor(min(pred_grid$rich_predator)), ceiling(max(pred_grid$rich_predator)), by = 1),
                          stroke = 0.0, color = "gray20", check_overlap = FALSE, skip = 0) +
  geom_point(data = delta_vars, aes(x = imp_reach_2100, y = imp_basin_2100), inherit.aes = FALSE, color = "firebrick") +
  geom_point(data = delta_vars, aes(x = imp_reach_2023, y = imp_basin_2023), inherit.aes = FALSE, color = "gray10") +
  coord_fixed(ratio = 1)






# Assuming imp_basin values are inevitable in 2100, estimate per site what value of imp_reach
# is needed to maintain the same amount of richness as in 2023

# Filter to the 2023 richness and projected imp_basin in 2100
sites_2100 <- data_pred %>%
  filter(year %in% c("2023", "2100")) %>%
  select(site_id, scenario, year, imp_reach, canopy_reach, imp_basin, bibi, rich_predator) %>%
  pivot_wider(
    names_from = year,
    values_from = c(imp_reach, canopy_reach, imp_basin, bibi, rich_predator),
    names_sep = "_"
  )

# Function that returns richness given imp_reach and site-specific values
richness_for_site <- function(imp_try, imp_basin_2100, m_canopy_reach, m_bibi, m_rich_predator) {
  # predict canopy
  canopy_try <- (sin(predict(m_canopy_reach, newdata = data.frame(
    imp_reach_tr = asin(sqrt(imp_try))
  ))))^2
  
  # predict bibi
  bibi_try <- predict(m_bibi, newdata = data.frame(
    canopy_reach = canopy_try,
    imp_basin = imp_basin_2100
  ))
  
  # predict richness
  predict(m_rich_predator, newdata = data.frame(
    bibi = bibi_try,
    canopy_reach = canopy_try,
    imp_reach = imp_try
  ), type = "response")
}

sites_2100 <- sites_2100 %>% # Find imp reach target values 
  rowwise() %>%
  mutate(
    imp_reach_target = tryCatch({
      uniroot(
        function(x) richness_for_site(x, imp_basin_2100, m_canopy_reach, m_bibi, m_rich_predator) - rich_predator_2023,
        lower = 0, upper = 1
      )$root
    }, error = function(e) NA)
  ) %>%
  ungroup()

sites_2100 <- sites_2100 %>% # Find canopy reach target values
  rowwise() %>%
  mutate(
    canopy_reach_target = (sin(predict(m_canopy_reach, newdata = data.frame(
      imp_reach_tr = asin(sqrt(imp_reach_target))
    ))))^2
  ) %>%
  ungroup()

sites_2100 = sites_2100 %>% arrange(imp_reach_2023) %>%
  mutate(site_id = factor(site_id, levels = site_id))

sites_2100_long = sites_2100 %>%
  pivot_longer(cols = c(imp_reach_2023, imp_reach_target, imp_reach_2100),
               names_to = "type", values_to = "imp_reach") %>%
  mutate(type = factor(type,
                       levels = c("imp_reach_2023", "imp_reach_target", "imp_reach_2100"),
                       labels = c("Current (2023)", "Restoration target", "No action (2100)")))

ggplot(sites_2100_long, aes(y = site_id)) +
  geom_segment(data = sites_2100,
               aes(x = imp_reach_2023, xend = imp_reach_target, y = site_id, yend = site_id),
               color = "gray40") +
  geom_point(aes(x = imp_reach, color = type, shape = type)) +
  scale_color_manual(values = c("Current (2023)" = "black", "Restoration target" = "forestgreen", "No action (2100)" = "grey40")) +
  scale_shape_manual(values = c("Current (2023)" = 19, "Restoration target" = 19, "No action (2100)" = 1)) +
  coord_flip() +
  labs(x = "Reach-scale imperviousness", y = "Site ID", color = "", shape = "",
       title = "Restoration efforts to maintain current predator richness by 2100 assuming unimpeded watershed urbanization") +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))

summary(sites_2100 %>% mutate(imp_reach_target_delta = imp_reach_2023 -imp_reach_target) %>% pull(imp_reach_target_delta))

# TODO: Canopy as well
# ggplot(sites_2100_long, aes(y = site_id)) +
#   geom_segment(data = sites_2100,
#                aes(x = canopy_reach_2023, xend = canopy_reach_target, y = site_id, yend = site_id),
#                color = "gray40") +
#   geom_point(aes(x = canopy_reach_2023)) +
#   geom_point(aes(x = canopy_reach_target), color = "forestgreen") +
#   geom_point(aes(x = canopy_reach_2100), color = "gray40", shape = 1) +
#   coord_flip() +
#   labs(x = "Reach-scale canopy cover", y = "Site ID", color = "", shape = "",
#        title = "Restoration efforts to maintain current predator richness by 2100 assuming unimpeded watershed urbanization") +
#   theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1))
