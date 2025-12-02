# 6_predictions.R ---------------------------------------------------------------------------------------------
# Use ICLUS impervious surface projections to predict future changes in B-IBI and predator richness across sites

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

# Predict bibi from predicted canopy_reach and projected imp_basin
data_pred = data_pred %>% left_join(projections_watershed %>%
                                    rename(imp_basin = sum_proportion) %>%
                                    select(-percent_change_sum_prop),
                                    by = c("scenario", "site_id", "year"))

data_pred = data_pred %>% mutate(bibi = predict(m_bibi, newdata = data.frame(
                        canopy_reach = canopy_reach,
                        imp_basin = imp_basin
                      ))) # TODO: arcsinsqrt transformation?

# Predicted change in bibi across sites
ggplot(data_pred, aes(x = year, y = bibi, group = site_id)) +
  geom_line(aes(color = site_id), alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Predicted B-IBI")

delta_bibi_2100 = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, year, bibi) %>%
  pivot_wider(names_from = year, values_from = bibi, names_prefix = "bibi_") %>%
  mutate(delta_bibi = bibi_2100 - bibi_2023)
summary(delta_bibi_2100$delta_bibi)

delta_bibi_2100_long = delta_bibi_2100 %>%
  mutate(site_id = factor(site_id, levels = delta_bibi_2100$site_id[order(delta_bibi_2100$delta_bibi)])) %>%
  pivot_longer(c(bibi_2023, bibi_2100), names_to = "year", values_to = "bibi")

ggplot(delta_bibi_2100_long, aes(x = site_id, y = bibi, group = site_id)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "#FFCCCC", alpha = 0.1) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "#FFE5CC", alpha = 0.1) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "#FFFFCC", alpha = 0.1) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "#CCFFCC", alpha = 0.1) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "#CCE5FF", alpha = 0.1) +
  geom_line() +
  geom_point(aes(color = year)) +
  scale_color_manual(values = c("bibi_2023" = "black", "bibi_2100" = "red")) +
  scale_y_continuous(breaks = c(0, 20, 40, 60, 80, 100))

# Predict rich_predator from predicted bibi, predicted canopy_reach, and projected imp_reach
data_pred = data_pred %>% mutate(rich_predator = predict(m_rich_predator, newdata = data.frame(
  bibi = bibi,
  canopy_reach = canopy_reach,
  imp_reach = imp_reach
), type = "response"))

ggplot(data_pred, aes(x = year, y = rich_predator, group = site_id)) +
  geom_line(aes(color = site_id), alpha = 0.5) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", color = "black", size = 1.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Year", y = "Predicted predator richness")

delta_rich_predator_2100 = data_pred %>% filter(year %in% c("2023", "2100")) %>%
  select(site_id, year, rich_predator) %>%
  pivot_wider(names_from = year, values_from = rich_predator, names_prefix = "rich_predator_") %>%
  mutate(delta_rich_predator = rich_predator_2100 - rich_predator_2023)
summary(delta_rich_predator_2100$delta_rich_predator)

delta_rich_predator_2100_long = delta_rich_predator_2100 %>%
  mutate(site_id = factor(site_id, levels = delta_rich_predator_2100$site_id[order(delta_rich_predator_2100$delta_rich_predator)])) %>%
  pivot_longer(c(rich_predator_2023, rich_predator_2100), names_to = "year", values_to = "rich_predator")

ggplot(delta_rich_predator_2100_long, aes(x = site_id, y = rich_predator, group = site_id)) +
  geom_line() +
  geom_point(aes(color = year)) +
  scale_color_manual(values = c("rich_predator_2023" = "black", "rich_predator_2100" = "red")) +
  scale_y_continuous(breaks = seq(0, 11, by = 2), limits = c(0, 11))

# TODO: What amount of impervious mitigation is necessary to constrain predator losses to < 1 species?

# Figure
# Heatmap of range in local imperviousness against range in basin imperviousness in 2050
# z color is richness change?
# dots for actual site locations


# Figure:
# Estimated % change in bird richness per site in 2050 as a function of % change in imp basin (A), imp local (B), canopy local (C)
#
# x variable change % (relative to original)
# y variable change % (relative to original)

