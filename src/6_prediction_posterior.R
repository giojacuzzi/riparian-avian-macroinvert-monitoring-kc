source("src/global.R")

# Select a scenario (a1, a2, b1, b2)
s = "a2"

# Load projections
path_projections_reach = "data/cache/3_calculate_vars/imp_projections_550m.rds"
path_projections_watershed = "data/cache/3_calculate_vars/imp_projections_5000m.rds"
projections_reach     = readRDS(path_projections_reach)
projections_watershed = readRDS(path_projections_watershed)

# Get projected imp_reach and imp_basin

data_project = projections_reach %>% filter(scenario == s) %>%
  mutate(imp_reach = sum_proportion)
data_project = data_project %>% left_join(projections_watershed %>%
                                      rename(imp_basin = sum_proportion) %>%
                                      select(-percent_change_sum_prop),
                                    by = c("scenario", "site_id", "year"))

sem_diet_draws = readRDS("data/cache/4_sem/sem_diet_draws.rds")
draws = length(sem_diet_draws)

summary_imp_reach = vector("list", draws)
summary_imp_basin = vector("list", draws)
summary_tcc_reach = vector("list", draws)
summary_bibi = vector("list", draws)
summary_rich_predator = vector("list", draws)

message("Calculating projection stats across all ", draws, " draws")
pb = progress_bar$new(format = "[:bar] :percent :elapsedfull (ETA :eta)", total = draws, clear = FALSE)
for (draw in 1:draws) {
  sem = sem_diet_draws[[draw]] # readRDS("data/cache/4_sem/sem_diet.rds")
  m_canopy = sem[[1]]
  m_bibi   = sem[[2]]
  m_rich   = sem[[3]]
  data     = sem[[4]]
  
  # Only retain sites with projections
  data = data %>% filter(site_id %in% projections_reach$site_id) %>% arrange(site_id)
  
  data_pred = data_project
  
  # imp_reach and imp_basin changes -------------------------------------------------------------------
  stats_imp_reach = data_pred %>% filter(year %in% c("2023", "2100")) %>% select(site_id, year, imp_reach) %>%
    pivot_wider(names_from = year, values_from = imp_reach, names_prefix = "imp_reach_") %>%
    mutate(
      imp_reach_delta_2100  = imp_reach_2100 - imp_reach_2023,
      imp_reach_factor_2100 = imp_reach_2100 / imp_reach_2023
    )
  stats_imp_reach[is.nan(stats_imp_reach$imp_reach_factor_2100), 'imp_reach_factor_2100'] = 1
  summary(stats_imp_reach$imp_reach_delta_2100)
  summary(stats_imp_reach$imp_reach_factor_2100)
  
  stats_imp_basin = data_pred %>% filter(year %in% c("2023", "2100")) %>% select(site_id, year, imp_basin) %>%
    pivot_wider(names_from = year, values_from = imp_basin, names_prefix = "imp_basin_") %>%
    mutate(
      imp_basin_delta_2100  = imp_basin_2100 - imp_basin_2023,
      imp_basin_factor_2100 = imp_basin_2100 / imp_basin_2023
    )
  stats_imp_basin[is.nan(stats_imp_basin$imp_basin_factor_2100), 'imp_basin_factor_2100'] = 1
  summary(stats_imp_basin$imp_basin_delta_2100)
  summary(stats_imp_basin$imp_basin_factor_2100)
  
  # Predict canopy_reach from imp_reach -------------------------------------------------------------------
  data_pred = data_pred %>% mutate(
    tcc_reach = predict(m_canopy, newdata = data.frame(imp_reach = imp_reach))
  )
  
  stats_tcc_reach = data_pred %>% filter(year %in% c("2023", "2100")) %>% select(site_id, year, tcc_reach) %>%
    mutate(tcc_reach = plogis(tcc_reach)) %>%
    pivot_wider(names_from = year, values_from = tcc_reach, names_prefix = "tcc_reach_") %>%
    mutate(
      tcc_reach_delta_2100  = tcc_reach_2100 - tcc_reach_2023,
      tcc_reach_factor_2100 = tcc_reach_2100 / tcc_reach_2023
    )
  summary(stats_tcc_reach$tcc_reach_delta_2100)
  summary(stats_tcc_reach$tcc_reach_factor_2100)
  
  # Predict bibi from predicted canopy_reach and projected imp_basin --------------------------------------
  data_pred = data_pred %>% mutate(
    bibi = predict(m_bibi, newdata = data.frame(tcc_reach = tcc_reach, imp_basin = imp_basin))
  )
  
  stats_bibi = data_pred %>% filter(year %in% c("2023", "2100")) %>% select(site_id, year, bibi) %>%
    mutate(bibi = plogis(bibi)) %>%
    pivot_wider(names_from = year, values_from = bibi, names_prefix = "bibi_") %>%
    mutate(
      bibi_delta_2100  = bibi_2100 - bibi_2023,
      bibi_factor_2100 = bibi_2100 / bibi_2023
    )
  summary(stats_bibi$bibi_delta_2100)
  summary(stats_bibi$bibi_factor_2100)
  
  # Predict rich_predator from predicted bibi, predicted canopy_reach, and projected imp_reach ------------
  data_pred = data_pred %>% mutate(
    rich_predator = predict(m_rich, newdata = data.frame(
      bibi = bibi,
      tcc_reach = tcc_reach,
      imp_reach = imp_reach
    )))
  
  stats_rich_predator = data_pred %>% filter(year %in% c("2023", "2100")) %>% select(site_id, year, rich_predator) %>%
    pivot_wider(names_from = year, values_from = rich_predator, names_prefix = "rich_predator_") %>%
    mutate(
      rich_predator_delta_2100  = rich_predator_2100 - rich_predator_2023,
      rich_predator_factor_2100 = rich_predator_2100 / rich_predator_2023
    )
  summary(stats_rich_predator$rich_predator_delta_2100)
  summary(stats_rich_predator$rich_predator_factor_2100)
  
  # Store results for the draw -----------------
  summary_imp_reach[[draw]] = stats_imp_reach %>% select(imp_reach_delta_2100, imp_reach_factor_2100) %>%
    summarise(across(everything(), ~ list(summary(.x)))) %>%
    pivot_longer(everything(), names_to = "stat") %>% unnest_wider(value) %>% mutate(draw = draw)
  
  summary_imp_basin[[draw]] = stats_imp_basin %>% select(imp_basin_delta_2100, imp_basin_factor_2100) %>%
    summarise(across(everything(), ~ list(summary(.x)))) %>%
    pivot_longer(everything(), names_to = "stat") %>% unnest_wider(value) %>% mutate(draw = draw)
  
  summary_tcc_reach[[draw]] = stats_tcc_reach %>% select(tcc_reach_delta_2100, tcc_reach_factor_2100) %>%
    summarise(across(everything(), ~ list(summary(.x)))) %>%
    pivot_longer(everything(), names_to = "stat") %>% unnest_wider(value) %>% mutate(draw = draw)
  
  summary_bibi[[draw]] = stats_bibi %>% select(bibi_delta_2100, bibi_factor_2100) %>%
    summarise(across(everything(), ~ list(summary(.x)))) %>%
    pivot_longer(everything(), names_to = "stat") %>% unnest_wider(value) %>% mutate(draw = draw)
  
  summary_rich_predator[[draw]] = stats_rich_predator %>% select(rich_predator_delta_2100, rich_predator_factor_2100) %>%
    summarise(across(everything(), ~ list(summary(.x)))) %>%
    pivot_longer(everything(), names_to = "stat") %>% unnest_wider(value) %>% mutate(draw = draw)
  
  pb$tick()
}

summary_imp_reach = bind_rows(summary_imp_reach)
summary_imp_basin = bind_rows(summary_imp_basin)
summary_tcc_reach = bind_rows(summary_tcc_reach)
summary_bibi = bind_rows(summary_bibi)
summary_rich_predator = bind_rows(summary_rich_predator)

# Calculate 95% BCIs and save along with scenario column

# Relative to baseline conditions, both reach- and catchment-level ISC are projected to increase on average by roughly an additional 10% of surface area cover by the year 2100 (a mean increase by a factor of approximately 1.6 across all streams).
summary_imp_reach %>% clean_names() %>% filter(stat == "imp_reach_delta_2100") %>%
  summarise(posterior_mean = mean(.data$mean), bci_lower = quantile(.data$mean, 0.025), bci_upper = quantile(.data$mean, 0.975))
summary_imp_basin %>% clean_names() %>% filter(stat == "imp_basin_delta_2100") %>%
  summarise(posterior_mean = mean(.data$mean), bci_lower = quantile(.data$mean, 0.025), bci_upper = quantile(.data$mean, 0.975))
summary_imp_reach %>% clean_names() %>% filter(stat == "imp_reach_factor_2100") %>%
  summarise(posterior_mean = mean(.data$mean), bci_lower = quantile(.data$mean, 0.025), bci_upper = quantile(.data$mean, 0.975))
summary_imp_basin %>% clean_names() %>% filter(stat == "imp_basin_factor_2100") %>%
  summarise(posterior_mean = mean(.data$mean), bci_lower = quantile(.data$mean, 0.025), bci_upper = quantile(.data$mean, 0.975))

# At the most dramatically altered stream reach-level ISC is expected to increase fourfold while catchment-level ISC will double.
summary_imp_reach %>% clean_names() %>% filter(stat == "imp_reach_factor_2100") %>%
  summarise(posterior_mean = mean(.data$max), bci_lower = quantile(.data$max, 0.025), bci_upper = quantile(.data$max, 0.975))
summary_imp_basin %>% clean_names() %>% filter(stat == "imp_basin_factor_2100") %>%
  summarise(posterior_mean = mean(.data$max), bci_lower = quantile(.data$max, 0.025), bci_upper = quantile(.data$max, 0.975))

# Consequently, TCC within the reach was predicted to decrease by as much as 22%, or roughly half of baseline canopy cover
summary_tcc_reach %>% clean_names() %>% filter(stat == "tcc_reach_delta_2100") %>%
  summarise(posterior_mean = mean(.data$min), bci_lower = quantile(.data$min, 0.025), bci_upper = quantile(.data$min, 0.975))
summary_tcc_reach %>% clean_names() %>% filter(stat == "tcc_reach_factor_2100") %>%
  summarise(posterior_mean = mean(.data$min), bci_lower = quantile(.data$min, 0.025), bci_upper = quantile(.data$min, 0.975))

# These environmental changes induce shifts in B-IBI (Figure 5A-B), with all but the most pristine rural forested streams degrading by as much as 32 points by year 2100. 
summary_bibi %>% clean_names() %>% filter(stat == "bibi_delta_2100") %>%
  summarise(posterior_mean = mean(.data$min), bci_lower = quantile(.data$min, 0.025), bci_upper = quantile(.data$min, 0.975))

# Collectively, these reductions in habitat and prey are predicted to reduce predator richness by as many as 2 species (–2.36, –1.75 95% BCI)
summary_rich_predator %>% clean_names() %>% filter(stat == "rich_predator_delta_2100") %>%
  summarise(posterior_mean = mean(.data$min), bci_lower = quantile(.data$min, 0.025), bci_upper = quantile(.data$min, 0.975))

