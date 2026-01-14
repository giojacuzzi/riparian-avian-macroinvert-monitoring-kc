# 5_msom_results.R ===================================================================
# Visualize multi-species occupancy modeling results
#
# Input
msom_paths = list(
  "all"     = "data/cache/models/msom_all.rds",
  "diet"    = "data/cache/models/msom_diet.rds",
  "forage"  = "data/cache/models/msom_forage.rds",
  "migrant" = "data/cache/models/msom_migrant.rds"
)

source("src/global.R")

# Load data for multi-species occupancy models --------------------------------------------------

model_data = list()
for (m in names(msom_paths)) {
  message("Loading data for model '", m, "'")
  model_data[[m]] = readRDS(msom_paths[[m]])
}

# msom_summary = model_data$msom_summary
# msom = model_data$msom
# sites = model_data$sites
# species = model_data$species

# Compare occupancy parameters --------------------------------------------------

message("Summarizing occupancy parameters")

groups           = list()
param_alpha_data = list()
occ_coef_summary = list()
occ_effect_sizes = list()
for (m in names(msom_paths)) {
  groups[[m]]           = model_data[[m]]$groups %>% arrange(common_name)
  param_alpha_data[[m]] = model_data[[m]]$param_alpha_data
  
  occ_coef_summary[[m]] = model_data[[m]]$msom_summary %>%
    filter(Reduce(`|`, lapply(
      c("mu.alpha", "sigma.alpha"),
      \(p) str_starts(param, p)
    ))) %>%
    mutate(
      group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
      param = str_remove(param, "\\[\\d+\\]")
    ) %>%
    select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
    left_join(groups[[m]] %>% distinct(group_idx, group), by = "group_idx")
  
  occ_effect_sizes[[m]] = full_join(occ_coef_summary[[m]] %>% filter(str_starts(param, "mu")),
                               param_alpha_data[[m]] %>% mutate(param = paste0("mu.", param), model = m), by='param')
}

occ_effect_sizes = bind_rows(occ_effect_sizes)
occ_effect_sizes$group = str_to_title(addNA(occ_effect_sizes$group))
occ_effect_sizes$model = str_to_title(occ_effect_sizes$model)

ggplot(occ_effect_sizes %>% filter(name == "bibi"),
       aes(x = mean, y = fct_rev(model), group = fct_rev(group), color = as.factor(overlap0) == 0)) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = c("gray20", "orange")) +
  labs(x = "B-IBI coefficient estimate", y = "Functional group") +
  theme(legend.position = "none")

plots_occ = list()
for (param_name in unique(occ_effect_sizes$name)) {
  plots_occ[[param_name]] = ggplot(occ_effect_sizes %>% filter(name == param_name, group != "other"),
             aes(x = mean, y = as.factor(model), group = as.factor(group), shape = as.factor(group), color = as.factor(overlap0) == 0)) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("gray20", "orange")) +
    labs(subtitle = param_name, x = "", y = "")
}

(plots_occ[["bibi"]] + labs(y = "Group") + theme(legend.position = "none")) +
(plots_occ[["canopy_reach"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
(plots_occ[["canopy_basin"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
(plots_occ[["imp_reach"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
(plots_occ[["imp_basin"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()))

# Compare detection parameters --------------------------------------------------

message("Summarizing detection parameters")

groups           = list()
param_beta_data  = list()
det_coef_summary = list()
det_effect_sizes = list()
for (m in names(msom_paths)) {
  groups[[m]]           = model_data[[m]]$groups %>% arrange(common_name)
  param_beta_data[[m]] = model_data[[m]]$param_beta_data
  
  det_coef_summary[[m]] = model_data[[m]]$msom_summary %>%
    filter(Reduce(`|`, lapply(
      c("mu.beta", "sigma.beta"),
      \(p) str_starts(param, p)
    ))) %>%
    mutate(
      group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
      param = str_remove(param, "\\[\\d+\\]")
    ) %>%
    select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
    left_join(groups[[m]] %>% distinct(group_idx, group), by = "group_idx")
  
  det_effect_sizes[[m]] = full_join(det_coef_summary[[m]] %>% filter(str_starts(param, "mu")),
                                    param_beta_data[[m]] %>% mutate(param = paste0("mu.", param), model = m), by='param')
}

det_effect_sizes = bind_rows(det_effect_sizes)
det_effect_sizes$group = addNA(det_effect_sizes$group)

ggplot(det_effect_sizes %>% filter(name == "yday"),
       aes(x = mean, y = as.factor(model), color = as.factor(group))) +
  geom_vline(xintercept = 0, color = "gray") +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
  labs(title = "Community level effect sizes for detection covariates",
       subtitle = m, x = "Coefficient estimate", y = "Parameter")

plots_det = list()
for (param_name in unique(det_effect_sizes$name)) {
  plots_det[[param_name]] = ggplot(det_effect_sizes %>% filter(name == param_name, group != "other"),
                               aes(x = mean, y = as.factor(model), group = as.factor(group), shape = as.factor(group), color = as.factor(overlap0) == 0)) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1, position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("gray20", "orange")) +
    labs(subtitle = param_name, x = "", y = "")
}

(plots_det[["yday"]] + labs(y = "Parameter") + theme(legend.position = "none")) +
  (plots_det[["canopy_reach"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank())) +
  (plots_det[["imp_reach"]] + theme(legend.position = "none", axis.text.y = element_blank(), axis.ticks.y = element_blank()))

# Marginal responses --------------------------------------------------------------------

plots_marg = list()
for (m in names(msom_paths)) {
  
  message("Calculating marginal responses for model '", m, "'")

  param_name = "bibi"
  param_data = param_alpha_data[[m]] %>% filter(name == param_name)
  species = model_data[[m]]$species
  groups  = model_data[[m]]$groups %>% arrange(common_name)
  
  # Mean marginal probabilities of occurrence for the metacommunity in relation to alpha
  coefs = model_data[[m]]$msom_summary %>% filter(stringr::str_starts(param, param_data$param)) %>% arrange(mean) %>%
    mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(common_name = species[species_idx])
  param_scaled_data = param_data %>% pull(scaled)
  p_mu = attr(param_scaled_data[[1]], "scaled:center") # to transform between z-scale and pred_range_original scale
  p_sd = attr(param_scaled_data[[1]], "scaled:scale")
  bound_low  = min(param_scaled_data[[1]]) * p_sd + p_mu
  bound_high = max(param_scaled_data[[1]]) * p_sd + p_mu
  pred_range_original = seq(bound_low, bound_high, by = 0.01) # range of possible alpha values
  pred_range_standardized = (pred_range_original - p_mu) / p_sd
  # species-specific occurrence intercepts u[i]
  intercepts = model_data[[m]]$msom_summary %>% filter(str_starts(param, "u")) %>%
    mutate(species_idx = as.integer(str_extract(param, "\\d+")), common_name = species[species_idx]) %>%
    select(common_name, u_i = mean)
  # species-specific coefficients
  intercepts_and_coeffs = coefs %>% rename(alpha6_i = mean) %>% select(common_name, alpha6_i) %>%
    left_join(intercepts, by = "common_name")
  # predict species-specific occurrence probabilities
  predictions = intercepts_and_coeffs %>%
    rowwise() %>% do({
      i <- .
      tibble(
        common_name = i$common_name,
        idx = pred_range_original,
        psi = plogis(i$u_i + i$alpha6_i * pred_range_standardized)
      )
    }) %>% bind_rows()
  predictions = predictions %>% left_join(groups %>% select(common_name, group), by = "common_name")
  
  mu_u_samples      = as.matrix(model_data[[m]]$msom$sims.list[["mu.u"]])
  mu_param_samples  = as.matrix(model_data[[m]]$msom$sims.list[[paste0("mu.", param_data$param)]])
  n_groups = ncol(mu_u_samples)
  
  meta_preds = map_dfr(1:nrow(mu_u_samples), function(i) {
    map_dfr(1:n_groups, function(g) {
      tibble(
        group_idx = g,
        idx = pred_range_original,
        psi = plogis(
          mu_u_samples[i, g] +
            mu_param_samples[i, g] * pred_range_standardized
        ),
        draw = i
      )
    })
  })
  meta_summary = meta_preds %>% # calculate means and 95% BCIs
    group_by(group_idx, idx) %>% summarize(
      psi_mean = mean(psi),
      psi_lower = quantile(psi, 0.025),
      psi_upper = quantile(psi, 0.975)
    )
  meta_summary = meta_summary %>%
    left_join(
      groups %>% select(group_idx, group) %>% distinct(),
      by = "group_idx"
    )
  
  # ggplot() +
  #   geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  #   geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  #   scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  #   labs(x = param_data$name, y = "Occurrence probability") +
  #   theme_sleek()
  
  plots_marg[[m]] = ggplot(predictions %>% filter(group != "other"), aes(x = idx, y = psi, group = common_name)) +
    geom_line(aes(color = group), alpha = 0.3) +
    geom_ribbon(data = meta_summary %>% filter(group != "other"), aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group), alpha = 0.2, inherit.aes = FALSE) +
    geom_line(data = meta_summary %>% filter(group != "other"), aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
    facet_wrap(~str_to_title(group)) +
    scale_x_continuous(limits = c(bound_low, bound_high)) +
    scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0)) + theme(
      legend.position = "none"
    )
}

(plots_marg[["all"]] + labs(x = "B-IBI", y = "Occupancy probability") +
  theme(legend.position = "none") + scale_color_manual(values = c("gray20")) + scale_fill_manual(values = c("gray20"))) +
(plots_marg[["diet"]] + labs(x = "B-IBI", y = "") +
  theme(legend.position = "none") + scale_color_manual(values = c("orange")) + scale_fill_manual(values = c("orange"))) +
(plots_marg[["forage"]] + labs(x = "B-IBI", y = "") +
  theme(legend.position = "none") + scale_color_manual(values = rep("gray20", 4)) + scale_fill_manual(values = rep("gray20", 4))) +
(plots_marg[["migrant"]] + labs(x = "B-IBI", y = "") +
  theme(legend.position = "none") + scale_color_manual(values = c("orange")) + scale_fill_manual(values = c("orange")))
