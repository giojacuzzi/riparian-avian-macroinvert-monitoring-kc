source("src/global.R")

path = "data/cache/models/reach_invert_predator.rds"

message("Loading data for model ", path)
model_data = readRDS(path)
groups = model_data$groups

community_baselines = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.u", "sigma.u", "mu.v", "sigma.v",
      "mu.w", "sigma.w", "mu.b", "sigma.b"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, prob, prob_lower95, prob_upper95)

community_baselines = community_baselines %>% left_join(groups %>% distinct(group_idx, group), by = "group_idx")

message("Baseline occurrence probability:")
print(subset(community_baselines, param == "mu.u"))

message("Baseline detection probability:")
print(subset(community_baselines, param == "mu.v"))

ggplot(community_baselines, aes(x = prob, y = param, color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = prob_lower95, xmax = prob_upper95), width = 0.2, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Parameter",
       title = "Hyperparameter estimates across models (mean and 95% BCI)",
       subtitle = "
u - occupancy probability across sites
v - true positive detection probability given presence
") + theme(legend.position = "bottom")

species = model_data$species

species_baselines = model_data$msom_summary %>%
  filter(grepl("^[uvwb]\\[", param)) %>%
  select(param, prob, prob_lower95, prob_upper95) %>%
  mutate(species_idx = str_extract(param, "(?<=\\[)[0-9]+") %>% as.integer()) %>%
  mutate(common_name = species[species_idx]) %>%
  left_join(groups, by = "common_name")

message("Species occurrence probability range:")
species_baselines %>% filter(startsWith(param, "u")) %>% summarise(
  prob_min = min(prob),
  species_min = common_name[which.min(prob)],
  prob_max = max(prob),
  species_max = common_name[which.max(prob)],
  .groups = "drop"
)

ggplot(species_baselines %>% filter(startsWith(param, "u")), aes(x = prob, y = reorder(common_name, prob), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific occurrence probability `u` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

ggplot(species_baselines %>% filter(startsWith(param, "v")), aes(x = prob, y = reorder(common_name, prob), color = group)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific TP detection probability given presence `v` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

# TODO
# Nsite_posterior = model_data$msom_summary %>% filter(stringr::str_starts(param, "Nsite")) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
#   mutate(site_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(site = sites[site_idx])
# Nsite_mean = mean(Nsite_posterior$mean)
# message("Mean estimated species richness across all sites: ", round(Nsite_mean,1), " (range ", round(min(Nsite_posterior$mean),1), "–", round(max(Nsite_posterior$mean),1), ")")
# ggplot(Nsite_posterior, aes(x = as.factor(plot_order), y = mean)) +
#   geom_hline(yintercept = Nsite_mean, linetype = "solid", color = "blue") +
#   geom_point() + geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0) +
#   scale_x_discrete(labels = Nsite_posterior$site) + 
#   labs(title = "Estimated species richness per site", x = "Site", y = "Estimated species richness") +
#   coord_flip()
# ggplot(Nsite_posterior, aes(x = mean)) +
#   geom_histogram(binwidth = 1) +
#   geom_vline(xintercept = Nsite_mean, color = "blue") +
#   labs(title = "Estimated species richness per site", x = "Number of species estimated", y = "Number of sites")

# Compare coefficients on occupancy -------------------------------------------------------------------------

message("Comparing model coefficients")

# occ_coef_summary = model_data$msom_summary %>%
#   filter(str_detect(param, "^(mu|sigma)\\.alpha")) %>% arrange(param) %>%
#   select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

occ_coef_summary = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.alpha", "sigma.alpha"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
  left_join(groups %>% distinct(group_idx, group), by = "group_idx")


# detect_coef_summary = model_data$msom_summary %>%
#   filter(str_detect(param, "^(mu|sigma)\\.beta")) %>% arrange(param) %>%
#   select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

detect_coef_summary = model_data$msom_summary %>%
  filter(Reduce(`|`, lapply(
    c("mu.beta", "sigma.beta"),
    \(p) str_starts(param, p)
  ))) %>%
  mutate(
    group_idx = str_extract(param, "(?<=\\[)\\d+(?=\\])") %>% as.integer(),
    param = str_remove(param, "\\[\\d+\\]")
  ) %>%
  select(param, group_idx, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0) %>%
  left_join(groups %>% distinct(group_idx, group), by = "group_idx")

param_alpha_data = model_data$param_alpha_data
param_beta_data = model_data$param_beta_data
  
param_occ_data = param_alpha_data

occ_effect_sizes = full_join(occ_coef_summary %>% filter(str_starts(param, "mu")), param_occ_data %>% mutate(param = paste0("mu.", param)), by='param')

occ_effect_sizes = occ_effect_sizes %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

# Compare species level effects of each covariate on occurrence
ggplot(occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point() +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`), width = 0, linewidth = 1) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`), width = 0) +
    # scale_color_manual(values = c("black", "gray")) +
    # xlim(c(-1, 1)) +
    labs(title = "Community level effect sizes for occurrence covariates",
         subtitle = tools::file_path_sans_ext(basename(path)),
         x = "Coefficient estimate", y = "Parameter")

species_effects = param_occ_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])

species_effects = species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")

species_effects = species_effects %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

plt = ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  # geom_point(data = species_effects, aes(x = coef_mean, y = name, shape = coef_overlap0), position = position_jitter(height = 0.2), size = 3, alpha = 0.95) +
  # geom_point(data = species_effects, aes(x = coef_mean, y = name, color = habitat_association, shape = coef_overlap0), position = position_jitter(height = 0.2), size = 3, alpha = 0.95) +
  geom_dotplot(data = species_effects, aes(x = coef_mean, y = name, fill = name, alpha = coef_overlap0), color = "transparent", binwidth = 0.03, stackdir = "center") +
  # scale_color_manual(values = c("darkgray", "#d8c18a", "#b2675e", "#3c8273")) +
  # scale_color_viridis_c() +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), shape = overlap0), size = 3.5) +
  # xlim(c(-2.5, 2.5)) +
  labs(title = "Effect sizes for occurrence covariates at community and species levels",
       subtitle = tools::file_path_sans_ext(basename(path)),
       x = "Coefficient estimate", y = "Parameter") +
  geom_text_repel(
    data = species_effects %>% filter(coef_overlap0 == 0),
    aes(x = coef_mean, y = name, label = common_name),
    size = 3, nudge_x = 0.05, direction = "y", hjust = 0.05, max.overlaps = 10
  ) + theme(legend.position = "bottom"); print(plt)

library(ggbeeswarm)
plt = ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_quasirandom(data = species_effects, aes(x = coef_mean, y = name, color = coef_overlap0),
    size = 3, width = 0.3, alpha = 0.75
  ) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), shape = overlap0), size = 3.5) +
  # xlim(c(-2.5, 2.5)) +
  labs(title = "Effect sizes for occurrence covariates at community and species levels",
       subtitle = tools::file_path_sans_ext(basename(path)),
       x = "Coefficient estimate", y = "Parameter") +
  geom_text_repel(
    data = species_effects,
    aes(x = coef_mean, y = name, label = common_name),
    size = 2, nudge_x = 0.02, direction = "y", hjust = 0.02, max.overlaps = 30
  ) + theme(legend.position = "bottom"); print(plt)

ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_beeswarm(data = species_effects, aes(x = coef_mean, y = name, color = group, shape = coef_overlap0), cex = 2, priority = "density") +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group)), size = 3.5) +
  scale_shape_manual(values = c(16, 1)) +
  labs(color = "Predator", shape = "Significance")
  geom_text_repel(
    data = species_effects,
    aes(x = coef_mean, y = name, label = common_name, color = group),
    size = 2, nudge_x = 0.02, direction = "y", hjust = 0.05, max.overlaps = 30
  )
  
  
label_species <- c("black-headed grosbeak", "belted kingfisher", "swainson's thrush", "western wood-pewee", "wilson's warbler", "house sparrow", "bewick's wren", "pacific wren", "black-throated gray warbler", "pileated woodpecker")
  
# ===========================================================================================================
p_occ = ggplot(species_effects, aes(x = coef_mean, y = name, group = group_label, color = group_label)) +
  geom_vline(xintercept = 0, color = "gray80") +
  geom_vline(xintercept = 0.5, color = "gray90", linetype = "dashed") +
  geom_vline(xintercept = -0.5, color = "gray90", linetype = "dashed") +
  geom_beeswarm(aes(shape = coef_overlap0), dodge.width = 0.5, cex = 1, priority = "density", size = 2, alpha = 0.6) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), xmin = `2.5%`, xmax = `97.5%`), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), xmin = `25%`, xmax = `75%`), size = 1, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), group = as.factor(group_label)), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("dodgerblue2", "gray20")) +
  labs(x = "Occupancy coefficient mean", y = "Occupancy predictor", color = "Diet group") +
  theme_sleek() + guides(shape = "none"); print(p_occ)


param_detect_data = param_beta_data

detect_effect_sizes = full_join(detect_coef_summary %>% filter(str_starts(param, "mu")), param_detect_data %>% mutate(param = paste0("mu.", param)), by='param')
detect_effect_sizes = detect_effect_sizes %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

detect_species_effects = param_detect_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])
detect_species_effects = detect_species_effects %>% left_join(species_traits, by ="common_name") %>%
  left_join(groups, by = "common_name")
detect_species_effects = detect_species_effects %>% mutate(group_label = ifelse(group == "invert_predator", "Invertivore", "Other species"))

p_detect = ggplot(detect_species_effects, aes(x = coef_mean, y = name, group = group_label, color = group_label)) +
  geom_vline(xintercept = 0, color = "gray80") +
  # geom_vline(xintercept = 0.5, color = "gray90", linetype = "dashed") +
  # geom_vline(xintercept = -0.5, color = "gray90", linetype = "dashed") +
  geom_beeswarm(aes(shape = coef_overlap0), dodge.width = 0.5, cex = 1, priority = "density", size = 2, alpha = 0.6) +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), xmin = `2.5%`, xmax = `97.5%`), size = 0.5, width = 0, position = position_dodge(width = 0.5)) +
  geom_errorbar(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), xmin = `25%`, xmax = `75%`), size = 1, width = 0, position = position_dodge(width = 0.5)) +
  geom_point(data = detect_effect_sizes, aes(x = mean, y = as.factor(name), color = as.factor(group_label), group = as.factor(group_label)), size = 3, position = position_dodge(width = 0.5)) +
  scale_shape_manual(values = c(19, 1)) +
  scale_color_manual(values = c("dodgerblue2", "gray20")) +
  labs(x = "Detection coefficient mean", y = "Detection predictor", color = "Diet group") +
  theme_sleek() + guides(shape = "none"); print(p_detect)

##################

param_name = "alpha2"
param_data   = param_alpha_data %>% filter(param == param_name)

# Mean marginal probabilities of occurrence for the metacommunity in relation to alpha
coefs = msom_summary %>% filter(stringr::str_starts(param, param_name)) %>% arrange(mean) %>% mutate(plot_order = 1:nrow(.)) %>%
  mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", param))) %>% mutate(common_name = species[species_idx])
param_scaled_data = param_data %>% pull(scaled)
p_mu = attr(param_scaled_data[[1]], "scaled:center") # to transform between z-scale and pred_range_original scale
p_sd = attr(param_scaled_data[[1]], "scaled:scale")
bound_low  = min(param_scaled_data[[1]]) * p_sd + p_mu
bound_high = max(param_scaled_data[[1]]) * p_sd + p_mu
pred_range_original = seq(bound_low, bound_high, by = 1) # range of possible alpha values
pred_range_standardized = (pred_range_original - p_mu) / p_sd
# species-specific occurrence intercepts u[i]
intercepts = msom_summary %>% filter(str_starts(param, "u")) %>%
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
meta_summary$group = factor(meta_summary$group_idx)
# meta_summary = meta_summary %>% left_join(groups %>% distinct(group_idx, group, grouping), by = "group_idx")

p_bibi = ggplot() +
  # geom_line(data = predictions, aes(x = idx, y = psi, group = common_name, color = group), alpha = 0.2) +
  geom_ribbon(data = meta_summary, aes(x = idx, ymin = psi_lower, ymax = psi_upper, fill = group, group = group), alpha = 0.2, inherit.aes = FALSE) +
  # geom_line(data = meta_summary, aes(x = idx, y = psi_mean, color = group, group = group), linewidth = 1.2, inherit.aes = FALSE) +
  # scale_x_continuous(limits = c(bound_low, bound_high)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.25, 0.5, 0.75, 1.0)) +
  labs(x = param_data$name, y = "Occurrence probability") +
  theme_sleek(); print(p_bibi)





ggplot(predictions, aes(x = mango, y = psi, group = common_name)) +
  geom_line() +
  gghighlight::gghighlight(use_direct_label = FALSE, unhighlighted_params = list(color = alpha("black", 0.1))) +
  facet_wrap(~reorder(common_name, psi)) +
  scale_x_continuous(limits = c(bound_low, bound_high)) +
  scale_y_continuous(limits = c(0.0, 1.0), breaks = c(0, 0.5, 1.0)) +
  labs(title = "", x = param_data$name, y = "Occurrence probability") +
  theme(legend.position = "none")













############################

ggplot(species_effects %>%
         mutate(group_shape = paste(group, coef_overlap0, sep = "_")), aes(x = coef_mean, y = name, color = group, group = group_shape)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_beeswarm(aes(shape = coef_overlap0), dodge.width = 0.5, cex = 1, priority = "density", size = 1.5)


ggplot(data = species_effects %>% filter(param == "alpha1"),
       aes(x = coef_mean, y = reorder(common_name, coef_mean))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(aes(color = coef_overlap0)) + labs(title = "bibi")

ggplot(data = species_effects %>% filter(param == "alpha2"),
       aes(x = coef_mean, y = reorder(common_name, coef_mean))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(aes(color = coef_overlap0)) + labs(title = "canopy")

ggplot(data = species_effects %>% filter(param == "alpha3"),
       aes(x = coef_mean, y = reorder(common_name, coef_mean))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
  geom_errorbar(aes(xmin = `coef_2.5%`, xmax = `coef_97.5%`)) +
  geom_point(aes(color = coef_overlap0)) + labs(title = "imp")

