library(tidyverse)

path = "data/cache/models/reach_all.rds"
in_path_avonet_traits   = "data/traits/AVONET Supplementary dataset 1.xlsx"
in_path_species_list    = "data/pam/species_list.txt"

message("Loading data for model ", path)
model_data = readRDS(path)

species_names = read_lines(in_path_species_list) %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>% mutate(
    common_name = tolower(common_name), scientific_name = tolower(scientific_name)
  ) %>% filter(common_name %in% sort(unique(model_data$species)))

avonet = readxl::read_xlsx(in_path_avonet_traits, sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  mutate(scientific_name = tolower(scientific_name)) %>%
  filter(scientific_name %in% species_names$scientific_name) %>%
  select(scientific_name, family, order, mass, habitat, habitat_density, migration, trophic_level, trophic_niche, primary_lifestyle)
species_metadata = left_join(species_names, avonet, by = "scientific_name")

community_baselines = model_data$msom_summary %>%
  filter(param %in% c("mu.u", "sigma.u",
                      "mu.v", "sigma.v",
                      "mu.w", "sigma.w",
                      "mu.b", "sigma.b")) %>%
  select(param, prob, prob_lower95, prob_upper95)

message("Baseline occurrence probability:")
print(subset(community_baselines, param == "mu.u"))

message("Baseline detection probability:")
print(subset(community_baselines, param == "mu.v"))

ggplot(community_baselines, aes(x = prob, y = param)) +
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
  mutate(species_name = species[species_idx])

message("Species occurrence probability range:")
species_baselines %>% filter(startsWith(param, "u")) %>% summarise(
  prob_min = min(prob),
  species_min = species_name[which.min(prob)],
  prob_max = max(prob),
  species_max = species_name[which.max(prob)],
  .groups = "drop"
)

ggplot(species_baselines %>% filter(startsWith(param, "u")), aes(x = prob, y = reorder(species_name, prob))) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbarh(aes(xmin = prob_lower95, xmax = prob_upper95), height = 0.1, position = position_dodge(width = 0.5)) +
  labs(x = "Posterior probability", y = "Species index",
       title = "Species-specific occurrence probability `u` across models (mean and 95% BCI)"
  ) + theme(legend.position = "bottom")

ggplot(species_baselines %>% filter(startsWith(param, "v")), aes(x = prob, y = reorder(species_name, prob))) +
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

occ_coef_summary = model_data$msom_summary %>%
  filter(str_detect(param, "^(mu|sigma)\\.alpha")) %>% arrange(param) %>%
  select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

detect_coef_summary = model_data$msom_summary %>%
  filter(str_detect(param, "^(mu|sigma)\\.beta")) %>% arrange(param) %>%
  select(param, mean, sd, `2.5%`, `97.5%`, `25%`, `75%`, overlap0)

param_alpha_data = model_data$param_alpha_data
param_beta_data = model_data$param_beta_data
  
param_occ_data = param_alpha_data

occ_effect_sizes = full_join(occ_coef_summary %>% filter(str_starts(param, "mu")), param_occ_data %>% mutate(param = paste0("mu.", param)), by='param')

# Compare species level effects of each covariate on occurrence
ggplot(occ_effect_sizes, aes(x = mean, y = as.factor(name))) +
    geom_vline(xintercept = 0, color = "gray") +
    geom_point(aes(color = overlap0)) +
    geom_errorbar(aes(xmin = `25%`,  xmax = `75%`,   color = overlap0), width = 0, linewidth = 1) +
    geom_errorbar(aes(xmin = `2.5%`, xmax = `97.5%`, color = overlap0), width = 0) +
    scale_color_manual(values = c("black", "gray")) +
    xlim(c(-1, 1)) +
    labs(title = "Community level effect sizes for occurrence covariates",
         subtitle = tools::file_path_sans_ext(basename(path)),
         x = "Coefficient estimate", y = "Parameter")

species_effects = param_occ_data %>%
  mutate(coef = map(param, ~ model_data$msom_summary %>% filter(str_detect(param, paste0("^", .x, "(?!\\d)\\["))))) %>%
  unnest(coef, names_sep = "_") %>% mutate(species_idx = as.integer(gsub(".*\\[(\\d+)\\]", "\\1", coef_param))) %>% mutate(common_name = species[species_idx])

species_effects = species_effects %>% left_join(species_metadata, by ="common_name")

# TODO: join with species metadata

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
  geom_beeswarm(data = species_effects, aes(x = coef_mean, y = name, color = primary_lifestyle), cex = 2, priority = "density") +
  geom_errorbar(data = occ_effect_sizes, aes(x = mean, y = as.factor(name), xmin = `2.5%`, xmax = `97.5%`), width = 0) +
  geom_point(data = occ_effect_sizes, aes(x = mean, y = as.factor(name)), size = 3.5)


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

