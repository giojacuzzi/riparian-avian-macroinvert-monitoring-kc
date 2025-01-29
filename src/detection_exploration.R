# Calculate naive species richness per site from detection histories

# INPUT
path_input_dir = "data/processed/predicted_detection_histories"

path_site_data = "/Users/Steve/Documents/site_data.csv"
############################################################################################################


theme_set(theme_classic())

paths_site = list.files(path_input_dir, full.names = TRUE)

# Calculate species presence and naive alpha diversity (species richness) by site
message('Calculating species presence and naive alpha diversity (species richness) by site...')
site_richness = tibble()
site_species = tibble()

for (path_site in paths_site) {
  # Load site detection history
  site_detections = read_csv(path_site, show_col_types = FALSE)
  site_id = sub('.csv', '', basename(path_site))
  message('Site ', site_id)
  # Determine which species were detected during at least one visit at the site
  species_detected = site_detections %>% filter(rowSums(select(., -species)) >= 1) %>% select(species)
  site_species = bind_rows(site_species, tibble(
    site = site_id,
    species = as.list(species_detected)
  ))
  # Calculate richness
  site_richness = bind_rows(site_richness, tibble(
    site  = site_id,
    alpha = nrow(species_detected)
  ))
}

message('Species by site:')
for (i in 1:nrow(site_species)) {
  message(site_species$site[[i]], ' (', length(site_species$species[[i]]), ')')
  print(site_species$species[[i]])
}
message('Alpha diversity (species richness) by site:')
print(site_richness)
p = ggplot(site_richness, aes(x = alpha)) +
  geom_histogram() +
  labs(title='Alpha diversity (species richness) distribution across sites'); print(p)

############################################################################################################
readline('Proceed with site data analysis? [enter]')

# Load data per site
site_data = na.omit(read_csv(path_site_data, show_col_types = FALSE))

# Join with species richness data
site_data = left_join(site_richness, site_data, by = 'site')


# Visualize impervious coverage across sites
p = ggplot(site_data, aes(x = Impervious)) +
  geom_histogram() +
  labs(title='Impervious coverage across sites'); print(p)

message('Analyzing impervious coverage x mean B-IBI...')

# Calculate Pearson's correlation coefficient for impervious coverage x mean B-IBI
correlation_pearson = cor(site_data$Impervious, site_data$mean_BIBI, method = "pearson")
message("Pearson's correlation coefficient", round(correlation_pearson, 2))

# Fit a linear regression model
model = lm(site_data$mean_BIBI ~ site_data$Impervious)
summary(model)

# Visualize impervious coverage x mean B-IBI
p = ggplot(site_data, aes(x = Impervious, y = mean_BIBI)) +
  geom_point() +
  geom_line(aes(y = model$fitted.values)) + # linear regression
  labs(title='Impervious coverage x mean B-IBI'); print(p)

# TODO: Explore relationship between impervious coverage and species richness (functional riparian insectivore guild)
ggplot(site_data, aes(x = Impervious, y = alpha)) +
  geom_smooth(method='lm') +
  geom_point()

correlation_pearson = cor(site_data$Impervious, site_data$alpha, method = "pearson")
message("Pearson's correlation coefficient ", round(correlation_pearson, 2))

# TODO: Explore relationship between B-IBI and species richness (functional riparian insectivore guild)

ggplot(riparian_BIBI, aes(x = mean_BIBI, y = species_count)) + 
  geom_point() +
  geom_smooth(method = "lm")
correlation_pearson = cor(riparian_BIBI$species_count, riparian_BIBI$mean_BIBI, method = "pearson")
message("Pearson's correlation coefficient ", round(correlation_pearson, 2))
model <-lm(species_count ~ mean_BIBI, data = riparian_BIBI)
summary(model)

# Land Data vs Riparian association

ggplot(riparian_imp, aes(x = totalimp, y = species_count)) +
  geom_point() +
  geom_smooth(method = "lm")
correlation_pearson = cor(riparian_imp$species_count, riparian_imp$totalimp, method = "pearson")
message("Pearson's correlation coefficient ", round(correlation_pearson, 2))

# Multiple regression model
riparian_bibi_imp <- riparian_BIBI %>%
  left_join(riparian_imp, by = "site") %>%
  rename(species_count = species_count.x) %>%
  select(-species_count.y)


model <-lm(species_count.x ~ mean_BIBI + totalimp, data = riparian_bibi_imp)
summary(model)

# Stacked barchart

# library
library(ggplot2)


# Stacked barchart

riparian_bibi_imp <- riparian_bibi_imp %>%
  mutate(bibi_category = case_when(
    mean_BIBI >= 0 & mean_BIBI <= 40 ~ "low",
    mean_BIBI > 40 & mean_BIBI <= 60 ~ "mid",
    mean_BIBI > 60 & mean_BIBI <= 100 ~ "high")) %>%
  mutate(imp_category = case_when(
    totalimp >= 0 & totalimp <= 10 ~ "sensitive",
    totalimp > 10 & totalimp <= 25 ~ "impacted",
    totalimp > 25 & totalimp <= 100 ~ "non-supporting"
  ))

barchart_table <- filtered_species_split %>%
  group_by(site, species) %>%
  summarize(present = 1) %>%
  pivot_wider(names_from = site, values_from = present, values_fill = list(present = 0))

barchart_table <- barchart_table %>%
  pivot_longer(cols = -species, names_to = "site")


species_counts <- barchart_table %>%
  left_join(riparian_bibi_imp, by = "site")

species_counts$bibi_category <- factor(species_counts$bibi_category, 
                                       levels = c("high", "mid", "low"))
species_counts$imp_category <- factor(species_counts$imp_category, 
                                       levels = c("sensitive", "impacted", "non-supporting"))

species_counts$species <- factor(species_counts$species, 
                                 levels = species_counts %>%
                                   group_by(species) %>%
                                   summarise(total_count = sum(value)) %>%
                                   arrange(desc(total_count)) %>%
                                   pull(species))

# BIBI Categories
ggplot(species_counts, aes(y = species, x = value, fill = bibi_category)) +
  geom_bar(stat = "identity", position = "stack") + labs(title = "Species Counts by Site, BIBI", x = "Number of sites detected", y = "Species") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))

# Impervious Categories
ggplot(species_counts, aes(y = species, x = value, fill = imp_category)) +
  geom_bar(stat = "identity", position = "stack") + labs(title = "Species Counts by Site, Impervious Percent", x = "Number of sites detected", y = "Species") +
  theme(axis.text.y = element_text(angle = 0, hjust = 1))
