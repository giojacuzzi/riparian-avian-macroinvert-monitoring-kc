# Calculate naive species richness per site from detection histories

# INPUT
path_input_dir = "data/processed/predicted_detection_histories"

############################################################################################################

library(tidyverse)

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
