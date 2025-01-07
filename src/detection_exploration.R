# Calculate naive species richness per site from detection histories

# INPUT
path_input_dir = "data/processed/predicted_detection_histories"

############################################################################################################

paths_site = list.files(path_input_dir, full.names = TRUE)

# Calculate naive species richness by site
message('Calculating naive species richness by site...')
site_richness = tibble()

for (path_site in paths_site) {
  # Load site detection history
  site_detections = read_csv(path_site, show_col_types = FALSE)
  site_id = basename(path_site)
  message('Site ', site_id)
  # Determine which species were detected during at least one visit at the site
  species_detected = site_detections %>% filter(rowSums(select(., -species)) >= 1) %>% select(species)
  # Calculate richness
  site_richness = bind_rows(site_richness, tibble(
    site  = site_id,
    alpha = nrow(species_detected)
  ))
}
message('Species richness by site:')
print(site_richness)
