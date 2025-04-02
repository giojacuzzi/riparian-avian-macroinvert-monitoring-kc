# TODO: Load all land cover and benthos data for each site and join with the site metadata to create a "master" data_habitat dataframe
library(tidyverse)
library(ggplot2)

#SET UP PATHS#
path_input_dir = "data/processed/predicted_detection_histories"
path_site_data = "/Users/Steve/Documents/site_data.csv"
paths_site = list.files(path_input_dir, full.names = TRUE)
path_land_cover = "data/processed/land_cover/land_cover.csv" # processed land cover data from GIS
path_benthos_raw = "data/raw/benthos/ScoresByYear.txt" # raw PSSB benthos data

# Calculate species presence and naive alpha diversity (species richness) by site
message('Calculating species presence and naive alpha diversity (species richness) by site...')
site_richness = tibble()
site_species = tibble()

# Read all predictions and load all species seen at each site #
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

# Load data per site and join with total species richness
site_data = na.omit(read_csv(path_site_data, show_col_types = FALSE))

site_data = left_join(site_richness, site_data, by = 'site')


#SET UP DATA#
data_land_cover = read.csv(path_land_cover) #LAND COVER DATA
data_benthos <- read.delim(path_benthos_raw, header = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE) #BENTHOS DATA
metadata_site = read.csv('data/site_metadata.csv') #SITE DATA
data_habitat = full_join(data_land_cover, metadata_site, by = c("site" = "aru_site")) #JOIN SITE WITH LAND COVER
data_habitat = left_join(data_habitat, data_benthos, by = c("benthos_site" = "Site.ID")) #JOIN SITE & LAND COVER WITH BENTHOS


#LOAD SPECIES LIST#
data_specieslist = read_csv('data/processed/Species_Habitat_List.csv')

#Create list of riparian-dependent insectivores#
riparian_species <- data_specieslist %>%
  filter(riparian_dependent_breeding == "Yes") %>%
  filter(insectivore == "Yes") %>%
  select(species)

#Create list of all non-insectivore species#
noninsectiv <- data_specieslist %>%
  filter(is.na(insectivore)) %>%
  select(species) %>%
  filter(!is.na(species))

##COUNT OF RIPARIAN SPECIES PER SITE, Insectivores in one and Non-insectivores in another##
site_species_split <- unnest(site_species, cols = "species")

filtered_species_split <- site_species_split %>%
  filter(species %in% riparian_species$species)

noninsectfiltered_species_split <- site_species_split %>%
  filter(species %in% noninsectiv$species) 

site_riparian_count <- filtered_species_split %>%
  group_by(site) %>%
  summarize(species_count = n_distinct(species))

site_noninsect_count <- noninsectfiltered_species_split %>%
  group_by(site) %>%
  summarize(species_count = n_distinct(species))

# Join Mean BIBI and Riparian Count #

site_bibi <- site_data %>%
  select(site, mean_BIBI)

riparian_BIBI <- site_riparian_count %>%
  left_join(site_bibi, by = "site")


# Table with riparian dependent count and total imp percent #

data_land_cover <- data_land_cover %>%
  mutate(totalimp = Impervious + Impervious..Covered.by.Tree.Canopy)

site_imp <- data_land_cover %>%
  select(site, totalimp)

riparian_imp <- site_riparian_count %>%
  left_join(site_imp, by = "site")

# Table with all noninsectivores #

noninsect_BIBI <- site_noninsect_count %>%
  left_join(site_bibi, by = "site")

noninsect_imp <- site_noninsect_count %>%
  left_join(site_imp, by = "site")

## Master riparian species count, mean bibi, and habitat data table ##
riparian_bibi_imp <- riparian_BIBI %>%
  left_join(riparian_imp, by = "site") %>%
  rename(species_count = species_count.x) %>%
  select(-species_count.y)

richness_BIBI_landcover <- riparian_BIBI %>%
  left_join(data_land_cover, by = "site")
