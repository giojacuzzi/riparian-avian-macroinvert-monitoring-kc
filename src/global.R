# TODO: Load all land cover and benthos data for each site and join with the site metadata to create a "master" data_habitat dataframe

path_land_cover = "data/processed/land_cover/land_cover.csv" # processed land cover data from GIS
path_benthos_raw = "data/raw/benthos/ScoresByYear.txt" # raw PSSB benthos data

data_land_cover = read.csv(path_land_cover)
data_benthos <- read.delim(path_benthos_raw, header = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE)

metadata_site = read.csv('data/site_metadata.csv')
  
data_habitat = full_join(data_land_cover, metadata_site, by = c("site" = "aru_site"))
data_habitat = left_join(data_habitat, data_benthos, by = c("benthos_site" = "Site.ID"))

data_riparianobligate = read_csv('data/processed/Species_Habitat_List.csv')

# Create count of riparian dependent insectivores in each site
riparian_species <- data_riparianobligate %>%
  filter(riparian_dependent_breeding == "Yes") %>%
  filter(insectivore == "Yes") %>%
  select(species)

site_species_split <- unnest(site_species, cols = "species")

filtered_species_split <- site_species_split %>%
  filter(species %in% riparian_species$species)

site_riparian_count <- filtered_species_split %>%
  group_by(site) %>%
  summarize(species_count = n_distinct(species))

# Table with riparian dependent count and mean BIBI by site

site_bibi <- site_data %>%
  select(site, mean_BIBI)

riparian_BIBI <- site_riparian_count %>%
  left_join(site_bibi, by = "site")

# Table with riparian dependent count and total imp percent

data_land_cover <- data_land_cover %>%
  mutate(totalimp = Impervious + Impervious..Covered.by.Tree.Canopy)

site_imp <- data_land_cover %>%
  select(site, totalimp)

riparian_imp <- site_riparian_count %>%
  left_join(site_imp, by = "site")
