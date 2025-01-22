# TODO: Load all land cover and benthos data for each site and join with the site metadata to create a "master" data_habitat dataframe

path_land_cover = "data/processed/land_cover/land_cover.csv" # processed land cover data from GIS
path_benthos_raw = "data/raw/benthos/ScoresByYear.txt" # raw PSSB benthos data

data_land_cover = read.csv(path_land_cover)
data_benthos <- read.delim(path_benthos_raw, header = TRUE, na.strings = c("", "NA"), stringsAsFactors = FALSE)

metadata_site = read.csv('data/site_metadata.csv')
  
data_habitat = full_join(data_land_cover, metadata_site, by = c("site" = "aru_site"))
data_habitat = left_join(data_habitat, data_benthos, by = c("benthos_site" = "Site.ID"))

data_riparianobligate = read_csv('data/processed/Species_Habitat_List.csv')
