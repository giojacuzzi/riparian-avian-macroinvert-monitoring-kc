## Simple Capstone SEM

library(piecewiseSEM)

richness_BIBI_landcover

urban_bibi = lm(mean_BIBI ~ totalimp, richness_BIBI_landcover)
urban_bibi_bird = lm(species_count ~ mean_BIBI + totalimp, richness_BIBI_landcover)

model = psem(urban_bibi, urban_bibi_bird, data = richness_BIBI_landcover)
model
plot(model)
fisherC(model)

## Validated Data ##
site_ripoblig_validated

urban_bibi = lm(mean_BIBI.y ~ totalimp, site_ripoblig_validated)
urban_bibi_bird = lm(species_count_val ~ mean_BIBI.y + totalimp, site_ripoblig_validated)
model = psem(urban_bibi, urban_bibi_bird, data = site_ripoblig_validated)
plot(model)


