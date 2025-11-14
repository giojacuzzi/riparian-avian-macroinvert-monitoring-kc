library(tidyverse)
library(ggrepel)
library(janitor)
library(sf)
library(DHARMa)
library(mapview)
library(piecewiseSEM)
theme_set(theme_minimal())

in_cache_detections     = "data/cache/1_preprocess_agg_pam_data/detections_calibrated_0.5.rds" # detections_calibrated_0.75.rds
in_path_species_list    = "data/pam/species_list.txt"
in_path_avonet_traits   = "data/traits/AVONET Supplementary dataset 1.xlsx"
in_path_elton_traits    = "data/traits/EltonTraits 1.0/BirdFuncDat.csv"

exclude_agri_sites = FALSE

# Load site variable data ------------------------------------------------------------

site_data_reach = readRDS("data/cache/3_calculate_vars/NEW_site_data_550m.rds")
site_data_basin = readRDS("data/cache/3_calculate_vars/NEW_site_data_5000m.rds")

# Urbanization x B-IBI gradient
ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = bibi)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.01) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "orange", alpha = 0.01) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "yellow", alpha = 0.01) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "green", alpha = 0.01) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "dodgerblue", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Canopy x B-IBI gradient
ggplot(site_data_reach, aes(x = rast_usfs_canopycover_sum_proportion, y = bibi)) +
  geom_rect(aes(ymin = 0, ymax = 20, xmin = -Inf, xmax = Inf), fill = "red", alpha = 0.01) +
  geom_rect(aes(ymin = 20, ymax = 40, xmin = -Inf, xmax = Inf), fill = "orange", alpha = 0.01) +
  geom_rect(aes(ymin = 40, ymax = 60, xmin = -Inf, xmax = Inf), fill = "yellow", alpha = 0.01) +
  geom_rect(aes(ymin = 60, ymax = 80, xmin = -Inf, xmax = Inf), fill = "green", alpha = 0.01) +
  geom_rect(aes(ymin = 80, ymax = 100, xmin = -Inf, xmax = Inf), fill = "dodgerblue", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Urbanization x Canopy gradient
ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = rast_usfs_canopycover_sum_proportion)) +
  geom_rect(aes(xmin = 0.0, xmax = 0.20, ymin = -Inf, ymax = Inf), fill = "pink", alpha = 0.01) +
  geom_rect(aes(xmin = 0.20, xmax = 0.50, ymin = -Inf, ymax = Inf), fill = "tomato", alpha = 0.01) +
  geom_rect(aes(xmin = 0.50, xmax = 0.80, ymin = -Inf, ymax = Inf), fill = "red", alpha = 0.01) +
  geom_point() + geom_text_repel(aes(label = site_id))

# Load species detection history data ------------------------------------------------
message("Loading species detection history data")

# Load species detection history data
detections = readRDS(in_cache_detections)
detections$long$common_name = tolower(detections$long$common_name)
detections$wide$common_name = tolower(detections$wide$common_name)

species_names = read_lines(in_path_species_list) %>% as_tibble() %>%
  separate(value, into = c("scientific_name", "common_name"), sep = "_") %>% mutate(
    common_name = tolower(common_name), scientific_name = tolower(scientific_name)
  ) #%>% filter(common_name %in% sort(unique(detections$long$common_name)))

# TODO: assess adequacy of sampling effort to detect total species richness via sample size-based rarefaction and extrapolation sampling curves for species richness (iNEXT)? Quantify sample completeness to determine the estimated proportion of species detected from the predicted species pool by plotting sample coverage with respect to the number of sampling units.
if (FALSE) {
  library(iNEXT)
  
  zero_detections <- detections$long %>%
    group_by(common_name) %>%
    summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
    filter(total_detections == 0) %>%
    pull(common_name)
  
  detections_list <- detections$wide %>%
    group_split(site_id) %>%                                 # Split by site
    set_names(unique(detections$wide$site_id)) %>%            # Name list elements by site_id
    map(~ {
      df <- .x %>%
        select(-site_id) %>%                                  # Remove site_id column
        tibble::column_to_rownames("common_name")             # Set rownames as common_name
      return(df)
    })
  detections_incidence <- detections_list %>%
    map(~ as.data.frame((.x > 0) * 1))
  
  out.raw <- iNEXT(detections_incidence, q = 0, datatype="incidence_raw", endpoint=100)
  out.raw
  ggiNEXT(out.raw, type = 1) # sample size-based rarefaction/extrapolation curve
  ggiNEXT(out.raw, type = 2) # sample completeness curve
  ggiNEXT(out.raw, type = 3) # coverage-based rarefaction/extrapolation curve
  
}

# Load AVONET species trait metadata
avonet = readxl::read_xlsx(in_path_avonet_traits, sheet = "AVONET2_eBird") %>%
  janitor::clean_names() %>%
  rename(scientific_name = species2, family = family2, order = order2) %>%
  mutate(scientific_name = tolower(scientific_name)) %>%
  filter(scientific_name %in% species_names$scientific_name) %>%
  select(scientific_name, family, order, mass, habitat, habitat_density, migration, trophic_level, trophic_niche, primary_lifestyle)

# Load EltonTraits 1.0 species trait metadata
eltontraits = read_csv(in_path_elton_traits) %>% clean_names() %>% rename(scientific_name = scientific, common_name = english) %>%
  mutate(
    common_name = str_to_lower(common_name), scientific_name = str_to_lower(scientific_name)
  ) %>% mutate(
    common_name = case_when(
      common_name == "american treecreeper" ~ "brown creeper",
      common_name == "winter wren"          ~ "pacific wren",
      common_name == "black-throated grey warbler" ~ "black-throated gray warbler",
      common_name == "common teal" ~ "green-winged teal",
      TRUE ~ common_name   # keep all others unchanged
    )
  ) %>% filter(common_name %in% species_names$common_name) %>%
  select(common_name, diet_inv, diet_5cat, starts_with("for_strat"), body_mass_value)

setdiff(species_names$common_name, eltontraits %>% filter(common_name %in% species_names$common_name) %>% pull(common_name))

# Add both trait data sources to species metadata
species_metadata = left_join(species_names, avonet, by = "scientific_name")
species_metadata = left_join(species_metadata, eltontraits, by = "common_name")

# NOTE: Very few aerial specialist foragers
eltontraits %>% filter(for_strat_aerial >= 10) %>% pull(common_name)

# Get AVONET invertivore community subset
# "Invertivore = species obtaining at least 60% of food resources from invertebrates in terrestrial systems, including insects, worms, arachnids, etc."
invertivores_avonet = species_metadata %>% filter(trophic_niche == "Invertivore") %>% pull(common_name) %>% sort()
# "Aquatic Predator = species obtaining at least 60% of food resources from vertebrate and invertebrate animals in aquatic systems, including fish, crustacea, molluscs, etc."
species_metadata %>% filter(trophic_niche == "Aquatic predator") %>% pull(common_name)

# Get Eltontraits invertivore community subset
# "Percent use of: Invertebrates-general, aquatic invertebrates, shrimp, krill, squid, crustacaeans, molluscs, cephalapod, polychaetes, gastropods, orthoptera, terrestrial Invertebrates, ground insects, insect larvae, worms, orthopterans, flying insects"
species_invert_ETdietGt10 = species_metadata %>% filter(diet_inv >= 10) %>% pull(common_name) %>% sort() # nearly all species
# "Assignment to the dominant among five diet categories based on the summed scores of constituent individual diets."
invertivores_eltontraits = species_metadata %>% filter(diet_5cat == "Invertebrate") %>% pull(common_name) %>% sort()

setdiff(invertivores_eltontraits, invertivores_avonet)
setdiff(invertivores_avonet, invertivores_eltontraits)

# Load guild data
species_guilds = read_csv("data/traits/species_guilds.csv") %>% clean_names() %>% mutate(common_name = tolower(common_name))

# Insectivores
sp_invert = species_metadata %>% filter(diet_inv >= 10) %>% pull(common_name) %>% sort() # most species
sp_invert_primary = species_metadata %>% filter(diet_5cat == "Invertebrate") %>% pull(common_name) %>% sort()

# Foraging guild: Aerial insectivores (e.g. swallows, swifts, flycatchers)
# primarily capture prey while they are flying in the air
sp_g_aerial_invert = species_guilds %>% filter(foraging_guild_cornell %in% c("aerial forager", "flycatching")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

sp_g_aerial_invert_primary = species_guilds %>% filter(foraging_guild_cornell %in% c("aerial forager", "flycatching")) %>%
  filter(common_name %in% sp_invert_primary) %>% pull(common_name)

# Foraging guild: Foliage gleaners (e.g. warblers, vireos)
# typically capture insects located on vegetation or woody substrate
sp_g_foliage_invert = species_guilds %>% filter(foraging_guild_cornell %in% c("foliage gleaner")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Foraging guild: Ground foragers (e.g. american robin)
# procure prey within forest leaf-litter and at the soil surface
sp_g_ground_invert = species_guilds %>% filter(foraging_guild_cornell %in% c("ground forager")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Foraging guild: Bark-probers (e.g. brown creeper, red-breated nuthatch, woodpeckers and sapsuckers)
# extract prey from under bark or by boring into wood
sp_g_bark_invert = species_guilds %>% filter(foraging_guild_cornell %in% c("bark forager")) %>%
  filter(common_name %in% sp_invert) %>% pull(common_name)

# Riparian associates
sp_ripasso = species_guilds %>% filter(rip_asso_rich2002 == "X") %>% pull(common_name) %>% sort()
# Riparian obligats
sp_ripobl  = species_guilds %>% filter(rip_obl_rich2002 == "X")  %>% pull(common_name) %>% sort()

# A priori list of species that:
# - Are primarily insectivorous
# - Are reported to forage on aquatic insects
#
# https://pugetsoundstreambenthos.org/About-BIBI.aspx
# B-IBI is calculated from stream site samples of primarily EPT species...
# - Ephemeroptera, mayflies
# - Plecoptera, stoneflies
# - Trichoptera, caddisflies
#
# ...but also:
# - Diptera, true flies and midges (with aquatic larvae e.g. Chironomidae, Simuliidae, Tipulidae)
# - Bivalvia, freshwater clams (e.g. Unionidae, Sphaeriidae)
# - Gastropoda, freshwater snails (e.g. Lymnaeidae, Planorbidae)
# - aquatic worms
#
# Diets for each species obtained from Birds of the World, unless otherwise commented
sp_apriori = c(
  ## Primary invertebrate diet (Eltontraits Diet-5cat) -----------------------------------------------
  
  # Emergent aquatic invertebrate (EPT) consumers:
  "american dipper",             # Ephemeroptera, Plecoptera, Trichoptera, Diptera (Chironomidae, Simuliidae, Tipulidae)
  "american pipit",              # Ephemeroptera, Plecoptera, Trichoptera, Diptera (Chironomidae, Tipulidae), Odonata
  "bufflehead",                  # Ephemeroptera, Trichoptera, Diptera (Chironomidae, Chaoborinae), Odonata, Hemiptera (Corixidae)
  "cliff swallow",               # Ephemeroptera, Diptera (Chironomidae, Simuliidae), Odonata
  "common yellowthroat",         # Ephemeroptera, Trichoptera, Diptera, https://www.jstor.org/stable/2426510?seq=1
  "golden-crowned kinglet",      # Plecoptera, Trichoptera, Diptera (Chironomus, Tipulidae, Culicidae)
  # "hermit thrush",             # marginal evidence of mayfly predation, https://etd.ohiolink.edu/acprod/odb_etd/ws/send_file/send?accession=osu1236799366&disposition=inline            
  # "house wren",                # marginal evidence of mayfly predation
  "purple martin",               # Ephemeroptera, Trichoptera, Diptera (Chironomidae, Tipulidae), Odonata
  "ruddy duck",                  # Ephemeroptera, Trichoptera, Diptera (Chironomidae, Culicidae), Odonata, Corixidae
  "spotted sandpiper",           # Ephemeroptera, Diptera
  "tree swallow",                # Ephemeroptera, Plecoptera, Trichoptera, Diptera, Bivalvia, Odonata
  "vaux's swift",                # Ephemeroptera, Diptera
  "western wood-pewee",          # Ephemeroptera, Diptera, Odonata
  "wilson's warbler",            # Ephemeroptera, Plecoptera, Trichoptera, Diptera
  "yellow-rumped warbler",       # Trichoptera, Diptera (Tipulidae, Muscidae)

  # Potential aquatic invertebrate consumers (Diptera and/or Odonata):
  # "american kestrel",              # Odonata
  # "barn swallow"                   # Diptera (Tipulidae), Odonata
  # "black-capped chickadee"         # marginal evidence of Diptera (Chironomidae) https://www.jstor.org/stable/2426510?seq=1
  "black-headed grosbeak",           # Unspecified, https://doi.org/10.1002/ecs2.3148
  # "brewer's blackbird",            # Diptera
  # "chestnut-backed chickadee",     # Diptera (potential aquatic larvae Syrphidae, Tabanidae, Tipulidae)
  # "hammond's flycatcher",          # Diptera (Tipulidae)
  # "killdeer",                      # Diptera, Odonata, Hemiptera (Corixidae)
  # "marsh wren",                    # Diptera (Tipulidae), Odonata, Coleoptera (Dytiscidae)
  # "northern pygmy-owl",            # Odonata
  # "northern rough-winged swallow", # Diptera
  # "northern shoveler",             # Diptera (Chironomidae)
  "olive-sided flycatcher",          # Diptera, Odonata, https://doi.org/10.3389/fevo.2021.633160
  # "pacific wren",                  # Diptera (Stratiomyidae)
  # "pacific-slope flycatcher",      # Diptera
  # "pied-billed grebe",             # Odonata, Bivalvia, Gastropoda
  "red-eyed vireo",                  # Diptera, Odonata
  # "ruby-crowned kinglet",          # Diptera
  # "swainson's thrush",             # Diptera (Tipulidae, Bibionidae)
  # "violet-green swallow",          # Diptera
  "warbling vireo",                  # Unspecified, https://doi.org/10.1002/ecs2.3148
  "western tanager",                 # Diptera (Tipulidae), Odonata
  "willow flycatcher",               # Diptera (Tabanidae, Syrphidae), Odonata
  "yellow warbler"                  # Diptera (Chironomidae), https://doi.org/10.1002/ecs2.3148
  
  # Terrestrial or Marine invertebrate consumers (no evidence of freshwater aquatic invert prey):
  
  # "bewick's wren",
  # "black-throated gray warbler",
  # "brown creeper",             # Tricoptera (only one study with sample size of 1, very likely incidental)
  # "bushtit",
  # "cassin's vireo",
  # "downy woodpecker",
  # "hairy woodpecker",
  # "hutton's vireo",
  # "macgillivray's warbler",
  # "marbled murrelet",
  # "northern flicker",
  # "orange-crowned warbler",
  # "pileated woodpecker",
  # "red-breasted sapsucker",
  # "townsend's warbler",
  
  
  ## Invertebrate diet (Eltontraits >= 10%) -----------------------------------------------
  
  # TODO: Emergent aquatic invertebrate (EPT) consumers:

  # ## Primary VertFishScav
  # american bittern
  # barn owl
  # barred owl
  # belted kingfisher         # NOTE: Secondary aquatic invert consumer
  # common merganser
  # common raven
  # glaucous-winged gull
  # great blue heron
  # great horned owl
  # green heron
  # merlin                    # NOTE: Secondary aquatic invert consumer
  # northern harrier
  # red-tailed hawk
  # ring-billed gull

  # ## Primary omnivore
  # american crow
  # american robin
  # california gull
  # fox sparrow
  # hooded merganser
  # lincoln's sparrow
  # mallard
  # red-breasted nuthatch
  # red-winged blackbird
  # savannah sparrow
  # song sparrow              # Ephemeroptera, Diptera (Chironomidae, Tipulidae), Gastropoda https://doi.org/10.1002/ecs2.3148
  # spotted towhee
  # steller's jay
  # virginia rail
  # western meadowlark
  # western screech-owl

  # ## Primary frugivore, granivore
  # anna's hummingbird
  # cedar waxwing
  # rufous hummingbird
  # varied thrush
  # american coot
  # american goldfinch
  # blue-winged teal
  # brown-headed cowbird
  # chipping sparrow
  # cinnamon teal
  # dark-eyed junco
  # eurasian collared-dove
  # evening grosbeak
  # golden-crowned sparrow
  # green-winged teal
  # house finch
  # house sparrow
  # lesser scaup
  # pine siskin
  # purple finch
  # red crossbill
  # ring-necked duck
  # rock pigeon
  # white-crowned sparrow
  # white-throated sparrow
  
) %>% sort()

# Invertivores that are reported to forage on aquatic insects / NOT bark foragers
# Non-bark foragers (possible aquatic insect predators)
sp_aqinv = c(sp_invert_primary[!sp_invert_primary %in% c(sp_g_bark_invert)],
             "belted kingfisher"
) %>% sort()

sp_aerialfoliage = sp_invert_primary[!sp_invert_primary %in% c(sp_g_bark_invert, sp_g_ground_invert)] %>% sort()

# Riparian associate invertivores
sp_ripasso_inv = c(sp_invert[sp_invert %in% c(sp_ripasso)])
sp_ripobl_inv  = c(sp_invert[sp_invert %in% c(sp_ripobl)])

# Exclude certain sites from analysis ------------------------------------------------

# Inspect total detections per site
# total_detections_by_site = detections$long %>% group_by(site_id) %>%
#   summarise(total_detections = sum(n_detections, na.rm = TRUE)) %>%
#   arrange(desc(total_detections))
# ggplot(total_detections_by_site, aes(x = reorder(site_id, total_detections), y = total_detections)) + geom_col()

# # Exclude sites 257 and 259 that are dominated by agriculture
if (exclude_agri_sites) {
  sites_to_exclude = c("257", "259")
  message("Excluding agricultural site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
  site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
  site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
  detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)
}

# Exclude sites with incomplete surveys
sites_to_exclude = setdiff(unique(site_data_reach$site_id), unique(detections$long$site_id))
sites_to_exclude = c(sites_to_exclude, "150", "3097") # ARU at sites 150 and 3097 destroyed mid-survey by water damage
message("Excluding incomplete survey site(s) ", paste(sites_to_exclude, collapse = ", "), " from analysis")
site_data_reach = site_data_reach %>% filter(!site_id %in% sites_to_exclude)
site_data_basin = site_data_basin %>% filter(!site_id %in% sites_to_exclude)
detections$long = detections$long %>% filter(!site_id %in% sites_to_exclude)

message("Retaining ", length(unique(site_data_reach$site_id)), " sites")

# Visualize joint data ----------------------------------------------------------------
message("Visualizing joint data")

presence_absence = detections$long %>% group_by(site_id, common_name) %>%
  summarise(presence = if_else(sum(n_detections, na.rm = TRUE) > 0, 1, 0), .groups = "drop")

# Summarize richness of different groups by site
site_group_richness =  presence_absence %>%
  group_by(site_id) %>%
  summarise(
    rich_all          = sum(presence),
    rich_inv          = sum(presence[common_name %in% sp_invert]),
    rich_inv_primary  = sum(presence[common_name %in% sp_invert_primary]),
    # Foraging guilds
    rich_aerial_inv  = sum(presence[common_name %in% sp_g_aerial_invert]),
    rich_aerial_inv_primary  = sum(presence[common_name %in% sp_g_aerial_invert_primary]),
    rich_foliage_inv = sum(presence[common_name %in% sp_g_foliage_invert]),
    rich_ground_inv  = sum(presence[common_name %in% sp_g_ground_invert]),
    rich_bark_inv    = sum(presence[common_name %in% sp_g_bark_invert]),
    rich_apriori        = sum(presence[common_name %in% sp_apriori]),
    rich_aqinv       = sum(presence[common_name %in% sp_aqinv]),
    rich_aerialfoliage = sum(presence[common_name %in% sp_aerialfoliage]),
    # Riparian habitat association
    rich_ripasso_inv     = sum(presence[common_name %in% sp_ripasso_inv]),
    rich_ripasso     = sum(presence[common_name %in% sp_ripasso]),
    rich_ripobl_inv      = sum(presence[common_name %in% sp_ripobl_inv]),
    rich_ripobl      = sum(presence[common_name %in% sp_ripobl])
  )

# TODO: Rarefied species richness?

## Indicator species analysis
library(indicspecies)
spmat = presence_absence %>% pivot_wider(names_from = common_name, values_from = presence)
d = left_join(spmat, site_data_reach, by = "site_id")
spmat = spmat %>% select(-site_id) %>% as.data.frame()

bibi_excellent = ifelse(d$bibi >= 80, 1, 0)

indval = multipatt(spmat, bibi_excellent, control = how(nperm=999)) 
summary(indval)

# ## NMDS visualization
# library(vegan)
# nmds = metaMDS(presence_absence %>%
#                  pivot_wider(names_from = common_name, values_from = presence, values_fill = 0) %>%
#                  column_to_rownames("site_id"), distance = "bray", k = 3, trymax = 100)
# nmds$stress # stress > 0.2 suggests weak fit
# nmds_scores = as.data.frame(scores(nmds, display = "sites"))
# nmds_scores$site_id = rownames(nmds_scores)
# nmds_plot_data = nmds_scores %>% left_join(site_metadata, by = "site_id")
# species_scores = as.data.frame(scores(nmds, display = "species"))
# species_scores$species = rownames(species_scores)
# 
# ggplot(nmds_plot_data %>% left_join(site_data, by = "site_id"), aes(NMDS1, NMDS2, color = bibi)) +
#   geom_text(data = species_scores, aes(x = NMDS1, y = NMDS2, label = species), color = "gray", size = 3, vjust = -0.5, check_overlap = TRUE) +
#   geom_point(size = 3, alpha = 0.8) +
#   geom_text(aes(label = site_id), vjust = -0.5, size = 3, color = "black") +
#   scale_color_continuous(type = "viridis") +
#   labs(color = "Impervious %")
# 
# ggplot(species_scores, aes(x = 0, y = 0)) +
#   geom_segment(aes(xend = NMDS1, yend = NMDS2),
#                arrow = arrow(length = unit(0.25, "cm")),
#                color = "darkgray", alpha = 0.8) +
#   geom_text(aes(x = NMDS1, y = NMDS2, label = species),
#             color = "black", size = 3, vjust = -0.5, check_overlap = TRUE) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "gray70") +
#   scale_x_continuous(limits = c(-2, 2)) +
#   scale_y_continuous(limits = c(-2, 2)) +
#   labs(title = "Species vectors")

# Richness across AVONET guilds per site
richness_by_trophic_niche = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, trophic_niche) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_primary_lifestyle = presence_absence %>% left_join(species_metadata, by = "common_name") %>%
  group_by(site_id, primary_lifestyle) %>% summarise(count = sum(presence), .groups = "drop")
richness_by_ripasso = presence_absence %>% mutate(ripasso = common_name %in% c(sp_ripasso)) %>% group_by(site_id, ripasso) %>% summarize(count = sum(presence))
richness_by_apriori = presence_absence %>% mutate(apriori = common_name %in% c(sp_apriori)) %>% group_by(site_id, apriori) %>% summarize(count = sum(presence))

ggplot(richness_by_trophic_niche, aes(x = count, y = reorder(site_id, count), fill = trophic_niche)) + geom_col()
ggplot(richness_by_primary_lifestyle, aes(x = count, y = reorder(site_id, count), fill = primary_lifestyle)) + geom_col()
ggplot(richness_by_ripasso, aes(x = count, y = reorder(site_id, count), fill = ripasso)) + geom_col()
ggplot(richness_by_apriori, aes(x = count, y = reorder(site_id, count), fill = apriori)) + geom_col()

# Join with site data
site_data_reach = full_join(site_data_reach, site_group_richness, by = "site_id")
site_data_basin = full_join(site_data_basin, site_group_richness, by = "site_id")

# site_data_reach = left_join(site_data_reach, richness, by = "site_id")
# site_data_reach = left_join(site_data_reach, richness_invert_et, by = "site_id")
# site_data_reach = left_join(site_data_reach, richness_invert_avo, by = "site_id")
# site_data_reach = left_join(site_data_reach, richness_ripasso, by = "site_id")
# site_data_reach = left_join(site_data_reach, richness_apriori, by = "site_id")
# 
# site_data_reach = left_join(site_data_reach, richness_aqinv, by = "site_id")
# 
# site_data_basin = left_join(site_data_basin, richness, by = "site_id")
# site_data_basin = left_join(site_data_basin, richness_invert_et, by = "site_id")
# site_data_basin = left_join(site_data_basin, richness_invert_avo, by = "site_id")
# site_data_basin = left_join(site_data_basin, richness_ripdep, by = "site_id")
# site_data_basin = left_join(site_data_basin, richness_apriori, by = "site_id")

# Richness as a function of different predictors
ggplot(site_data_reach, aes(x = rast_usfs_canopycover_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Canopy cover") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_reach, aes(x = rast_nlcd_impervious_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Reach imperviousness") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_basin, aes(x = rast_nlcd_impervious_sum_proportion, y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Basin imperviousness") + geom_text_repel(aes(label = site_id))

ggplot(site_data_reach, aes(x = bibi, y = rich_inv)) + geom_point() + geom_smooth(method = "lm") +
  labs(title = "BIBI") + geom_text_repel(aes(label = site_id)) 

ggplot(site_data_reach, aes(x = bibi, y = rich_aerial_inv)) + geom_point() + geom_smooth(method = "lm") +
  labs(title = "BIBI") + geom_text_repel(aes(label = site_id))

ggplot(site_data_reach, aes(x = (density_roads_paved), y = rich_inv)) + geom_point() + geom_smooth() +
  labs(title = "Density paved roads") + geom_text_repel(aes(label = site_id)) 

mapview(site_data_reach %>% select(site_id, rich_inv), zcol = "rich_inv")

# Richness of guilds as a function of...
ggplot(left_join(richness_by_trophic_niche, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = trophic_niche, fill = trophic_niche)) +
  geom_point() + geom_smooth(aes(group = trophic_niche), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_primary_lifestyle, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = primary_lifestyle, fill = primary_lifestyle)) +
  geom_point() + geom_smooth(aes(group = primary_lifestyle), method = "lm", se = FALSE)

ggplot(left_join(richness_by_ripasso, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = ripasso, fill = ripasso)) +
  geom_point() + geom_smooth(aes(group = ripasso), method = "lm", se = FALSE)
ggplot(left_join(richness_by_ripasso, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = ripasso, fill = ripasso)) +
  geom_point() + geom_smooth(aes(group = ripasso), method = "lm", se = FALSE)

ggplot(left_join(richness_by_apriori, site_data_reach, by = "site_id"),
       aes(x = bibi, y = count, color = apriori, fill = apriori)) +
  geom_point() + geom_smooth(aes(group = apriori), method = "lm", se = FALSE)
ggplot(left_join(richness_by_apriori, site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = count, color = apriori, fill = apriori)) +
  geom_point() + geom_smooth(aes(group = apriori), method = "lm", se = FALSE)

site_group_richness_long = site_group_richness %>% left_join(site_data_reach) %>% pivot_longer(
  cols = starts_with("rich_"),   # select all richness columns
  names_to = "richness_type",
  values_to = "richness_value"
)
ggplot(site_group_richness_long, aes(x = bibi, y = richness_value, color = richness_type)) +
  geom_point() +
  facet_wrap(~richness_type) +
  geom_smooth(aes(group = richness_type), method = "lm", se = FALSE)

# Species-specific presence/absence as a function of BIBI
ggplot(left_join(presence_absence %>% filter(common_name == "wilson's warbler"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "black-throated gray warbler"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "western wood-pewee"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

ggplot(left_join(presence_absence %>% filter(common_name == "violet-green swallow"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "pacific wren"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "swainson's thrush"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "black-headed grosbeak"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

ggplot(left_join(presence_absence %>% filter(common_name == "golden-crowned kinglet"), site_data_reach, by = "site_id"),
       aes(x = bibi, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE)

# Urban adapted species
ggplot(left_join(presence_absence %>% filter(common_name == "bewick's wren"), site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

ggplot(left_join(presence_absence %>% filter(common_name == "song sparrow"), site_data_reach, by = "site_id"),
       aes(x = rast_nlcd_impervious_sum_proportion, y = presence)) + geom_point() + geom_smooth(method = "glm", method.args = list(family = "binomial"), se = TRUE) + geom_text_repel(aes(label = site_id))

# Structural equation modeling -----------------------------------------------------------------------------

# Calculate pairwise collinearity among predictors
pairwise_collinearity = function(vars, threshold = 0.7) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}

# Multiscale model
{
  # 550 m represents riparian zone within the local reach, and the 90% dispersal distance
  data_reach   = site_data_reach
  # 5 km represents the catchment landscape (roughly basin)
  data_basin = site_data_basin
  d = data.frame(
    "rich_all"         = data_reach$rich_all,
    "rich_inv"         = data_reach$rich_inv,
    "rich_inv_primary" = data_reach$rich_inv_primary,
    "rich_aerial_inv"  = data_reach$rich_aerial_inv,
    "rich_aerial_inv_primary"  = data_reach$rich_aerial_inv_primary,
    "rich_foliage_inv" = data_reach$rich_foliage_inv,
    "rich_ground_inv"  = data_reach$rich_ground_inv,
    "rich_bark_inv"    = data_reach$rich_bark_inv,
    "rich_apriori"        = data_reach$rich_apriori,
    "rich_ripasso"     = data_reach$rich_ripasso,
    "rich_ripobl"      = data_reach$rich_ripobl,
    "rich_ripasso_inv"     = data_reach$rich_ripasso_inv,
    "rich_ripobl_inv"      = data_reach$rich_ripobl_inv,
    "rich_aqinv"       = data_reach$rich_aqinv,
    "rich_aerialfoliage"       = data_reach$rich_aerialfoliage,
    "bibi"          = data_reach$bibi,
    "imp_reach"     = data_reach$rast_nlcd_impervious_sum_proportion,
    "imp_basin"     = data_basin$rast_nlcd_impervious_sum_proportion,
    "fhd_reach"     = data_reach$rast_gedi_fhd_mean,
    "fhd_basin"     = data_basin$rast_gedi_fhd_mean,
    "canopy_reach"  = data_reach$rast_usfs_canopycover_sum_proportion,
    "canopy_basin"  = data_basin$rast_usfs_canopycover_sum_proportion,
    "ed_reach"      = data_reach$edge_density,
    "ed_basin"      = data_basin$edge_density,
    "forest_reach"  = data_reach$nlcd_forest,
    "forest_basin"  = data_basin$nlcd_forest,
    "site_id"       = data_reach$site_id
  )
  pairwise_collinearity(d %>% select(where(is.numeric)))
  
  m_bibi = lm(bibi ~ canopy_reach + imp_basin, d)
  
  # All species
  m_all = glm(rich_all ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_all); plot(sem); print(summary(sem))
  
  # Riparian associates/obligates
  m_ripasso_inv = glm(rich_ripasso_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_ripasso_inv); plot(sem); print(summary(sem))
  
  m_ripobl_inv = glm(rich_ripobl_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_ripobl_inv); plot(sem); print(summary(sem))
  
  # All invertivores (primary and >10% diet)
  m_inv = glm(rich_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_inv); plot(sem); print(summary(sem))
  
  m_inv_primary = glm(rich_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_inv_primary); plot(sem); print(summary(sem))
  
  # Invertivorous foraging guilds
  m_aerial_inv = glm(rich_aerial_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_aerial_inv); plot(sem); print(summary(sem))
  {
    m_aerial_inv_primary = glm(rich_aerial_inv_primary ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(m_bibi, m_aerial_inv_primary); plot(sem); print(summary(sem))
  }
  
  m_foliage_inv = glm(rich_foliage_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_foliage_inv); plot(sem); print(summary(sem))
  
  m_ground_inv = glm(rich_ground_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_ground_inv); plot(sem); print(summary(sem))
  
  m_bark_inv = glm(rich_bark_inv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_bark_inv); plot(sem); print(summary(sem))

  m_apriori = glm(rich_apriori ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_apriori); plot(sem); print(summary(sem))
  
  m_aqinv = glm(rich_aqinv ~ bibi + canopy_reach + imp_reach, d, family = poisson)
  sem = psem(m_bibi, m_aqinv); plot(sem); print(summary(sem))
  
  # Use forest cover in lieu of canopy
  {
    # m_aerialfoliage = glm(rich_aerialfoliage ~ bibi + canopy_reach + imp_reach, d, family = poisson)
    sem = psem(
      lm(bibi ~ forest_reach + imp_basin, d), 
      glm(rich_apriori ~ bibi + forest_reach + imp_reach, d, family = poisson)); plot(sem); print(summary(sem))
  }
  
  # Check overdispersion -- if overdispersed, fit negative binomial
  # simres_pois = simulateResiduals(m_apriori, n = 1000) # plot(simres_pois)
  # testDispersion(simres_pois)
}

print(presence_absence %>% group_by(common_name) %>% summarize(n_sites = sum(presence)) %>% arrange(n_sites) %>% filter(common_name %in% sp_aqinv), n = 100)

species_metadata %>% filter(common_name %in% c(
  "green-winged teal", "mallard", "gadwall", "wood duck", "common merganser", "blue-winged teal", "cackling goose", "caspian tern", "great blue heron", "red-winged blackbird", "american bittern", "green heron"
)) %>% select(common_name, trophic_level, trophic_niche, diet_5cat)
