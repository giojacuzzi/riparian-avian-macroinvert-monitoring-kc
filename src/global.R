# global.R ------------------------------------------------------------------------------
#
# Global packages, data, and helper functions sourced by multiple scripts 

# Load required packages (automatically install any missing) ----------------------------

if (!exists("pkgs", envir = .GlobalEnv)) {
  message("Loading required packages (automatically installing any missing)")
  pkgs = c(
    # data manipulation
    "arrow",            # cache data compression
    "janitor",          # data cleaning
    "tidyverse",        # general purpose
    "units",            # unit standardization
    # geospatial data
    "geosphere",        # distance metrics
    "leafsync",         # synchronized mapview panels
    "mapview",          # interactive geospatial visualization
    "sf",               # vector data
    "terra",            # raster data
    "tidyterra",        # raster data manipulation
    "tigris",           # political boundaries
    # visualization and plotting
    "ggplot2",
    "ggbeeswarm",       # beeswarm figures
    "ggnewscale",       # multiple scales
    "ggrepel",          # annotations
    "ggeffects",        # marginal effects
    "metR",             # contours
    "patchwork",        # multipanel plots
    # statistics
    "car",              # variance inflation factors
    "DHARMa",           # modeling diagnostics
    "jagsUI",           # hierarchical bayesian MSOM
    "landscapemetrics", # landscape metrics
    "piecewiseSEM",     # structural equation modeling
    "PRROC",            # classifier performance evaluation
    # utility
    "crayon",           # console warnings
    "progress"          # dynamic progress bar
  )
  print(sapply(pkgs, function(pkg) {
    if (!pkg %in% installed.packages()[, "Package"]) install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
    as.character(packageVersion(pkg)) # print package version
  }))
}

# Set ggplot theme -----------------------------------------------------------------------

theme_sleek <- function(base_size = 11, base_family = "") {
  half_line <- base_size/2
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.ticks.length = unit(half_line / 2.2, "pt"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "grey30"),
      strip.text.y = element_text(colour = "grey30"),
      axis.text = element_text(colour = "grey30"),
      axis.title = element_text(colour = "grey30"),
      legend.title = element_text(colour = "grey30", size = rel(0.9)),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 1),
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = rel(0.7), colour = "grey30"),
      legend.key = element_rect(colour = NA, fill = NA),
      legend.background = element_rect(colour = NA, fill = NA),
      plot.title = element_text(colour = "grey30", size = rel(1)),
      plot.subtitle = element_text(colour = "grey30", size = rel(.85))
    )
}

theme_set(theme_sleek())

# Load species trait data ----------------------------------------------------------------

if (!exists("species_names", envir = .GlobalEnv)) {
  species_names = read_lines("data/pam/species_list.txt") %>% as_tibble() %>%
    separate(value, into = c("scientific_name", "common_name"), sep = "_") %>% mutate(
      common_name = tolower(common_name), scientific_name = tolower(scientific_name)
    )
}
  
# Load AVONET species trait data
if (!exists("avonet", envir = .GlobalEnv)) {
  avonet = readxl::read_xlsx("data/traits/AVONET Supplementary dataset 1.xlsx", sheet = "AVONET2_eBird") %>%
    janitor::clean_names() %>%
    rename(scientific_name = species2, family = family2, order = order2) %>%
    mutate(scientific_name = tolower(scientific_name)) %>%
    filter(scientific_name %in% species_names$scientific_name) %>%
    select(scientific_name, family, order, trophic_level, trophic_niche, primary_lifestyle, migration, mass)
}

# Load EltonTraits 1.0 species trait data
if (!exists("eltontraits", envir = .GlobalEnv)) {
  eltontraits = read_csv("data/traits/EltonTraits 1.0/BirdFuncDat.csv", show_col_types = FALSE) %>% clean_names() %>% rename(scientific_name = scientific, common_name = english) %>%
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
}

# setdiff(species_names$common_name, eltontraits %>% filter(common_name %in% species_names$common_name) %>% pull(common_name))

# Load species diet data
if (!exists("diets", envir = .GlobalEnv)) {
  diets = read_csv("data/traits/species_diets.csv", show_col_types = FALSE) %>% clean_names() %>%
    mutate(across(where(is.character), ~ na_if(.x, "")))
}

if (!exists("guilds", envir = .GlobalEnv)) {
  guilds = read_csv("data/traits/species_guilds.csv", show_col_types = FALSE) %>% clean_names() %>%
    mutate(common_name = tolower(common_name)) %>% select(common_name, foraging_guild_cornell, rip_asso_rich2002, rip_obl_rich2002)
}

# Combine all trait data sources
if (!exists("species_traits", envir = .GlobalEnv)) {
  message("Assembling species trait data")
  species_traits = left_join(species_names,  avonet,      by = "scientific_name")
  species_traits = left_join(species_traits, eltontraits, by = "common_name")
  species_traits = left_join(species_traits, diets,       by = c("common_name", "scientific_name"))
  species_traits = left_join(species_traits, guilds,      by = "common_name")

  # Invertivorous species with diets containing macroinvert taxa and foraging behavior either directly in the stream, gleaning from foliage, or in the air (plus documented predation of benthic macroinverts, i.e. dipper, sandpiper, kingfisher), which are all foraging groups with documented predation of aquatic macroinvertebrates.
  # 1) species with documented predation of aquatic macroinvertebrates, as well as primarily invertivorous species with diets containing B-IBI macroinvert taxa and foraging behavior directly in the stream, gleaning from foliage, or in the air.
  # These criteria include species that forage in-stream for aquatic macroinvertebrate larvae as a primary source of nutrition (X), as well as aerial insectivores and foliage gleaners that supplement their diets with emergent macroinvertebrates and have shown positive responses to aquatic subsidies (Y, Z), while excluding primarily terrestrial, bark-probing, and generalist foragers that have shown little evidence for responses to aquatic subsidies (A).
  
  # All grouping (full community):
  species_traits = species_traits %>% mutate(
    group_all = "all"
  )
  
  # Migrant group:
  species_traits = species_traits %>% mutate(
    group_migrant = case_when(
      trophic_niche == "Invertivore" & migration == 3 ~ "migrant",
      TRUE                                            ~ "other"
    )
  )
  
  # Diet group:
  species_traits = species_traits %>% mutate(
    group_diet = case_when(
      # Documented predators of aquatic insect larvae and/or emergent adult forms...
      (common_name %in% c("american dipper", "belted kingfisher", "killdeer", "marsh wren", "pacific wren", "spotted sandpiper")) |
      # ...OR invertivore with diet containing B-IBI taxa
        # https://pugetsoundstreambenthos.org/About-BIBI.aspx
        # B-IBI is calculated from stream site samples of the following insect orders...
        # - Ephemeroptera, mayflies
        # - Plecoptera, stoneflies
        # - Trichoptera, caddisflies
        # - Diptera, true flies and midges (with aquatic larvae e.g. Chironomidae, Simuliidae, Tipulidae)
        # ... and classes of mollusks:
        # - Bivalvia, freshwater clams (e.g. Unionidae, Sphaeriidae)
        # - Gastropoda, freshwater snails (e.g. Lymnaeidae, Planorbidae)
        # - aquatic worms
      trophic_niche == "Invertivore" & benthic_macroinverts == "predator" &
      # ...foraging in stream or within riparian vegetation
      (primary_lifestyle %in% c("Aerial", "Insessorial", "Aquatic")) &
      (!foraging_guild_cornell %in% c("bark forager", "ground forager")) ~ "diet",
      TRUE                                                               ~ "other"
    )
  )
  
  # Foraging group:
  species_traits = species_traits %>% mutate(
    group_forage = case_when(
      trophic_niche == "Invertivore" & (foraging_guild_cornell %in% c("aerial forager", "soaring", "flycatching", "hovering", "aerial dive")) ~ "aerial",
      trophic_niche == "Invertivore" & (foraging_guild_cornell %in% c("foliage gleaner"))                                                     ~ "gleaner",
      trophic_niche == "Invertivore" & (foraging_guild_cornell %in% c("ground forager"))                                                      ~ "ground",
      trophic_niche == "Invertivore" & (foraging_guild_cornell %in% c("bark forager"))                                                        ~ "bark",
      trophic_niche == "Invertivore" & (foraging_guild_cornell %in% c("dabbler", "probing", "stalking", "surface dive"))                      ~ "aquatic",
      TRUE                                                                                                                                    ~ "other"
    )
  )
}

# Geospatial data ----------------------------------------------------------------

crs_standard = "EPSG:32610" # shared coordinate reference system (metric)

# Load site location and survey metadata
site_metadata = read_csv("data/site_metadata.csv", show_col_types = FALSE) %>% clean_names() %>% mutate(site_id = as.character(site_id))

# Create sf points for site locations (reference original crs 4326)
sites_aru = site_metadata %>%
  st_as_sf(coords = c("long_aru", "lat_aru"), crs = 4326)
sites_pssb = site_metadata %>%
  st_as_sf(coords = c("long_pssb", "lat_pssb"), crs = 4326)

# Site exclusions --------------------------------------------------------
sites_to_exclude = c(
  "257", "259",  # Exclude sites 257 and 259 that are dominated by agriculture
  "150", "3097", # ARU at sites 150 and 3097 destroyed mid-survey by water damage
  "155"          # Exclude site 155 to prevent spatial autocorrelation with nearby site 159
)

# Pairwise collinearity ---------------------------------------------
pairwise_collinearity = function(vars, threshold = 0.0) {
  cor_matrix = cor(vars, use = "pairwise.complete.obs", method = "pearson")
  cor_matrix[lower.tri(cor_matrix, diag = TRUE)] = NA
  return(collinearity_candidates = subset(as.data.frame(as.table(cor_matrix)), !is.na(Freq) & abs(Freq) >= threshold))
}
