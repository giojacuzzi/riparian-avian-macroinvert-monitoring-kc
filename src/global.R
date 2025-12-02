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
    "ggnewscale",       # multiple scales
    "ggrepel",          # annotations
    "ggeffects",        # marginal effects
    "patchwork",        # multipanel plots
    # statistics
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

theme_set(theme_minimal())

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
  species_traits = species_traits %>% mutate(
    invert_predator = ifelse(
      benthic_macroinverts == "predator" & diet_5cat == "Invertebrate",
      "invert_predator", "NA"
    ))
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
