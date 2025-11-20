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
    mutate(
      benthic_macroinverts = replace_na(benthic_macroinverts, FALSE),
      across(where(is.character), ~ na_if(.x, ""))
    )
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
}

