###############################################################################
# Process and prepare disturbance metrics data across spatial scales 
# Author: Sofía Miguel Romero
###############################################################################

#### Prepare session ####
rm(list = ls(all.names = TRUE)) 
pacman::p_unload(pacman::p_loaded(), character.only = TRUE)

pacman::p_load(tidyverse,
               here,
               terra,
               tidyterra,
               scales,
               janitor,
               conflicted,
               fs)

## Solve package conflicts:
conflicts_prefer(dplyr::filter, dplyr::select, purrr::map)

# Load functions
source(here("lib/00_functions.R"))
list.scales <- c("_25km", "_50km", "_100km")

#### Process disturbance data by grid scale  ####


for (scl in 1:3) {
  grid.scale <- list.scales[scl]
  n.periods <- 3
  
  message("--------------------------------------------------")
  message("Processing scale: ", grid.scale)
  message("--------------------------------------------------")
  
  years.rename <-  setNames(1985:2023, 1:39)
  if (n.periods == 1) {
    decades.list <- list(c(1985, 2023))
    reference.period <- 1
  } else if (n.periods == 2) {
    decades.list <- list(c(1985, 2004), c(2005, 2023))
    reference.period <- 2
    predict.period <- c(1)
  } else if (n.periods == 3) {
    decades.list <- list(c(1985, 1997), c(1998, 2010), c(2011, 2023))
    reference.period <- 3
    predict.period <- c(2, 1)
  }
  
  europe <- vect(here(
    "data",
    "reference",
    "grids",
    paste0("grids_europe", grid.scale, ".shp")
  ))
  
  grid.id <- europe %>%
    as_tibble() %>%
    select(grid_id2, grid_id) %>%
    mutate(year = list(c(1985:2023)), agent = list(as.factor(c(
      "fire", "wind_barkbeetle"
    )))) %>%
    unnest(year) %>%
    unnest(agent)
  
  #### 1.Prepare or load patch-level data  ####

  ##### Option 1: Load and preprocess patch level data #####
  patch.data <- dir_ls(
    here(
      "data",
      "reference",
      "patch_summary",
      paste0("grid_level", grid.scale)
    ),
    recurse = FALSE,
    regexp = "*patch_summary_grid*"
  ) %>%
    map( ~ readRDS(.x)) %>%
    bind_rows()
  
  patch.data <- patch.data %>%
    mutate(
      year = years.rename[year],
      biome = as.factor(biome),
      country = as.factor(str_to_lower(country)),
      agent = as.factor(
        case_when(
          agent == 1  ~ "wind_barkbeetle",
          agent == 2  ~ "fire",
          agent == 3  ~ "harvest",
          agent == 4  ~ "mixed_agents"
        )
      )
    ) %>%
    filter(!agent %in% c("harvest" , "mixed_agents")) %>%
    mutate(severity_mean =  scales::rescale(-severity_mean, to = c(0, 1))) %>%
    select(
      "biome",
      "country",
      "grid_id",
      "year",
      "agent",
      "patch",
      "count",
      "area_ha",
      "severity_mean",
      "severity_sd",
      "grid_ha",
      "forest_ha"
    )
  
  saveRDS(patch.data,
          here(
            "data",
            "reference",
            "patch_summary",
            paste0("all_patch_data", grid.scale, ".rds")
          ))
  
  
  ##### Option 2: Once we have all_data created: #####
  patch.data <- readRDS(here(
    "data",
    "reference",
    "patch_summary",
    paste0("all_patch_data", grid.scale, ".rds")
  )) %>%
    select(-biome, -country)
  
  patch.data.join <- patch.data %>%
    right_join(grid.id, by = c("grid_id", "year", "agent")) %>%
    group_by(grid_id2) %>%
    mutate(
      grid_ha2 = sum(unique(grid_ha), na.rm = TRUE),
      forest_ha2 =  sum(unique(forest_ha), na.rm = TRUE)
    )
  
  #### 2. Create selectors ####

  # Create a selector to filter out grids with small than 10 disturbance
  # years areas of forest
  
  selector.forest <- patch.data.join  %>%
    group_by(grid_id2) %>%
    summarize(
      forest_ha = first(na.omit(forest_ha2)),
      grid_ha = first(na.omit(grid_ha2)),
      percentage_forest = (forest_ha * 100) / grid_ha
    ) %>%
    filter(percentage_forest < 8) %>%
    select(grid_id2) %>%
    rename(grid_id = grid_id2)
  
  # Filter grid with latitudes > 67ºN:
  centroids <- centroids(europe)
  centroids_latlon <- project(centroids, "EPSG:4326")
  coords <- crds(centroids_latlon)
  keep <- coords[, 2] <= 66
  europe.filtered <- europe[keep, ]
  selector.nothern.latitudes <- europe %>%
    mutate(selector = if_else(!grid_id2 %in% europe.filtered$grid_id2, FALSE, NA)) %>%
    filter(selector == FALSE) %>%
    select(grid_id2) %>%
    rename(grid_id = grid_id2) %>%
    as.data.frame()
  
  # Filter grids within grasslands and taiga:
  selector.biomes <- europe %>%
    filter(biome %in% c("grasslands", "Tundra")) %>%
    select(grid_id2) %>%
    rename(grid_id = grid_id2) %>%
    as.data.frame()
  
  selector.all <- selector.forest %>%
    bind_rows(selector.nothern.latitudes) %>%
    bind_rows(selector.biomes) %>%
    distinct() %>%
    rename(grid_id2 = grid_id) %>%
    mutate(selector = FALSE)
  saveRDS(selector.all,
          here(
            "data",
            "reference",
            "patch_summary",
            paste0("selector_forest", grid.scale, ".rds")
          ))
  
  grids.noforest <- europe %>%
    left_join(selector.all, by = "grid_id2") %>%
    filter(selector != FALSE)
  
  saveRDS(grids.noforest,
          here("data/reference/grids/grids_noforest.rds"))
  selector.forest <- selector.all
  
  #### 3. Summarise data at grid-level ####
  
  patch.data.filter <- patch.data.join   %>%
    left_join(selector.forest, by = c("grid_id2")) %>%
    dplyr::filter(is.na(selector)) %>%
    select(-selector)
  
  ##### Summarise yearly #####
  
  annual.regime <- summarise_annual_regime(
    patch.data.filter,
    vars = c("grid_id2", "agent"),
    var_forest = "forest_ha2",
    var_grid = "grid_ha2"
  )
  
  selector.trend <- annual.regime %>%
    mutate(count = if_else(rate != 0, 1, 0)) %>%
    group_by(grid_id2, agent) %>%
    summarize(years_disturbed = sum(count, na.rm = TRUE)) %>%
    filter(years_disturbed <= 10) %>%
    mutate(selector = FALSE) %>%
    select(-years_disturbed)
  
  saveRDS(selector.trend,
          here("data", "reference", "patch_summary", "selector_trend.rds"))
  
  
  annual.regime.filter <- annual.regime # %>%
  # left_join(selector.trend, by = c("grid_id2", "agent")) %>%
  # dplyr::filter(is.na(selector)) %>%
  # select(-selector)
  
  period.regime <- summarise_period(
    annual.regime.filter,
    vars = c("grid_id2", "agent"),
    decades.list =  decades.list
  )
  
  ##### Summarise by period #####
  period.regime.filter <- period.regime %>%
    map2(., 1:n.periods, ~ .x %>%
           ungroup() %>%
           mutate(period = .y)) %>%
    bind_rows() %>%
    group_by(agent) %>%
    mutate(
      across(
        contains("size") |
          contains("sev") |
          contains("frequency") | contains("rate") | contains("n_years") ,
        ~ log(.),
        .names = "log_{.col}"
      )
    ) %>%
    group_by(agent) %>%
    mutate(across(
      where(is.numeric) &
        -contains("sev") &
        !all_of("grid_id2"),
      ~ scales::rescale(., to = c(0, 1))
    )) %>%
    ungroup()
  
  period.regime.filter <- split(period.regime.filter, period.regime.filter$period) %>%
    unname()
  period.regime.filter <- period.regime.filter %>%
    map(. %>%
          select(-period))
  
  #### 4. Save output ####
  grid_regime_tibble <- tibble(
    period.regime.filter = list(period.regime.filter),
    annual.regime.filter = list(annual.regime.filter)
  )
  
  saveRDS(grid_regime_tibble,
          here(
            "data",
            "reference",
            "patch_summary",
            paste0("regime_", n.periods, "_period", grid.scale, ".rds")
          ))
  
}
