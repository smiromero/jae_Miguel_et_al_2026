# Create variables

pca_select <- c("PC1", "PC2")
pc.combination <- list(c(1, 2))
plot.titles <- c("PCA")
years.rename <-  setNames(1985:2023, 1:39)
keep.vars <- c('log_sumfrequency', 'log_size', 'log_sev')


if (var.agent == 'wind_barkbeetle') {
  class.levels <- c('Frequent & small' ,
                    'Moderate' ,
                    'Severe & large',
                    'Mild & Rare' ,
                    "Residual")
} else if (var.agent == 'fire') {
  class.levels <- c("Moderate", "Severe & rare", "Large & frequent", "Residual")
}

if (var.agent == 'wind_barkbeetle') {
  var.g <- 4
  var.modelname <- 'VEV'
} else if (var.agent == 'fire') {
  var.g <- 3
  var.modelname <- 'VVE'
}


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

europe.full <- vect(here(
  "data",
  "reference",
  "grids",
  paste0("grids_europe", grid.scale, ".shp")
))
europe.orig <- europe.full %>%
  select(grid_id2, biome) %>%
  mutate(biome7 = as.factor(biome)) %>%
  select(grid_id2, biome7)
europe <- europe.full %>%
  select(-grid_id) %>%
  rename(grid_id = grid_id2) %>%
  aggregate(by = "grid_id") %>%
  select(-agg_n) %>%
  left_join(as.data.frame(europe.orig), by = c("grid_id" = "grid_id2")) %>%
  mutate(biome3 = as.factor(biome7), biome3 = as.factor(
    case_when(
      biome3 == "Tundra" ~ "Boreal",
      biome3 == "boreal" ~ "Boreal",
      biome3 %in% c("broadleafMixed", "coniferous", "grasslands") ~ "Temperate",
      biome3 == "mediterranean" ~ "Mediterranean",
      TRUE ~ biome3
    )
  )) %>%
  distinct(grid_id, .keep_all = TRUE)

europe_df <- europe %>%
  as.data.frame() %>%
  dplyr::select(grid_id, biome3, biome7)


selector.forest <- readRDS(here(
  "data",
  "reference",
  "patch_summary",
  paste0("selector_forest", grid.scale, "_TUM.rds")
)) %>%
  rename('grid_id' = 'grid_id2')


europe.boundaries <- vect(
  here(
    "data/reference/european_boundaries/european_boundaries_study_region_epsg3035.shp"
  )
)
world.boundaries <- vect(
  here(
    "data/reference/world-administrative-boundaries/world-administrative-boundaries_etrs89.shp"
  )
) %>%
  crop(ext(europe.boundaries)) %>%
  filter(continent == 'Europe')
biome.boundaries <- vect(here("data/reference/biomes/european_biomes_epsg3035.shp")) %>%
  aggregate('BIOME_NAME')
soviet.countries <- readRDS(here('data/reference/european_boundaries/soviet_countries.rds'))
grids.noforest <- readRDS(here('data/reference/grids/grids_noforest_TUM.rds'))

europe_wgs84 <- st_transform(st_as_sf(europe), crs = 4326) %>%
  st_make_valid(.)

europe_filter <- europe_wgs84[st_coordinates(st_centroid(europe_wgs84))[, "Y"] <= 66, ]

europe_filter <- st_transform(europe_filter, crs = 3035)
grids_filtered <- setdiff(europe$grid_id, europe_filter$grid_id)
grids_filtered_shp <- europe %>%
  filter(grid_id %in% selector.forest$grid_id)


regime.data <- readRDS(here(
  "data",
  "reference",
  "patch_summary",
  paste0("regime_", n.periods, '_period', grid.scale, "_TUM.rds")
)) %>%
  pull(period.regime.filter) %>%
  .[[1]] %>%
  map(
    ~ rename(.x, grid_id = grid_id2) %>%
      filter(agent == var.agent) %>%
      select(
        "grid_id",
        "agent",
        "log_size",
        "log_vcsize",
        "log_sev",
        "log_vcsev",
        "log_sumfrequency",
        "log_frequency",
        "log_vcfrequency",
        "log_rate",
        "log_sumrate",
        "log_sdrate",
        "log_vcrate",
        "log_sumrate",
        "log_n_years"
      ) %>%
      filter(across(everything(), is.finite)) # %>%
  )

if (n.periods == 1) {
  annual.regime <- readRDS(here(
    "data",
    "reference",
    "patch_summary",
    paste0("regime_", n.periods, '_period', grid.scale, "_TUM.rds")
  ))  %>%
    pull(annual.regime.filter) %>%
    .[[1]] %>%
    rename(grid_id = grid_id2) %>%
    filter(agent == var.agent) %>%
    mutate(period = case_when(between(year, 1985, 2023) ~ '1985-2023'))
} else if (n.periods == 2) {
  annual.regime <- readRDS(here(
    "data",
    "reference",
    "patch_summary",
    paste0("regime_", n.periods, '_period', grid.scale, "_TUM.rds")
  ))  %>%
    pull(annual.regime.filter) %>%
    .[[1]] %>%
    rename(grid_id = grid_id2) %>%
    filter(agent == var.agent) %>%
    mutate(period = case_when(
      between(year, 1985, 2004) ~ '1985-2004' ,
      between(year, 2005, 2023) ~ '2005-2023'
    ))
  
}

logo.ggdraw <- agent_icon(var.agent)
expression_ha_year <- expression(Frequency * "\n" * (No. ~ events ~ per ~
                                                       forest ~ ha ^ {
                                                         -1
                                                       } ~ year ^ {
                                                         -1
                                                       }))

if (var.agent == 'wind_barkbeetle') {
  var.colors <- c('#f9f9f9e5', "#88B680", '#9c89b8', "#F8D67B" , "#CD6463")
  var.colors2 <- c("#CD6463", "#88B680", "#F8D67B", '#f9f9f9e5', '#9c89b8')
  
} else if (var.agent == 'fire') {
  var.colors <- c('#f9f9f9e5', "#F8D67B", "#CD6463", "#88B680")
  var.colors2 <- c("#CD6463", "#F8D67B", '#f9f9f9e5', "#88B680")
  
}

if (var.agent == 'fire') {
  var.colors.boxplot <- c(
    "Western Mediterranean" = "#CD6463",
    'Large & frequent' = "#CD6463",
    "Central Eastern Mediterranean" = "#F8D67B",
    'Moderate' = "#F8D67B",
    "Temperate Boreal" = "#88B680",
    'Severe & rare' = "#88B680",
    "Residual" = "grey95"
  )
  var.colors.boxplot2 <- c(
    "Western Mediterranean" = "#CD6463",
    'Large & frequent' = "#CD6463",
    "Central Eastern Mediterranean" = "#F8D67B",
    'Moderate' = "#F8D67B",
    "Temperate Boreal" = "#88B680",
    'Severe & rare' = "#88B680",
    "Residual" = "#f9f9f9e5"
  )
} else {
  var.colors.boxplot <- c(
    "Central European" = "#CD6463" ,
    'Frequent & small' = "#CD6463" ,
    "Temperate Boreal" = "#F8D67B",
    'Moderate' = "#F8D67B",
    "Atlantic" = '#9c89b8',
    'Severe & large' = '#9c89b8',
    "Mediterranean" = "#88B680" ,
    'Mild & Rare' = "#88B680" ,
    "Residual" = "grey95"
  )
  
  
  var.colors.boxplot2 <- c(
    "Central European" = "#CD6463" ,
    'Frequent & small' = "#CD6463" ,
    "Temperate Boreal" = "#F8D67B",
    'Moderate' = "#F8D67B",
    "Atlantic" = '#9c89b8',
    'Severe & large' = '#9c89b8',
    "Mediterranean" = "#88B680" ,
    'Mild & Rare' = "#88B680" ,
    "Residual" = "#f9f9f9e5"
  )
  
}

biome.colors <- c(
  "Boreal" = "#532E57FF",
  # Azul
  "Temperate" = "#296656FF",
  # Verde
  "Mediterranean" = "#C1B178FF"
)  # Naranja
forest.colors <- c('#532E57FF', '#296656FF', '#C1B178FF')