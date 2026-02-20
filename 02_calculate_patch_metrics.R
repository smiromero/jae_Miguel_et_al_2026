###############################################################################
# This script processes the 'disturbance_agent' raster layers to calculate
# disturbance metrics
# Author: Sofía Miguel Romero
###############################################################################

## Load packages
library(tidyverse)
library(dplyr)
library(tidyterra)
library(terra)
library(here)
library(sf)
library(conflicted)
library(fs)
library(landscapemetrics)
library(mapview)
library(skimr)

## Load the shapefile containing Europe's grids
grids_shp <- vect( here("data", "reference", "grids", "grids_europe.shp"))

## Define local data directory and the path to disturbance layers
dist_dir <- paste0(here(),'/data/reference/Viana2024_forest_disturbance_atlas/europe/')

## define grid scale: 
list.scales <- c("_25km", "_50km", "_100km")
grid.scale <- list.scales[2]

## Preprocess disturbance agent layers

# Create a tibble with biome and country names
country_biomes_names <- tibble(
  biome_country = paste(  grids_shp$biome,  str_to_lower(grids_shp$country),   sep = '_'),
  country =  str_to_lower(grids_shp$country))

## Loop through all countries and biomes within the grids
all_grids <- map(1:nrow(country_biomes_names), function(country){
  
    # Print information about the current country for progress tracking
    print(paste('country:', country_biomes_names$country[[country]], 
                '- Number:', country, 'of', nrow(country_biomes_names)))
    country_name <- country_biomes_names$country[[country]]
    
    # Extract the grids for the current country
    country_grids <- grids_shp[country]
    
    # Load the disturbance agent raster layers and crop them to match the grids
    agent_grid <-  crop(rast(paste0(dist_dir,country_name,
                                    '/disturbance_agent_1985_2023_',
                                    country_name ,'.tif')),
                        country_grids, mask = TRUE)  %>% 
      clamp(., lower= 1, upper=4, values=FALSE)  #  by setting values=FALSE then values above/below are set to NA.
    
    # Load the disturbance severity layer and crop it
    severity_grid <- crop(rast(paste0(dist_dir,country_name,'/disturbance_severity_1985_2023_', country_name,'.tif')),
                            country_grids, mask = TRUE) 
    # Load the forest cover mask, calculate areas, and rename the column with forest area in hectares
    forest_grid<- crop(rast(paste0(dist_dir,country_name, '/forest_mask_', country_name,'.tif')),
                       country_grids, mask = TRUE) %>% 
      lsm_l_ta(.) %>% 
      select(value) %>%
      rename("forest_ha" = "value") 
   
    # Calculate the grid area in hectares
    grid_area <- expanse(country_grids, unit = 'ha')
   
    ## Summarize the data by year and calculate metrics for each disturbance patch
    patch_metrics_grid_year <- map(1:nlyr(agent_grid), function(yr) {
      return(
        tibble(
          'patch' = values(patches(agent_grid[[yr]]))[,1], # Use ladscapemetrics::get_patches()
          # or terra::patches()
          'agent' =  values(agent_grid[[yr]])[,1],  # Extract disturbance agent values
          'severity' = values(severity_grid[[yr]])[,1]) %>%  # Extract severity values
          dplyr::filter(!if_all(everything(), is.na)) %>%
          group_by(agent, patch) %>% 
          reframe(   
            biome = country_grids$biome, 
            country = country_grids$country,  
            grid_id = country_grids$grid_id,
            grid_ha = grid_area,
            forest_ha = forest_grid$forest_ha,
            year = yr,
            count = n(),
            area_ha = n()*0.09,
            severity_mean = mean(severity, na.rm = TRUE),
            severity_sd = sd(severity, na.rm = TRUE),
            ) %>% 
          ungroup())
        }) %>%  
        bind_rows() %>% 
        na.omit(patch)
        # Save patch metrics for the grid to an RDS file
    saveRDS(
      patch_metrics_grid_year,
      here( "data", "reference", "patch_summary",
            paste0("grid_level", grid.scale),
      paste0("patch_summary_grid_", country_grids$grid_id, ".rds"))
         )
    return(patch_metrics_grid_year)
    
  }) %>% 
  
  bind_rows()


# saveRDS(
#   all_grids,
#   here(
#     "data",
#     "reference",
#     "patch_summary",
#     paste0("all_grids_summary", country_grids$grid_id, ".rds")
#   )
# )

  
  
  