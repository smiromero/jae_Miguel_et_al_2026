###############################################################################
# This script generates the spatial grid used to extract and aggregate 
# disturbance data across the study area.
# Author: Sofía Miguel Romero
###############################################################################

#
library(tidyverse)
library(tidyterra)
library(terra)
library(here)
library(sf)
library(conflicted)
library(fs)
library(grid)
library(gridExtra)

#
conflicts_prefer(
  dplyr::filter,
  dplyr::select
)

# create grids for each biome and country 
create_grid<- function(grid_size,  base_map){
  
  grid_result <- st_make_grid(base_map, c(grid_size, grid_size), what = "polygons", square = FALSE) %>% 
    st_sf(.) %>%
    mutate(grid_id = 1:length(lengths(.)),
           grid01 = lengths(st_intersects(. , st_as_sf(base_map)))) %>% 
    tidyterra::filter(., grid01 > 0) %>% 
    mutate(grid_id = 1:nrow(.)) %>% 
    select(grid_id)
  
  return( terra::crop(vect(grid_result), base_map))
  
}

# files path outside the project folder
data_cmemory <- "" # write here the path to your project folder

biomes <- dir_ls(paste( data_cmemory,
                             "reference", "biomes", "biomes_divided",  sep = "/"),
                      regexp = '*.shp',
                      recurse = FALSE) %>% 
  map(~vect(.x) %>% 
        select(BIOME_NAME) %>% 
        rename('biome' = 'BIOME_NAME'))
#
countries <- vect(paste(data_cmemory, "reference", 
                         "european_boundaries", 
                        "european_boundaries_study_region_epsg3035.shp", sep = "/")) %>% 
  select(name) %>% 
  rename( 'country' = 'name')

europe <- aggregate(countries) 

europe_grid <- create_grid(50000, europe)      # 25 km, 50km and 100 km grid scales 

europe_grid <- biomes %>% 
  map( \(x) terra::intersect(europe_grid, x) %>% 
         terra::intersect(., countries))  %>% 
  vect() %>% 
  mutate(biome =  case_when( 
    biome == "Mediterranean Forests, Woodlands & Scrub"  ~ "mediterranean",
    biome == "Temperate Broadleaf & Mixed Forests"  ~ "broadleafMixed",
    biome == "Temperate Conifer Forests"    ~ "coniferous",
    biome == "Tundra"    ~ "Tundra",
    biome == "Boreal Forests/Taiga"   ~ "boreal",
    biome == "Temperate Grasslands, Savannas & Shrublands"  ~ "grasslands")) %>% 
  arrange(grid_id)
#

europe_grid <- europe %>% 
  distinct(geometry, .keep_all = TRUE) %>% 
  rename('grid_id2' = 'grid_id') %>% 
  mutate(grid_id = row_number())
#
writeVector(europe_grid, here("data", "reference", "grids", "grids_europe_50km.shp"), overwrite = TRUE)
