# This script generates a tibble containing data for all biomes. Additionally, for each biome, 
# it includes shapefiles of the countries that intersect with the biomes and their corresponding grids.
# Note: This script is used by other R scripts (e.g., 01_5_biome_country_grid_create_plots.r),
# so any changes made here will impact those scripts as well.

# Load packages
library(tidyverse)
library(tidyterra)
library(terra)
library(here)
library(sf)
library(conflicted)
library(fs)
library(grid)
library(gridExtra)

# create grids for each biome and country 
create_grid<- function(grid_size,  base_map){
  
  grid_result <- st_make_grid(base_map, c(grid_size, grid_size), what = "polygons", square = FALSE) %>% 
    #To sf and add grid ID
    st_sf(.) %>%
    # add grid ID
    mutate(grid_id = 1:length(lengths(.)),
           # Intersect both layers. Use st_as_sf if the base map is a SpatVect
           grid01 = lengths(st_intersects(. , st_as_sf(base_map)))) %>% 
    # Remove grids falling outside the basemap boundaries
    tidyterra::filter(., grid01 > 0) %>% 
    # Recalculate the IDs 
    mutate(grid_id = 1:nrow(.)) %>% 
    select(grid_id)
  
  return( terra::crop(vect(grid_result), base_map))
  
}

# files path outside the project folder
data_cmemory <- "C:/Users/SOFIA/Documents/novel_disturbances/data_cmemory"

biomes <- dir_ls(paste( data_cmemory,
                             "reference", "biomes", "biomes_divided",  sep = "/"),
                      regexp = '*.shp',
                      recurse = FALSE) %>% 
  map(~vect(.x) %>% 
        select(BIOME_NAME) %>% 
        rename('biome' = 'BIOME_NAME'))

countries <- vect(paste(data_cmemory, "reference", 
                         "european_boundaries", "european_boundaries_study_region_epsg3035.shp", sep = "/")) %>% 
  select(name) %>% 
  rename( 'country' = 'name')

europe <- aggregate(countries) # Unite all polygons(aka country boundaries) of the shapefile

europe_grid <- create_grid(25000, europe)       

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


europe_grid <- europe %>% 
  distinct(geometry, .keep_all = TRUE) %>% 
  rename('grid_id2' = 'grid_id') %>% 
  mutate(grid_id = row_number())

writeVector(europe_grid, here("data", "reference", "grids", "grids_europe_25km.shp"), overwrite = TRUE)
