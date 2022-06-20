# Utility/designate_bioregions.R
# 2022-06-14

# Map and raster
library(maps); library(mapdata)
# Read shapefiles of biome
library(sf); library(raster); library(rgdal)
library(dplyr); library(tidyr)
library(ggplot2)

# data --------------------------------------------------------------------

qgis_df <- read.csv('data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv')
alias <- read.csv('data/WWF_Australia/wwf_biome_name_alias.csv')
alias$BIOME <- as.character(alias$BIOME)

# Check which region they fall into ---------------------------------------

# Need to read this from directory where all the other files are... *.shp by itself won't work.
wwf_biome <- raster::shapefile("data/WWF_Australia/wwf_terr_ecos.shp")

# Make Spatial Points
rec_points <- sp::SpatialPoints(qgis_df[, c(7, 6)], proj4string = raster::crs(wwf_biome))

wwf_biome_records <- raster::extract(wwf_biome, rec_points, method='bilinear')

length(wwf_biome_records)
length(qgis_df$species)

# Join species names to each record
wwf_biome_records$species <- qgis_df$species
wwf_biome_records$BIOME <- as.character(wwf_biome_records$BIOME)

### Number of occurence in each eco region.
top_3_ecoregions <- wwf_biome_records %>% 
  dplyr::group_by(species, ECO_NAME) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::arrange(desc(prop)) %>% 
  dplyr::group_by(species) %>% 
  dplyr::slice(1:3)  # Choose top three eco regions for each group

ggplot(data = top_3_ecoregions, aes(fill = ECO_NAME, y = prop, x = species)) +
  geom_bar(position = 'stack', stat = 'identity')

sp_biomes <- wwf_biome_records %>% 
  dplyr::group_by(species, BIOME) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::arrange(desc(prop)) %>% 
  dplyr::group_by(species) %>% 
  # dplyr::slice(1:2) %>% # Choose top three eco regions for each group 
  dplyr::left_join(alias, by = 'BIOME')

# Pivot long for BioGeoBears format
output_sp_biomes <- sp_biomes %>% 
  select(-BIOME, -prop, -full_name) %>% 
  tidyr::pivot_wider(names_from = short_name, values_from = n, values_fill = 0) %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(species = stringr::str_replace(species, pattern = ' ', replacement = '_'))

write.csv(output_sp_biomes, file = 'data/intermediate_data/geo_file_precursor.csv', quote = FALSE, row.names = FALSE)
