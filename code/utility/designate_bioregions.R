# Utility/designate_bioregions.R
# 2022-06-14

# Map and raster
library(maps); library(mapdata)
# Read shapefiles of biome
library(sf); library(raster); library(rgdal)
library(dplyr)
library(ggplot2)

# data --------------------------------------------------------------------

qgis_df <- read.csv('data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv')

# Check which region they fall into ---------------------------------------

# Need to read this from directory where all the other files are... *.shp by itself won't work.
wwf_biome <- raster::shapefile("data/WWF_Australia/wwf_terr_ecos.shp")

# Make Spatial Points
rec_points <- SpatialPoints(qgis_df[, c(7, 6)], proj4string = raster::crs(wwf_biome))

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

top_2_biome <- wwf_biome_records %>% 
  dplyr::group_by(species, BIOME) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::arrange(desc(prop)) %>% 
  dplyr::group_by(species) %>% 
  dplyr::slice(1:2)  # Choose top three eco regions for each group

top_2_biome %>% 
  dplyr::distinct(BIOME, .keep_all = FALSE)

ggplot(data = top_3_biome, aes(fill = BIOME, y = prop, x = species)) +
  geom_bar(position = 'stack', stat = 'identity') +
  facet_wrap(~ species)

# Biome designation from WWF_terr_ecos.htm meta data
# 1, Tropical & Subtropical Moist Broadleaf Forests, Tropical Forest
# 2 Tropical & Subtropical Dry Broadleaf Forests, Lesser Sunda Islands
# 3 Tropical & Subtropical Coniferous Forests, 
# 4 Temperate Broadleaf & Mixed Forests
# 5 Temperate Conifer Forests
# 6 Boreal Forests/Taiga
# 7 Tropical & Subtropical Grasslands, Savannas & Shrublands, Tropical grasslands
# 8 Temperate Grasslands, Savannas & Shrublands, Temperate grasslands
# 9 Flooded Grasslands & Savannas
# 10 Montane Grasslands & Shrublands
# 11 Tundra
# 12 Mediterranean Forests, Woodlands & Scrub, Mediterranean
# 13 Deserts & Xeric Shrublands, Arid
# 14 Mangroves

