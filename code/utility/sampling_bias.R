# Utility/sampling_bias.R
# 2022-06-14

# Map and raster
library(maps); library(maptools)
# Read shapefiles of biome
library(sf); library(raster); library(rgdal)
library(sampbias)


# data --------------------------------------------------------------------

qgis_df <- read.csv('data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv')

# Check which region they fall into ---------------------------------------

# Map
qgis_df <- qgis_df %>% 
  tibble() %>% 
  filter(species != "",
         catalogue_number != "") %>% 
  arrange(species) 

# Sample bias
## Prepare
samptest <- qgis_df[, c('species', 'decimal_longitude_wgs84', 'decimal_latitude_wgs84')]
names(samptest) <- c('species', 'decimalLongitude', 'decimalLatitude')

# samptest$species <- gsub(pattern = ' ', replacement = '_', x = samptest$species)

# Export for Infomap Bioregions 
write.csv(x = samptest, file = 'data/20220218_occurence_ALA_cleaned.csv', row.names = FALSE)

# Running sampbias
samp.out <- sampbias::calculate_bias(x = samptest, res = 1, terrestrial = TRUE)

summary(samp.out)
plot(samp.out)

# Project bias on map
proj <- project_bias(samp.out)
map_bias(proj, type = "log_sampling_rate")
