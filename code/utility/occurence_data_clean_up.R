# Utility/occurence_data_clean_up.R
# 2021-10-08

# Rename columns and filter data from Atlas of the Living Australia

library(dplyr)
library(tidyr)
library(janitor)

# Read data ---------------------------------------------------------------

df <- read.csv("data/2021_ALA_blindsnake_occurence_data/records-2021-10-08.csv", header = TRUE) |> 
  clean_names()

df

df |> tibble() |> 
  colnames()


# Columns to subset
subset_cols <- c("genus", "species", "scientific_name_intepreted", "catalogue_number", "institution", 
                 "decimal_latitude_wgs84", "decimal_longitude_wgs84", 
                 "year", "month", "day", "record_id")

qgis_df <- df |>
  select(subset_cols) |> 
  filter(!is.na(decimal_latitude_wgs84)) |> 
  filter(species != "",
         catalogue_number != "") |> 
  arrange(species) 

# Write to csv
# write.csv(qgis_df, file = 'data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv')


# Map
qgis_df <- qgis_df |> 
  tibble() |> 
  filter(species != "",
         catalogue_number != "") |> 
  arrange(species) 

# Read shapefiles of biome
library(sf)
library(raster)

# Need to read this from directory where all the other files are... *.shp by itself won't work.
wwf_biome <- raster::shapefile("C:/Users/ST/Documents/GIS DataBase/WWF_Australia/wwf_terr_ecos.shp")

# Make Spatial Points
rec_points <- SpatialPoints(qgis_df[, c(7, 6)], proj4string = crs(wwf_biome))

wwf_biome_records <- extract(wwf_biome, rec_points, method='bilinear')


length(wwf_biome_records)
length(qgis_df$species)

# Join species names to each record
wwf_biome_records$species <- qgis_df$species

### Number of occurence in each eco region.
wwf_biome_records |> 
  group_by(species, ECO_NAME) |> 
  summarise(number = n())


