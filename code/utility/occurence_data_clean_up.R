# Utility/occurence_data_clean_up.R
# 2021-10-08

# Rename columns and filter data from Atlas of the Living Australia

library(dplyr); library(tidyr)
library(janitor)
# check sampling bias (require install via github as of 20220614 so need devtools)
library(sampbias)
library(ape); library(phytools)
library(ggplot2)


# Read data ---------------------------------------------------------------

df <- read.csv("data/2021_ALA_blindsnake_occurence_data/records-2021-10-08.csv", header = TRUE) %>% 
  clean_names()

df %>% tibble() %>% 
  colnames()

# Columns to subset
subset_cols <- c("genus", "species", "scientific_name_intepreted", "catalogue_number", "institution", 
                 "decimal_latitude_wgs84", "decimal_longitude_wgs84", 
                 "year", "month", "day", "record_id")

qgis_df <- df %>%
  dplyr::select(all_of(subset_cols)) %>% 
  dplyr::filter(!is.na(decimal_latitude_wgs84)) %>% 
  dplyr::filter(species != "",
         catalogue_number != "") %>% 
  dplyr::distinct(catalogue_number, .keep_all = TRUE) %>% # remove duplicate catalogue number
  dplyr::arrange(species)  

## SUPPLEMENTARY TABLE FOR OCCURENCE

# Write to csv
# write.csv(qgis_df, file = 'data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv', row.names = FALSE)


# Plot distribution to check ----------------------------------------------

distro.map <- leaflet(qgis_df) %>%
  addProviderTiles(providers$OpenTopoMap) %>%
  # addProviderTiles(providers$Esri.WorldPhysical) %>%
  addCircleMarkers(lng = ~decimal_longitude_wgs84, lat = ~decimal_latitude_wgs84,
                   radius = 0.5,
                   popup = ~glue("<h3>Species: {species}", "<br>",
                                 label = ~glue("{species}"),
                                 labelOptions = labelOptions(noHide = F, textsize = "15px"),
                                 stroke = FALSE, fillOpacity = 0.5)) %>% 
  setView(lng = 133, lat = -25, zoom = 4.4) %>% 
  addMiniMap(position = 'bottomleft') %>% 
  addScaleBar(position = "bottomright")

distro.map

# After this step, need to check locality records manually to make sure everything makes sense. 
# There may be some odd locality or ones that fall in the water. 
