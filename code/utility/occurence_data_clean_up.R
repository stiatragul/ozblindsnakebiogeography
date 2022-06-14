# Utility/occurence_data_clean_up.R
# 2021-10-08

# Rename columns and filter data from Atlas of the Living Australia

library(dplyr); library(tidyr)
library(janitor); library(stringr)
# check sampling bias (require install via github as of 20220614 so need devtools)
library(sampbias)


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
  dplyr::filter(stringr::str_detect(string = genus, pattern = "Anilios|Ramphotyphlops|Sundatyphlops")) %>% 
  dplyr::distinct(catalogue_number, .keep_all = TRUE) %>% # remove duplicate catalogue number
  dplyr::arrange(species)  

# Write to csv and inspect the data when mapped to QGIS. 
# write.csv(qgis_df, file = 'data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv', row.names = FALSE)
# write.csv(qgis_df, file = 'data/WWF_Australia/2021_blindsnake_distribution.csv', row.names = FALSE)

# After inspecting data, I decided to omit the following samples

## Omit all Anilios australis from Museums Victoria because probably mis-id as A. australis is from WA. MV records are from eastern Australia (longitude 135).
mv_australis <- qgis_df %>% 
  dplyr::filter(institution != "") %>%      # omit records without institution
  dplyr::filter(species == "Anilios australis" & decimal_longitude_wgs84 > 135) 



# omit records without institution
no_museum <- qgis_df %>% 
  dplyr::filter(institution == "")

# A. affinis records that are 
affinis_rec <-  c("R05278", "R804", "R22664")

affinis_no <- qgis_df %>%  
  dplyr::filter(catalogue_number %in% affinis_rec)

delete_rec <- rbind(mv_australis, no_museum, affinis_no)

qgis_df_refine <- qgis_df %>% dplyr::anti_join(delete_rec)

## SUPPLEMENTARY TABLE FOR OCCURENCE

# Write to csv and inspect the data when mapped to QGIS. 
# write.csv(qgis_df_refine, file = 'data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv', row.names = FALSE)
# write.csv(qgis_df_refine, file = 'data/WWF_Australia/2021_blindsnake_distribution.csv', row.names = FALSE)


# Plot distribution to check ----------------------------------------------

distro.map <- leaflet(qgis_df_refine) %>%
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

distro.maps