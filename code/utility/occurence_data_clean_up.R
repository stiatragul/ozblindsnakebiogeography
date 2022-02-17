# Utility/occurence_data_clean_up.R
# 2021-10-08

# Rename columns and filter data from Atlas of the Living Australia

library(dplyr)
library(tidyr)
library(janitor)
# Read shapefiles of biome
library(sf)
library(raster)
# check sampling bias
library(sampbias)
library(ape)
library(ggplot2)

# Read data ---------------------------------------------------------------

df <- read.csv("data/2021_ALA_blindsnake_occurence_data/records-2021-10-08.csv", header = TRUE) |> 
  clean_names()


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
  distinct(catalogue_number, .keep_all = TRUE) |> # remove duplicate catalogue number
  arrange(species)  

qgis_df 

## SUPPLEMENTARY TABLE FOR OCCURENCE

# Write to csv
# write.csv(qgis_df, file = 'data/2021_ALA_blindsnake_occurence_data/2021_blindsnake_distribution.csv', row.names = FALSE)



# Check which region they fall into ---------------------------------------

# Map
qgis_df <- qgis_df |> 
  tibble() |> 
  filter(species != "",
         catalogue_number != "") |> 
  arrange(species) 

# Sample bias
## Prepare
samptest <- qgis_df[, c('species', 'decimal_longitude_wgs84', 'decimal_latitude_wgs84')]
names(samptest) <- c('species', 'decimalLongitude', 'decimalLatitude')

# samptest$species <- gsub(pattern = ' ', replacement = '_', x = samptest$species)

# Export for Infomap Bioregions 
write.csv(x = samptest, file = 'data/20220218_occurence_ALA_cleaned.csv', row.names = FALSE)

# Running sampbias
# samp.out <- sampbias::calculate_bias(x = samptest, res = 1)

# summary(samp.out)
# plot(samp.out)

# Need to read this from directory where all the other files are... *.shp by itself won't work.
wwf_biome <- raster::shapefile("data/WWF_Australia/wwf_terr_ecos.shp")

# Make Spatial Points
rec_points <- SpatialPoints(qgis_df[, c(7, 6)], proj4string = crs(wwf_biome))

wwf_biome_records <- extract(wwf_biome, rec_points, method='bilinear')

length(wwf_biome_records)
length(qgis_df$species)

# Join species names to each record
wwf_biome_records$species <- qgis_df$species

### Number of occurence in each eco region.
top_3_ecoregions <- wwf_biome_records |> 
  dplyr::group_by(species, ECO_NAME) |> 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(prop = n / sum(n)) %>% 
  dplyr::arrange(desc(prop)) %>% 
  dplyr::group_by(species) %>% 
  dplyr::slice(1:3)  # Choose top three eco regions for each group
  
  
ggplot(data = top_3_ecoregions, aes(fill = ECO_NAME, y = prop, x = species)) +
  geom_bar(position = 'stack', stat = 'identity')
  





# tree --------------------------------------------------------------------

mt_tree <- ape::read.nexus('data/tree/MCC_10_meanh.tre')

# MANIPULATE TREE ---------------------------------------------------------
# Keep only tips of species

tip_keep <- c(
  'acuticauduso', 'suboculariso', 'affinis1',
  'wiedii92', 'ganei36', 'ligatus70',
  'kimberleyensis65', 'troglodytes85', 'polygrammicus81',
  'nigrescens75', 'silvia84', 'guentheri51',
  'howi62', 'unguirostris87', 'grypus43', 'leptosoma67',
  'leptosoma68', 'longissimus73', 'bicolor14',
  'pinguis80', 'bituberculatus19', 'proximus83',
  'australis13', 'endoterus35', 'hamatus60',
  'pilbarensis79', 'centralis21', 'waitii90',
  'ammodytes11', 'diversus26')

# sub tree gets rid of number and SH
sub_mt_tree <- mt_tree %>%
  keep.tip(., tip_keep)

sub_mt_tree$tip.label[which(sub_mt_tree$tip.label == 'leptosoma67')] <- "systenos24"

# Force ultrametric 
sub_mt_tree <- phytools::force.ultrametric(sub_mt_tree,"extend")

# Clean up tip label
sub_mt_tree$tip.label <- gsub("[0-9]+", '', sub_mt_tree$tip.label) # rename so it so we have tips labels without numbers
sub_mt_tree$tip.label <- gsub("^", 'Anilios ', sub_mt_tree$tip.label) # rename so it so we have tips labels without numbers


## Rename the outgroups
sub_mt_tree$tip.label[which(sub_mt_tree$tip.label == 'Anilios acuticauduso')] <- 'Ramphotyphlops acuticaudus'
sub_mt_tree$tip.label[which(sub_mt_tree$tip.label == 'Anilios suboculariso')] <- 'Acutotyphlops subocularis'


# FILTER DATA TO MATCH TREE TIPS -------------------------------------------

tree_tips_df <- as.data.frame(sub_mt_tree$tip.label)
names(tree_tips_df) <- 'species'

# tree_tips_df$species_in_tree <- sub('Anilios ', tree_tips_df$species, replacement = "", perl = TRUE)

# subset only species in tree

subset_samp <- samptest %>% 
  dplyr::semi_join(tree_tips_df, by = 'species') %>% 
  group_by(species) %>% 
  count()


# prune tree again for only tips that we have data for
shape_data_tips <- as.vector()

# PRUNED tree
subset_mt_tree <- sub_mt_tree %>%
  keep.tip(., subset_samp$species)

write.nexus(subset_mt_tree, file = 'data/tree/tree.nwk')
