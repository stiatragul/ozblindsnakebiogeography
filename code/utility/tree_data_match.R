# Utility/tree_data_match.R
# 2022-06-14

library(ape); library(phytools)

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
