# Utility/tree_data_match.R
# 2022-06-14

# For BioGeoBEARS analysis, we need a species-level time-calibrated phylogeny. 
# Will only include Australian species of Anilios as we don't have data on A. erycinus 

# libraries ---------------------------------------------------------------
library(ape); library(phytools)

# Read tree ---------------------------------------------------------------

# Time calibrated tree from SqCL probe and mcmctree analysis with fossils
# Full tree with fossil calibrations
fos_tree_full <- ape::read.nexus("data/tree/20220420mcmctree.tre")

# tips we want for blindsnakes
blindsnakes_samples <- c('Anilios_affinis_J93635','Anilios_ammodytes_R158097','Anilios_aspina_J91822','Anilios_australis_R115859',
                         'Anilios_bicolor_R55342','Anilios_bituberculatus_R63990','Anilios_broomi_J87805','Anilios_centralis_R14631',
                         'Anilios_diversus_R112027','Anilios_endoterus_R102725','Anilios_ganei_R165000', 
                         'Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272',
                         'Anilios_guentheri_R108431','Anilios_hamatus_R136277','Anilios_howi_R146381','Anilios_kimberleyensis_R164213', 
                         'Anilios_leptosoma_R119241','Anilios_leucoproctus_J87547', 
                         'Anilios_ligatus_R31019','Anilios_longissimus_R120049','Anilios_margaretae_R163269',
                         'Anilios_nigrescens_R31022','Anilios_obtusifrons_R146400','Anilios_pilbarensis_R108813',
                         'Anilios_pinguis_R146995','Anilios_polygrammicus_R98715','Anilios_proximus_R132486', 
                         'Anilios_silvia_J46128','Anilios_splendidus_R119900','Anilios_systenos_R114894', 
                         'Anilios_torresianus_J47640','Anilios_tovelli_R21063','Anilios_troglodytes_R146048', 
                         'Anilios_unguirostris_R21669','Anilios_waitii_R113302','Anilios_wiedii_J59020', 
                         'Anilios_yirrikalae_J85695','Ramphotyphlops_cfwaitii_R51244', 
                         'Anilios_ligatus_R19109'
)

# Subset the tree
fos_tree <- ape::keep.tip(phy = fos_tree_full, tip = blindsnakes_samples)

plotTree(fos_tree)
# Rename tips (for grypus we have three potential species)
fos_tree$tip.label[which(fos_tree$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272', 'Anilios_ligatus_R31019'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_grypusET_R55272', 'Anilios_ligatusE_R31019')

# Just keep the species name (without rego)
fos_tree$tip.label <- stringr::str_replace(string = fos_tree$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")







###############################
##### MITOCHONDRIAL TREE ######
###############################

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

# write.nexus(subset_mt_tree, file = 'data/tree/tree.nwk')
