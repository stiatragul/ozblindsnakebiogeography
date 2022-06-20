# Utility/tree_data_match.R
# 2022-06-14

# For BioGeoBEARS analysis, we need a species-level time-calibrated phylogeny. 
# Will only include Australian species of Anilios as we don't have data on A. erycinus 

### NOTE ####
# DEC, DEC+J, and most methods in use in historical biogeography implicitly assume that the observed tree is the complete tree, 
# and that lineage-level speciation and extinction are independent of the range evolution process*

# libraries ---------------------------------------------------------------
library(ape); library(phytools)
library(stringr)

# Read tree ---------------------------------------------------------------

# Time calibrated tree from SqCL probe and mcmctree analysis with fossils
# Full tree with fossil calibrations
fos_tree_full <- ape::read.nexus("data/tree/20220420mcmctree.tre")

fos_tree_full <- phytools::force.ultrametric(fos_tree_full,"extend")
phytools::plotTree(fos_tree_full)

fos_tree_full$tip.label

# tips we want for blindsnakes (the best we can do since we do not have all described species)
blindsnakes_samples <- c('Anilios_affinis_J93635','Anilios_ammodytes_R158097','Anilios_aspina_J91822','Anilios_australis_R115859',
                         'Anilios_bicolor_R55342','Anilios_bituberculatus_R63990','Anilios_broomi_J87805','Anilios_centralis_R14631',
                         'Anilios_diversus_R112027','Anilios_endoterus_R102725','Anilios_ganei_R165000', 
                         'Anilios_grypus_R108596',   # Western Australia (W)
                         'Anilios_grypus_R157297',   # Sister to leptosoma-complex (NW)
                         'Anilios_grypus_R55272',    # Eastern (ET)
                         'Anilios_guentheri_R108431','Anilios_hamatus_R136277','Anilios_howi_R146381','Anilios_kimberleyensis_R164213', 
                         'Anilios_leptosoma_R119241','Anilios_leucoproctus_J87547', 
                         'Anilios_ligatus_R31019','Anilios_longissimus_R120049','Anilios_margaretae_R163269',
                         'Anilios_nigrescens_R31022','Anilios_obtusifrons_R146400','Anilios_pilbarensis_R108813',
                         'Anilios_pinguis_R146995',
                         'Anilios_proximus_R132486', 
                         'Anilios_silvia_J46128','Anilios_splendidus_R119900','Anilios_systenos_R114894', 
                         'Anilios_torresianus_J47640','Anilios_tovelli_R21063','Anilios_troglodytes_R146048', 
                         'Anilios_unguirostris_R21669','Anilios_waitii_R113302','Anilios_wiedii_J59020', 
                         'Anilios_yirrikalae_J85695',
                         # 'Ramphotyphlops_cfwaitii_R51244', #omit
                         'Anilios_ligatus_R19109' # eastern
)

outgroups <- c('Ramphotyphlops_multillineatus_ABTC148379', 'Acutotyphlops_subocularis_R64768', 'Anilios_polygrammicus_R98715')

# Subset the tree
fos_tree <- ape::keep.tip(phy = fos_tree_full, tip = c(blindsnakes_samples, outgroups))

plotTree(fos_tree)
# Rename tips (for grypus we have three potential species and ligatus)
fos_tree$tip.label[which(fos_tree$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272', 'Anilios_ligatus_R31019'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_grypusET_R55272', 'Anilios_ligatusE_R31019')

# Change Anilios polygrammicus to Sundatyphlops
fos_tree$tip.label[which(fos_tree$tip.label %in% c('Anilios_polygrammicus_R98715', 'Ramphotyphlops_multillineatus_ABTC148379'))] <- c('Sundatyphlops_polygrammicus_R98715', 'Ramphotyphlops_multilineatus_R148379')

# Just keep the species name (without rego)
fos_tree$tip.label <- stringr::str_replace(string = fos_tree$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")
plotTree(fos_tree)


# write as newick for BioGeoBears -----------------------------------------

ape::write.tree(fos_tree, file = 'data/tree/anilios_newick.tre')


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
