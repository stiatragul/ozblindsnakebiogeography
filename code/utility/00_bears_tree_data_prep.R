# Utility/tree_data_match.R
# 2022-06-14

# For BioGeoBEARS analysis, we need a species-level time-calibrated phylogeny. 
# Will only include Australian species of Anilios as we don't have data on A. erycinus 

### NOTE ####
# DEC, DEC+J, and most methods in use in historical biogeography implicitly assume that the observed tree is the complete tree, 
# and that lineage-level speciation and extinction are independent of the range evolution process*

# The purpose of this script is to prepare the data required for BioGeoBEARS.
# --- Tree data -- SqCL exon-capture tree (two versions depends on prior)
# --- Anilios and close outgroups are considered

# libraries ---------------------------------------------------------------

library(dplyr); library(phytools); library(ape); library(stringr)

# phylogeny ---------------------------------------------------------------
# Time calibrated tree from SqCL probe and mcmctree analysis with fossils

# Full tree with fossils
phy_full <- ape::read.nexus("data/tree/anilios_v3_st.tre")
phy_full_b <- ape::read.nexus("data/tree/anilios_v3_b.tre")

# Rename cfwaitii
phy_full$tip.label[which(phy_full$tip.label == "Ramphotyphlops_cfwaitii_R51244")] <- "Anilios_centralis_R51244"
phy_full_b$tip.label[which(phy_full_b$tip.label == "Ramphotyphlops_cfwaitii_R51244")] <- "Anilios_centralis_R51244"

# Rename polygrammicus to Sundatyphlops
phy_full$tip.label[which(phy_full$tip.label == "Anilios_polygrammicus_R98715")] <- "Sundatyphlops_polygrammicus_R98715"
phy_full_b$tip.label[which(phy_full_b$tip.label == "Anilios_polygrammicus_R98715")] <- "Sundatyphlops_polygrammicus_R98715"


# Subset species tree to Acutotyphlops & Anilios ------------------------------------------

# Additional species to drop
# A. splendidus, we decided to synonymise this species with A. pinguis
# A. centralis R51244 is redundant
# only need one unguirostris at the moment. Will only keep 180002
to_drop <- c("Anilios_splendidus_R119900", "Anilios_centralis_R51244",
             "Anilios_unguirostris_R115861", "Anilios_unguirostris_R115861", "Anilios_unguirostris_R21669", "Anilios_unguirostris_R54430")

plot(phy_full);nodelabels()

phy_small <- phy_full %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

phy_small_b <- phy_full_b %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

phy_small <- phytools::force.ultrametric(phy_small,"extend")
phy_small_b <- phytools::force.ultrametric(phy_small_b,"extend")


# Rename tips (for grypus, ligatus we have three potential species)
phy_small$tip.label[which(phy_small$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297', 'Anilios_ligatus_R31019', 'Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_ligatusE_R31019', 'Anilios_grypusET_R55272')
phy_small_b$tip.label[which(phy_small_b$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297', 'Anilios_ligatus_R31019', 'Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_ligatusE_R31019', 'Anilios_grypusET_R55272')

# Rename R. multilineatus
phy_small$tip.label[which(phy_small$tip.label == 'Ramphotyphlops_multillineatus_ABTC148379')] <- 'Ramphotyphlops_multilineatus_A148379'
phy_small_b$tip.label[which(phy_small_b$tip.label == 'Ramphotyphlops_multillineatus_ABTC148379')] <- 'Ramphotyphlops_multilineatus_A148379'

# Just keep the species name (without rego)
phy_small$tip.label <- stringr::str_replace(string = phy_small$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")
phy_small_b$tip.label <- stringr::str_replace(string = phy_small_b$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")


# Check -------------------------------------------------------------------

plotTree(phy_small); plotTree(phy_small_b)

ape::write.tree(phy = phy_small, file = "data/intermediate_data/bears/blindsnake_st.tre")
ape::write.tree(phy = phy_small_b, file = "data/intermediate_data/bears/blindsnake_b.tre")
