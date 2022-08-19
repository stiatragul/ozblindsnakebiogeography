# 00_diversification_tree_data_prep.R
# Putter Tiatragul
# Started on 27 May 2022

# The purpose of this script is to prepare the data required for diversification analyses.
# --- Tree data -- SqCL exon-capture tree (two versions depends on prior)
# --- Only Anilios taxa are considered.

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
             "Anilios_unguirostris_R115861", "Anilios_unguirostris_R115861", "Anilios_unguirostris_R21669", "Anilios_unguirostris_R54430",
             "Acutotyphlops_subocularis_R64768", "Sundatyphlops_polygrammicus_R98715", "Ramphotyphlops_multillineatus_ABTC148379")

plot(phy_full);nodelabels()

phy_small <- phy_full %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

phy_small_b <- phy_full_b %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

# Rename tips (for grypus, ligatus we have three potential species)
phy_small$tip.label[which(phy_small$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297', 'Anilios_ligatus_R31019', 'Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_ligatusE_R31019', 'Anilios_grypusET_R55272')
phy_small_b$tip.label[which(phy_small_b$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297', 'Anilios_ligatus_R31019', 'Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_ligatusE_R31019', 'Anilios_grypusET_R55272')

# Rename R. multilineatus
phy_small$tip.label[which(phy_small$tip.label == 'Ramphotyphlops_multillineatus_ABTC148379')] <- 'Ramphotyphlops_multilineatus_A148379'
phy_small_b$tip.label[which(phy_small_b$tip.label == 'Ramphotyphlops_multillineatus_ABTC148379')] <- 'Ramphotyphlops_multilineatus_A148379'

# Just keep the species name (without rego)
phy_small$tip.label <- stringr::str_replace(string = phy_small$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")
phy_small_b$tip.label <- stringr::str_replace(string = phy_small_b$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")


# Write out the trees -----------------------------------------------------


# Convert to multiphylo, not essential but nice...
plotTree(phy_small); plotTree(phy_small_b)
obj <- list(phy_small, phy_small_b)
class(obj) <- "multiPhylo"

ape::write.tree(phy = obj, file = "data/intermediate_data/diversification_analyses/blindsnake.trees")
