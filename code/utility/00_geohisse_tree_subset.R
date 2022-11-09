# 05_subset_st_tree.R
# November 2022

# libraries ---------------------------------------------------------------

library(ape)
library(dplyr); library(tidyr)

# Data --------------------------------------------------------------------

trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

fos_tree <- phytools::force.ultrametric(trees[[1]],"extend")
fos_tree$edge.length <- fos_tree$edge.length * 100
phy <- fos_tree

# Plot to check dates
plot(phy); axisPhylo()

# Check to make sure tree is ultrametric
is.ultrametric(phy)

# Prune -------------------------------------------------------------------

# List of species to drop from analysis. In this case we only want Anilios
droplist <- c("Acutotyphlops_subocularis", "Sundatyphlops_polygrammicus",
              "Sundatyphlops_polygrammicus", "Ramphotyphlops_multilineatus",
              "Anilios_splendidus")

pruned.tree <- ape::drop.tip(phy = phy, tip = droplist)

plot(pruned.tree); ape::axisPhylo()

# Write to file
write.tree(pruned.tree, file = "data/tree/subset_anilios_newick_st.tre")

