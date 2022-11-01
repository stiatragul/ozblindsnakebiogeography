# 06_CLaDS_Anilios.R

# Cladogenetic diversification in rate shifts 
# Estimates lineage specific diversification rates on a phylogeny. 
# Accounts for diverse sources of variation in diversification rates that occur
# during the evolutionary history of clades. 

library(ape); library(phytools); library(RPANDA)

# Tree --------------------------------------------------------------------
s_tree <- ape::read.tree("data/tree/subset_anilios_newick_b.tre")

# ClaDS -------------------------------------------------------------------

# Estimate branch-specific diversification rates
# RPANDA v 2.0

sample_fraction = 37/47 # fraction of tips

# ClaDS0 -- no extinction, mu = 0

clads.out <- fit_ClaDS0(s_tree, name = 'data/intermediate_data/ClaDS/clads0.Rdata', 
                        pamhLocalName = "pamhLocal_ClaDS0",
                        iteration = 1e+07, thin = 20000, 
                        update = 1000, adaptation = 10, seed = 4821, nCPU = 3)

# ClaDS1 -- constant extinction rates across the tree.

clads1.out <- fit_ClaDS(s_tree, sample_fraction, 
                        iterations = 1e+07, thin = 20000,
                        # iterations = 1000, thin = 50, 
                        model_id = "ClaDS1", seed = 4821, nCPU = 3,
                        file_name = 'data/intermediate_data/ClaDS/clads1.Rdata')

# ClaDS2 -- turnover = ext/speciation, remains constant across the tree

clads2.out <- fit_ClaDS(s_tree, sample_fraction, 
                        iterations = 1e+07, thin = 20000,
                        # iterations = 1000, thin = 50, 
                        model_id = "ClaDS2", seed = 4821, nCPU = 3,
                        file_name = 'data/intermediate_data/ClaDS/clads2.Rdata')


save.image('data/intermediate_data/ClaDS/clads_models.Rdata')