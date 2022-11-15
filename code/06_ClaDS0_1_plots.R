# 06_ClaDS0_1_plots.R

# Cladogenetic diversification in rate shifts 
# Estimates lineage specific diversification rates on a phylogeny. 
# Accounts for diverse sources of variation in diversification rates that occur
# during the evolutionary history of clades. 

library(ape); library(phytools); library(RPANDA)
library(phyloch); data(strat2012)

# R fit ClaDS0 and ClaDS1
load('data/intermediate_data/ClaDS/clads_models.Rdata')

# ClaDS0 ------------------------------------------------------------------

plot_ClaDS0_chains(clads.out, burn = 1/2, thin = 1, 
                   param = c("sigma", "alpha", "l_0", "LP"))

clads0_rates <- getMAPS_ClaDS0(s_tree, clads.out, burn=0.5, thin=1)
clads0_rates

# 
pdf(file = 'output/ClaDS0_plot.pdf', width = 30, height = 40)
plot_ClaDS_phylo(s_tree, clads0_rates[-(1:3)], rates2 = NULL,
                 same.scale = T, main = NULL, 
                 lwd = 2, log = T, show.tip.label = F)
dev.off()


# ClaDS1 ------------------------------------------------------------------

plot_ClaDS_chains(clads1.out)
# Extract the maxima A posteriori for each parameter

maps = getMAPS_ClaDS(clads1.out, thin = 1)
print(paste0("sigma = ", maps[1], " ; alpha = ", 
             maps[2], " ; epsilon = ", maps[3], " ; l_0 = ", maps[4] ))

# Plot 
plot_ClaDS_phylo(phylo = s_tree, rates = maps[-(1:4)], 
                 main = "ClaDS1",
                 same.scale = T, lwd = 2, log = T, show.tip.label = F)
axisPhylo() # plots timescale

source('code/utility/plot_ClaDS_phylo_colour.R')
plot_ClaDS_phylo_colour(phylo = s_tree, rates = maps[-(1:4)], 
                        main = "ClaDS1",
                        same.scale = T, lwd = 2, log = T, show.tip.label = F)