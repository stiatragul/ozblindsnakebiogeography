# 06_CLaDS_plot_tips.R

# Cladogenetic diversification in rate shifts 
# Estimates lineage specific diversification rates on a phylogeny. 
# Accounts for diverse sources of variation in diversification rates that occur
# during the evolutionary history of clades. 

# Plot the MCMC chains obtained with fit_ClaDS.

library(ape); library(phytools); library(RPANDA)

# Load data ---------------------------------------------------------------
s_tree <- ape::read.tree("data/tree/subset_anilios_newick_b.tre")

# load data from 06_ClaDS_Anilios.R


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







rates <- getMAPS_ClaDS0(subset.tree.noout, clads.out, burn=0.5, thin=1)
rates

plot_ClaDS0_chains(clads.out, burn = 1/2, thin = 1, 
                   param = c("sigma", "alpha", "l_0", "LP"))


plot_ClaDS_phylo(subset.tree.noout, rates[-(1:3)], rates2 = NULL,same.scale = T, main = NULL, lwd = 2, log = T, show.tip.label = F)

plot_ClaDS_phylo(subset.tree, rates[-(1:3)], rates2 = NULL,same.scale = T, main = NULL, lwd = 2, log = T, show.tip.label = F)


#plot with magma palette
source("plot_ClaDS_phylo_colours.R")

plot_ClaDS_phylo_col(subset.tree.noout, rates[-(1:3)], rates2 = NULL,same.scale = T, main = NULL, lwd = 2, log = T, show.tip.label = F)
axisPhylo() # plots timescale


pdffn = "CLaDs_Sahul_magma1.pdf"

pdf(pdffn, width=30, height=40)

plot_ClaDS_phylo_col(subset.tree.noout, rates[-(1:3)], rates2 = NULL,same.scale = T, main = NULL, lwd = 10, log = T, show.tip.label = F)
axisPhylo() # plots timescale


dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it


#make scale bar
Colors = colorRampPalette(c("#ffffff","#ffffe0","#FEFCD7","#FCFFB2","#FBC17D","#FA8657","#ED504A","#C92D59","#981D69","#6B116F","#43006A","#1E0848","#080616")) 
Colors

color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  
  dev.new(width=1.75, height=5)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}


color.bar(colorRampPalette(c("#ffffff","#ffffe0","#FEFCD7","#FCFFB2","#FBC17D",
                             "#FA8657","#ED504A","#C92D59","#981D69","#6B116F",
                             "#43006A","#1E0848","#080616"))(100), -1)

