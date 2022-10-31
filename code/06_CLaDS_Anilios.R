# 06_CLaDS_Anilios.R

# Cladogenetic diversification in rate shifts 
# Estimates lineage specific diversification rates on a phylogeny. 
# Accounts for diverse sources of variation in diversification rates that occur
# during the evolutionary history of clades. 

library(ape)
library(phytools)
library(RPANDA)
# library(diversitree)

full.tree <-read.tree("finalthesis_2mill_run_FigTree_run4_newick.tre")
full.tree

keep <- read.table("keep_tips_mcmctree_names.txt", header = T, row.names = NULL)

tree_sp <- full.tree$tip.label
table_sp <- keep[,2]
table_sp
matched_tips <- na.omit(match(table_sp,tree_sp))

drop_tips <- tree_sp[-matched_tips]

subset.tree.noout <- drop.tip(full.tree,drop_tips)

#CLADS

clads.out <- fit_ClaDS0(subset.tree.noout, name = NULL, 
                        pamhLocalName = "pamhLocal_nooutgroup",
                        iteration = 1e+07, thin = 20000, 
                        update = 1000,adaptation = 10, seed = 4821, nCPU = 1)
clads.out


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


color.bar(colorRampPalette(c("#ffffff","#ffffe0","#FEFCD7","#FCFFB2","#FBC17D","#FA8657","#ED504A","#C92D59","#981D69","#6B116F","#43006A","#1E0848","#080616"))(100), -1)

#example
MAPS = getMAPS_ClaDS0(ClaDS0_example$tree, 
                      ClaDS0_example$Cl0_chains, 
                      thin = 10)

plot_ClaDS_phylo(ClaDS0_example$tree, 
                 ClaDS0_example$speciation_rates, 
                 MAPS[-(1:3)])




