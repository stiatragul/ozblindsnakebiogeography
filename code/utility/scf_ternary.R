# name: sct_ternary.R
# Code from Ian to plot Ternary plot in baseR 

library(Ternary)
library(patchwork)
library(dplyr)

# gCF is defined as the percentage of “decisive” gene trees containing that branch
# sCF is a new measure defined as the percentage of decisive sites supporting a branch in the reference tree
# gCF and sCF complement classical measures of branch support in phylogenetics by providing a full description of underlying disagreement among loci and sites



# read in the CF file and appropriate tree (MCMC)
support = read.delim("data/tree/concord.cf.stat", header=T, comment.char="#") # mcmctree tree
mcmc.tree <- read.tree("data/tree/concord.cf.tree")


# setwd("/Users/Ian/Documents/GitHub/Egernia_Evolution")
#gCF.obj <- read.delim("gCF_sCF/T545_ASTRAL_MACSE.tre.cf.stat", header=T, comment.char="#")

# Source the color.by.CF script (attached to the email)
source("code/utility/color_by_CF_Ian.R")

# for ASTRAL
# CF.obj <- color.by.CF(astral.tree, cf.file="gCF_sCF/ASTRAL/T545_ASTRAL_MACSE.tre.cf.stat", value.check = T, legend=T, CF="sCF")
# for MCMCtree
CF.obj <- color.by.CF(mcmc.tree, cf.file="data/tree/concord.cf.stat", value.check = T, legend=T, CF="sCF")


# plot sCF as a function of gCF
SxG <- CF_comparison(CF.obj, type="sCFxgCF")
# plot sCF as a function of Branch Length
SxBL <- CF_comparison(CF.obj, type="sCFxBL")
# plot gCF as a function of Branch Length
GxBL <- CF_comparison(CF.obj, type="gCFxBL")

SxG
SxBL
GxBL

# put the plots together with patchwork
SxG + SxBL + GxBL

# Plot site concordance factors in a ternary plot
# choose the site-concordance-factor columns
sCF <- dplyr::select(support, "sCF", "sDF1", "sDF2"); sCF <- sCF[2:nrow(sCF),]
# set the background colors for point density
plotCols <- colorRampPalette(brewer.pal(6, "RdYlBu"), alpha=T); point.colors <- rev(plotCols(100))
# plot the ternary background
Ternary::TernaryPlot(alab = 'sCF', blab = 'sDF1', clab = 'sDF2')
#ColourTernary(TernaryDensity(sCF, resolution = 10))
# plot the ternary point density colors
Ternary::ColourTernary(TernaryDensity(sCF, resolution = 10), spectrum=point.colors)

# plot the points, but start with ones where sCF > sDF1 & sDF2
Ternary::TernaryPoints(filter(sCF, sCF > sDF1 & sCF > sDF2), col = "white", pch = 16, cex=1)
# plot points where sCF < sDF1 & 2
Ternary::TernaryPoints(filter(sCF, sCF < sDF1 | sCF < sDF2), col = "black", pch = 16, cex=1)
# plot points where sCF~sDF1~sDF2
Ternary::TernaryPoints(filter(sCF, sCF < 40 & sCF > 30 &
                          sDF1 < 40 & sDF1 > 30 &
                          sDF2 < 40 & sDF2 > 30), col = "grey", pch = 16, cex=1)



dplyr::filter(sCF, sDF1 > 70)
dplyr::filter(sCF, sDF2 > 70)
dplyr::filter(sCF, sCF < 20)


# investigate which the nodes are where an alternate topology is supported by CFs
node_df <- dplyr::filter(scF, sCF < sDF1 | sCF < sDF2)


# investigate which nodes are equally supported in 3 topologies
node_df_3topo <- dplyr::filter(support, sCF < 40 & sCF > 30 &
         sDF1 < 40 & sDF1 > 30 &
         sDF2 < 40 & sDF2 > 30)

dplyr::filter(support, sDF1 > 70)
dplyr::filter(support, sDF2 > 70)
dplyr::filter(support, sCF < 30)


getDescendants(mcmc.tree, node = 348)
getDescendants(mcmc.tree, node = 248)

mcmc.tree$tip.label[getDescendants(mcmc.tree, node = 348)]
mcmc.tree$tip.label[getDescendants(mcmc.tree, node = 248)]
mcmc.tree$tip.label[getDescendants(mcmc.tree, node = 269)]
mcmc.tree$tip.label[getDescendants(mcmc.tree, node = 269)]
mcmc.tree$tip.label[getDescendants(mcmc.tree, node = 361)]



plot(mcmc.tree)
nodelabels(node = c(node_df_3topo$ID))

dev.off()
plotTree(mcmc.tree, fsize = 0.5, lwd = 0.6)
add.arrow(tip = c(node_df_3topo$ID), arrowl =1)

nodelabels(node = c(node_df_3topo$ID), 
           adj = c(1.2, -0.2), frame = "n", cex = 0.8)
# nodelabels(node = c(node_df$ID[2]))

nodelabels(node = c(node_df_3topo$ID[7]))


nrow(filter(support, sDF1 > sCF & sDF2 > sCF))
nrow(filter(support, sDF1 > sCF & sDF2 < sCF))

