library(dplyr)
library(patchwork)

source('code/utility/color_by_CF_Ian.R')

# read in the CF file and appropriate tree (ASTRAL)
support = read.delim("data/tree/concord.cf.stat", header=T, comment.char="#") # astral tree
astral.tree <- read.tree("data/tree/20220617_all_loci.tre")

# Visualize the concordance factors as branch colors
CF.obj <- color.by.CF(astral.tree, cf.file="data/tree/concord_ahe.cf.stat", value.check = T, legend=T, CF="gCF")

# Plot gene concordance factors as a result of branch lengths
CF_comparison(CF.obj, type="gCFxBL")
