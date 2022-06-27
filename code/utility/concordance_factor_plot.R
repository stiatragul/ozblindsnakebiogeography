library(dplyr)
library(patchwork)

source('code/utility/color_by_CF_Ian.R')

# read in the CF file and appropriate tree (ASTRAL)
support = read.delim("data/tree/concord.cf.stat", header=T, comment.char="#") # astral tree
astral.tree <- read.tree("data/tree/concord.cf.tree")
# astral.tree <- read.tree("data/tree/con")

# Visualize the concordance factors as branch colors
dev.off()

# s1 <- astral.tree$tip.label[stringr::str_detect(astral.tree$tip.label, pattern = "Anilios")]
# s2 <- astral.tree$tip.label[stringr::str_detect(astral.tree$tip.label, pattern = "Acutotyphlops")]
# s3 <- astral.tree$tip.label[stringr::str_detect(astral.tree$tip.label, pattern = "Ramphotyphlops")]
# s4 <- astral.tree$tip.label[stringr::str_detect(astral.tree$tip.label, pattern = "Indotyphlops")]
# subset_things <- c(s1, s2, s3, s4)


# CF.obj <- color.by.CF(astral.tree, cf.file="data/tree/concord.cf.stat", value.check = T, legend=T, CF="gCF", subset = TRUE, subset_tips = subset_things)
CF.obj <- color.by.CF(astral.tree, cf.file="data/tree/concord.cf.stat", value.check = T, legend=T, CF="gCF", subset = FALSE)


CF.scf <- color.by.CF(astral.tree, cf.file="data/tree/concord.cf.stat", value.check = T, legend=T, CF="sCF")
# CF.obj <- color.by.CF(mcmctree_anilios, cf.file="data/tree/concord_ahe.cf.stat", value.check = T, legend=T, CF="gCF")


dev.off()

length(astral.tree$tip.label)

# Plot gene concordance factors as a result of branch lengths
CF_comparison(CF.obj, type="gCFxBL")
CF_comparison(CF.obj, type="sCFxBL")
CF_comparison(CF.obj, type="sCFxgCF")



