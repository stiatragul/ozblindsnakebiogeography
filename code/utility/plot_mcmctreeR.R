# Plot MCMCTREE

library(MCMCtreeR)
library(ape)
library(phytools)


phy_full_mcmc <- readMCMCtree(inputPhy = 'data/tree/v3_01_FigTree.tre', from.file = TRUE)
MCMC.chain <- read.delim('data/tree/v3_01_combo.txt') # This text is very large. Need to find it from local computer (not committed to github)

phy_full <- phy_full_mcmc$apePhy 

# Rename tips (for grypus we have three potential species and ligatus)
phy_full$tip.label[which(phy_full$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_grypusET_R55272')

# Change Anilios polygrammicus to Sundatyphlops
phy_full$tip.label[which(phy_full$tip.label %in% c('Anilios_polygrammicus_R98715', 'Ramphotyphlops_multillineatus_ABTC148379'))] <- c('Sundatyphlops_polygrammicus_R98715', 'Ramphotyphlops_multilineatus_R148379')


# Rename some taxa in full tree
phy_full_mcmc$apePhy$tip.label <- phy_full$tip.label

# Full phylogeny with node ages uncertainty
# Include this in supplementary material
# dev.off()
# pdf(file = 'output/mcmctree_full_age.pdf', paper = "a4" )
# 
MCMC.tree.plot(phy_full_mcmc, analysis.type = "MCMCtree", cex.tips = 0.6,
               time.correction = 100, plot.type = "phylogram", lwd.bar = 2,
               scale.res = c("Eon", "Period"), node.method = "bar", col.age = "navy",
               no.margin = TRUE, label.offset = 4)
# 
# dev.off()


# Plot with posterior -----------------------------------------------------

MCMC.chain <- MCMC.chain[c(1, grep(x = colnames(MCMC.chain), pattern = 't_n'))]
phy.edge <- phy_full$edge

# extract ages with node age posteriors from column 2
MCMC.node.ages <- MCMC.chain[, 2:Ntip(phy_full)]

# extract the names from these data so they are numeric node
# labels that match the APE tree format
all.nodes <- as.numeric(gsub("t_n", "", colnames(MCMC.node.ages)))

# create a vector of names for each list element as internal
# nodes from APE tree, using phy$edge object.
node.ages.names <- c(Ntip(phy_full) + 1, phy.edge[which(phy.edge[, 2] > Ntip(phy_full)), 2])
# find where each posterior node age appears in APE edge
# object.
match.nodes <- match(all.nodes, as.numeric(node.ages.names))
# create a list extracting the information from the MCMC
# chain in APE node order.
node.ages <- lapply(match.nodes, function(uu) MCMC.node.ages[, 
                                                             uu])
# name each element in list.
names(node.ages) <- node.ages.names

dev.off()
pdf(file = 'output/mcmctree_full_age_distribution.pdf', paper = "a4" )
MCMC.tree.plot(phy = phy_full, node.ages = node.ages, analysis.type = "user", 
               cex.tips = 0.2, time.correction = 100, scale.res = c("Eon", 
                                                                    "Period"), plot.type = "distributions", cex.age = 0.4, 
               cex.labels = 0.5, relative.height = 0.08, col.tree = "grey40", 
               no.margin = TRUE)

dev.off()

# Subset species tree to Anilios ------------------------------------------

plot(phy_full)
nodelabels()

phy_small <- ape::extract.clade(phy = phy_full, 159)

plot(phy_small)
nodelabels()

# 
# node_ages <- as.data.frame(phy_full_mcmc$nodeAges)
# node_ages$row <- rownames(node_ages)


phy_small$edge
node.numbers <- c(47:91) # corresponds to 159:203 in the original phylogeny. This is without taking away any internal taxa.
node.numbers_old <- c(159:203)

# create a list to store each posterior sample for every node
node.posteriors <- vector(mode = "list", length = Nnode(phy_small))
names(node.posteriors) <- node.numbers

node.ages$`159`

#  list of posterior from full analysis
class(node.ages)

match.nodes_s <- match(node.numbers_old, as.numeric(node.ages.names))
node.ages_new <- lapply(match.nodes_s, function(uu) MCMC.node.ages[, uu])

names(node.ages_new) <- as.character(node.numbers)

names(node.ages_new)


# dev.off()
# pdf(file = 'output/mcmctree_anilios_age_bar.pdf', paper = "a4" )
MCMC.tree.plot(phy = phy_small, node.ages = node.ages_new, analysis.type = "user", 
               cex.tips = 0.6, time.correction = 100, 
               scale.res = c("Eon", "Period"), plot.type = "phylogram", 
               cex.age = 0.6, lwd.bar = 1, cex.labels = 0.6, 
               relative.height = 0.08, col.tree = "grey40", col.age = "navy",
               node.method = "bar",
               no.margin = TRUE)
# dev.off()


