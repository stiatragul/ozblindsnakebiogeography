# Plot MCMCTREE

library(MCMCtreeR)
library(ape)
library(phytools)


phy_full_mcmc <- readMCMCtree(inputPhy = 'data/tree/v3_01_FigTree.tre', from.file = TRUE)
MCMC.chain <- read.delim('data/tree/v3_01_combo.txt') 

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


# Subset species tree to Anilios ------------------------------------------

# tips we want for blindsnakes (the best we can do since we do not have all described species)
blindsnakes_samples <- c('Anilios_affinis_J93635','Anilios_ammodytes_R158097','Anilios_aspina_J91822','Anilios_australis_R115859',
                         'Anilios_bicolor_R55342','Anilios_bituberculatus_R63990','Anilios_broomi_J87805','Anilios_centralis_R14631',
                         'Anilios_diversus_R112027','Anilios_endoterus_R102725','Anilios_ganei_R165000', 
                         'Anilios_grypusW_R108596',   # Western Australia (W)
                         'Anilios_grypusNW_R157297',   # Sister to leptosoma-complex (NW)
                         'Anilios_grypusET_R55272',    # Eastern (ET)
                         # 'Anilios_grypus_R108596',   # Western Australia (W)
                         # 'Anilios_grypus_R157297',   # Sister to leptosoma-complex (NW)
                         # 'Anilios_grypus_R55272',    # Eastern (ET)
                         'Anilios_guentheri_R108431','Anilios_hamatus_R136277','Anilios_howi_R146381','Anilios_kimberleyensis_R164213', 
                         'Anilios_leptosoma_R119241','Anilios_leucoproctus_J87547', 
                         'Anilios_ligatus_R31019','Anilios_longissimus_R120049','Anilios_margaretae_R163269',
                         'Anilios_nigrescens_R31022','Anilios_obtusifrons_R146400','Anilios_pilbarensis_R108813',
                         'Anilios_pinguis_R146995', 'Anilios_proximus_R132486', 'Anilios_silvia_J46128','Anilios_splendidus_R119900','Anilios_systenos_R114894', 
                         'Anilios_torresianus_J47640','Anilios_tovelli_R21063','Anilios_troglodytes_R146048', 
                         'Anilios_unguirostris_R21669','Anilios_waitii_R113302','Anilios_wiedii_J59020', 
                         # 'Anilios_ligatus_R19109', # eastern
                         'Anilios_yirrikalae_J85695'
                         # 'Ramphotyphlops_cfwaitii_R51244', #omit
)

outgroups <- c('Sundatyphlops_polygrammicus_R98715', 'Acutotyphlops_subocularis_R64768', 'Ramphotyphlops_multilineatus_R148379')
# c('Ramphotyphlops_multillineatus_ABTC148379', 'Acutotyphlops_subocularis_R64768', 'Anilios_polygrammicus_R98715')

# Subset ------------------------------------------------------------------
# Subset the tree

phy_small <- ape::keep.tip(phy = phy_full, tip = c(blindsnakes_samples, outgroups))

head(MCMC.chain)
phy.edge <- phy_small$edge

MCMCtree.posterior

node_ages <- as.data.frame(phy_full_mcmc$nodeAges)
node_ages$row <- rownames(node_ages)
subset_node_age <- node_ages[which(node_ages$row >= 159), ]

# create a list to store each posterior sample for every node
node.posteriors <- vector(mode = "list", length = Nnode(phy_small))

node.numbers <- c(47:91)
names(node.posteriors) <- node.numbers
for (i in 1:Nnode(phy_small)) node.posteriors[[i]] <- rnorm(1000, 
                                                            mean = subset_node_age$mean[i], sd = 0.1)

node.posteriors

MCMC.tree.plot(phy = phy_small, analysis.type = "user", node.ages = node.posteriors,
               cex.tips = 0.7, time.correction = 100, 
               scale.res = c("Eon", "Period"), plot.type = "distributions")



