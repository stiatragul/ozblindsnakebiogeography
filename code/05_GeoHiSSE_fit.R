# 05_GeoHiSSE_fit.R
# October 2022

# Are rates equal between arid and non-arid zones? 
# How does it differ? Is it better explained by hidden states? 

# install -----------------------------------------------------------------

# library( devtools )
# install_github(repo = "thej022214/hisse", ref = "master")

# libraries ---------------------------------------------------------------

library(hisse); library(ape)
library(diversitree)
library(ape); library(RColorBrewer)
library(dplyr); library(tidyr)
library(ggplot2)

# Data --------------------------------------------------------------------

trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

fos_tree <- phytools::force.ultrametric(trees[[2]],"extend")
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
write.tree(pruned.tree, file = "data/tree/subset_anilios_newick_b.tre")

# Re read in the subset tree
phy <- read.tree("data/tree/subset_anilios_newick_b.tre")

taxon_list <- data.frame(taxon = sort(phy$tip.label))

# write.csv(taxon_list, file = "data/intermediate_data/geohisse/taxon_list.csv", row.names = FALSE)


# State data --------------------------------------------------------------
# State data indicates whether the species is classified as (state 1) arid, (state 2) non-arid, (state 0) widespread between both regions.
# Classification is based on distribution of 'species'. I plotted ALA distribution on biome map of Australia in QGIS. 

# Test whether arid habitation influences rates of diversification 
# and to reconstruct ancestral states while accounting for diversification rate heterogeneity.

state <- read.csv("data/intermediate_data/geohisse/arid_nonarid_both_states.csv", header = T)
state <- data.frame(taxon=state$taxon, ranges=as.numeric(state$state))
head(state)

state$new_tip_names <- paste(state$taxon, state$ranges, sep = "_")

# Rename tips

state_temp <- data.frame(tip.labels = phy$tip.label)

new_tip_names <- state_temp %>% 
  dplyr::left_join(state, by = c("tip.labels" = "taxon")) %>% 
  select(new_tip_names)

phy$tip.label <- new_tip_names$new_tip_names
plot(phy)

state <- data.frame(taxon = state$new_tip_names, ranges = state$ranges)
head(state)

## Proportion of species for input in `f` argument, estimated proportion of extant species in states
prop_df <- read.csv("data/intermediate_data/geohisse/2022_species_list_arid_nonarid_widespread.csv", header = T)

prop_df$state <- as.character(prop_df$state)

proportion_f <- prop_df %>% 
  dplyr::group_by(state, in_tree) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) %>% 
  dplyr::filter(in_tree == "yes")
  
proportion_f

# Model fitting -----------------------------------------------------------

# Setting up the models
# We will fit a total of four models. Two models with a range-independent diversification 
# process and two other models in which the range have an effect on the diversification rate 
# of the lineages (each with either one or two rate classes).

### ASSUMPTIONS and PARAMETERS ###

#turnover - turnover parameter tau00, tau11, tau01 
#f = proportion of extant species in the different state (areas). 
#eps = extinct fraction parameters. 

#porportion of extant species probably not fully sampled
fraction <- c(0.89, # arid
              0.71, # non-arid
              0.83) # widespread

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
#### Assume equal rates regardless of biogeographic region ###

turnover <- c(1,1,0) # widespread range is removed from the model. Diversification rates are independent range evolution
eps <- c(1,1) 
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <- GeoHiSSE(phy = phy, data = state, f=fraction,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod, 
                 turnover.upper=100, trans.upper=10,
                 assume.cladogenetic = TRUE)

## Model 2. Canonical GeoSSE model, range effect on diversification
#### Range evolution affects diversification by adding turnover rate for widespread range. 
# Three turnover parameters and two extinction fraction parameters. 
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod2 <- GeoHiSSE(phy = phy, data = state, f=fraction,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

## Hidden states
## Model 3. GeoHiSSE model with 1 hidden trait, no range-dependent diversification.
## Note below how parameters vary among hidden classes but are the same within each
##      hidden class.
# tau00A, tau11A, tau01A, tau00B, tau11B, and tau01B

turnover <- c(1,1,0,2,2,0)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1, make.null=TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod3 <- GeoHiSSE(phy = phy, data = state, f=fraction,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)

mod3$index.par

## Model 4. GeoHiSSE model with 1 hidden trait, range-dependent diversification.
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod4 <- GeoHiSSE(phy = phy, data = state, f=fraction,
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


## Model 5. MuSSE-like model with no hidden trait, no cladogenetic effects.
# Fit a complementary set of models that remove the cladogenetic effect entirely, such that all changes occur along branches (i.e., anagenetic change). 
# This requires the removal of the turnover rate for lineages in the widespread range and ensuring that range contraction is distinct from the extinction of endemics:
turnover <- c(1,2,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0, make.null=FALSE, 
                                    separate.extirpation = TRUE)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
trans.rate.mod <- ParEqual(trans.rate.mod, c(2,3))
mod5 <- GeoHiSSE(phy = phy, data = state, f=fraction,
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10, # need to change trans.upper to 11 instead of 10 for model to run. 
                 sann=FALSE, 
                 assume.cladogenetic = FALSE)


# Get AIC -----------------------------------------------------------------
# Akaike weights are important to evaluate the relative importance of each of
# the models to explain the variation observed in the data. 
# Accounts for penalties associated to the number of free parameters.
# Models with higher weight show better fit to the data and, as a result, 
# have more weight when performing model averaging.

hisse::GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3,
                          model4 = mod4, model5 = mod5), criterion="AICc")


### Model averaging and plotting

# Model average the results. Reflects Akaike model weights.

# Marginal reconstruction for each models
# Reconstructs hidden states at the nodes of the phylogeny. 

recon.mod1 <- MarginReconGeoSSE(phy = mod1$phy, data = mod1$data, f = mod1$f,
                                pars = mod1$solution, hidden.states = 1,
                                root.type = mod1$root.type, root.p = mod1$root.p,
                                AIC = mod1$AIC, n.cores = 1)
recon.mod2 <- MarginReconGeoSSE(phy = mod2$phy, data = mod2$data, f = mod2$f,
                                pars = mod2$solution, hidden.states = 1,
                                root.type = mod2$root.type, root.p = mod2$root.p,
                                AIC = mod2$AIC, n.cores = 1)
recon.mod3 <- MarginReconGeoSSE(phy = mod3$phy, data = mod3$data, f = mod3$f,
                                pars = mod3$solution, hidden.states = 2,
                                root.type = mod3$root.type, root.p = mod3$root.p,
                                AIC = mod3$AIC, n.cores = 1)
recon.mod4 <- MarginReconGeoSSE(phy = mod4$phy, data = mod4$data, f = mod4$f,
                                pars = mod4$solution, hidden.states = 2,
                                root.type = mod4$root.type, root.p = mod4$root.p,
                                AIC = mod4$AIC, n.cores = 1)
recon.mod5 <- MarginReconGeoSSE(phy = mod5$phy, data = mod5$data, f = mod5$f,
                                pars = mod5$solution, hidden.states = 1,
                                root.type = mod5$root.type, root.p = mod5$root.p,
                                AIC = mod5$AIC, n.cores = 1)

# Marginal reconstruction computed above is used to compute the model average.

recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4, recon.mod5)

# save.image(file = 'data/intermediate_data/geohisse/geohisse_three_states.Rdata')
load('data/intermediate_data/geohisse/geohisse_three_states.Rdata')


# See matrix with parameter estimates for each species averaged over all models.
# For each tip (as indicated)
# model.ave.rates <- hisse::GetModelAveRates(recon.models, type = "both", 
#                                     bound.par.matrix=cbind(c(0,0,0,0,0),
#                                                            c(10000000,10000000,10000000,10000000,10000000)) )
model.ave.rates <- hisse::GetModelAveRates(recon.models, type = "both")
head( model.ave.rates )

model.ave.rates$tips

# add new column for range and change species name
rates_df <- model.ave.rates$tips %>% 
  tidyr::separate(col = taxon, into = c("genus", "species", "range"), sep = "_")

rates_df$species <- paste("A.", rates_df$species, sep = "_")
rates_df$range[which(rates_df$range == 0)] <- "widespread"
rates_df$range[which(rates_df$range == 1)] <- "arid"
rates_df$range[which(rates_df$range == 2)] <- "non-arid"

statecolours <- c("#f03b20", "#ffeda0", "#feb24c")
# statecolours <- c("#e7564d",
                  # "#5ec284",
                  # "#fef6ff")

pdf(file = 'output/three_states.pdf')
boxplot(rates_df$net.div ~ rates_df$range, 
        rates_df,
        las = 1, ylab="Net Diversification", col=statecolours)

boxplot(rates_df$speciation ~ rates_df$range, 
        rates_df,
        las = 1, ylab="Speciation", col=statecolours)

boxplot(rates_df$extinction ~ rates_df$range, 
        rates_df, notch=FALSE , las = 1, ylab="Extinction",
        col=statecolours)

boxplot(rates_df$turnover ~ rates_df$range,
        rates_df, notch=FALSE , las = 1, ylab="Turnover", 
        col=statecolours)
dev.off()

# rates -------------------------------------------------------------------

# (state 1) arid, (state 2) non-arid (1), (state 0) widespread
# statecolours <- c("#EEA47FFF", "#00539CFF", "black")
# statecolours <- c("#f03b20", "#ffeda0", "#feb24c")
statecolours <- c("#e7564d",
                  "#5ec284",
                  "#fef6ff")
ratecolours <- colorRampPalette(brewer.pal(6, 'RdYlBu'))(100)

plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "phylogram",
                     show.tip.label = TRUE, legend = T,
                     legend.cex = 1,
                     # rate.colors = ratecolours,
                     # state.colors = statecolours, fsize = 0.8, 
                     outline.color = 'black')

model.ave.rates

## Plot tree and label clades
# phytools::plotTree(phy, type = "fan", ftype = "off")
plot.geohisse.states(x = recon.models, rate.param = "turnover", type = "fan",
                     show.tip.label = T, legend = T,
                     # rate.colors = ratecolours,
                     state.colors = statecolours,
                     fsize = 0.8)
# phytools::arc.cladelabels(text="clade A", node=phytools::findMRCA(phy, c("Anilios_pilbarensis_1", "Anilios_endoterus_1")),
                          # ln.offset=1.1, lab.offset=1.16, lwd = 6)

# # Net diversification based on range
# rates_df %>% 
#   ggplot(., aes(x = range, y = net.div)) +
#   geom_boxplot() +
#   theme_classic()
# 
# rates_df %>% 
#   ggplot(., aes(x = range, y = speciation)) +
#   geom_boxplot() +
#   theme_classic()


# AIC ---------------------------------------------------------------------

recon.mod1$AIC; recon.mod2$AIC; recon.mod3$AIC; recon.mod4$AIC; recon.mod5$AIC
mod1$AIC; mod2$AIC; mod3$AIC; mod4$AIC; mod5$AIC

mod1$AICc; mod2$AICc; mod3$AICc; mod4$AICc; mod5$AICc

recon.mod1$rates.mat
recon.mod2$rates.mat
recon.mod3$rates.mat
recon.mod4$rates.mat

