# geoHiSSE

# library( devtools )
# install_github(repo = "thej022214/hisse", ref = "master")

# library(hisse)
# library(diversitree)

phy <- read.tree("subset_2mill.tree")
library(ape)

# put the names of tips you want to drop in a vector

droplist<-c(
  "Chiropodomys_gliroides_Z25153_EC1",
  "Chiropodomys_karlkoopmani_NHM19782959_EC4",
  "Chiropodomys_muroides_BMNH290_EC4",
  "Haeromys_minahassae_Z23000_EC1",
  "Apomys_lubangensis_JAE169_EXOME",
  "Archboldomys_luzonensis_EAR1826_EXOME",
  "Chrotomys_mindorensis_JAE520_EXOME",
  "Rhynchomys_isarogensis_JAE2195_EXOME")

pruned.tree <- drop.tip(phy,phy$tip.label[match(droplist, phy$tip.label)])
write.tree(pruned.tree, file="subset_2mill_just_sahul.tree")


phy <- read.tree("subset_2mill_just_sahul.tree")
phy$tip.label

state <- read.csv("shelf_vs_island_binary.txt",header=F)
state  

state <- data.frame(taxon=state$V1, ranges=as.numeric(state$V2))
state


# Setting up the models
# We will fit a total of four models. Two models with a range-indendent diversification process and 
# two other models in which the range have an effect on the diversification rate of the lineages (each with either one or two rate classes).

#turnover - turnover parameter tau00, tau11, tau01
#f = estimated proportion of extant species in the different state (areas)
#eps = extinct fraction parameters. 

## Model 1 - Dispersal parameters vary only, no range-dependent diversification.
#### Assume equal rates regardless of biogeographic region. ###
turnover <- c(1,1,0)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod1 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod, 
                 turnover.upper=100, trans.upper=10)

## Model 2. Canonical GeoSSE model, range effect on diversification
#### Range evolution affects diversification by adding turnover rate for widespread range. 
# Three turnover parameters and two extinction fraction parameters. 
turnover <- c(1,2,3)
eps <- c(1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=0)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod2 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
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
mod3 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
                 turnover=turnover, eps=eps,
                 hidden.states=TRUE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10)


## Model 4. GeoHiSSE model with 1 hidden trait, range-dependent diversification.
turnover <- c(1,2,3,4,5,6)
eps <- c(1,1,1,1)
trans.rate <- TransMatMakerGeoHiSSE(hidden.traits=1)
trans.rate.mod <- ParEqual(trans.rate, c(1,2))
mod4 <- GeoHiSSE(phy = phy, data = sim.dat, f=c(1,1,1),
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
mod5 <- GeoHiSSE(phy = phy, data = state, f=c(1,1,1),
                 turnover=turnover, eps=eps,
                 hidden.states=FALSE, trans.rate=trans.rate.mod,
                 turnover.upper=100, trans.upper=10, sann=FALSE, 
                 assume.cladogenetic = FALSE)



load( "geohisse_new_vignette.Rsave" )
GetAICWeights(list(model1 = mod1, model2 = mod2, model3 = mod3, model4 = mod4), criterion="AIC")

### Model averaging and plotting

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




recon.models <- list(recon.mod1, recon.mod2, recon.mod3, recon.mod4)
#model.ave.rates <- GetModelAveRates(x = recon.models, type = "tips")

#head( model.ave.rates )

#plot.geohisse.states(x = recon.models, rate.param = "net.div", type = "fan",
#                    show.tip.label = FALSE, legend = FALSE)


recon.mod1$AIC
recon.mod2$AIC
recon.mod3$AIC
recon.mod4$AIC

recon.mod1$rates.mat
recon.mod2$rates.mat
recon.mod3$rates.mat
recon.mod4$rates.mat

