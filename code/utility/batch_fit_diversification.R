# batch_fit_diversification.R

# Fit diversification by environmental data rather than full list

# Library -----------------------------------------------------------------

library(RPANDA); library(tidyverse)
library(phytools); library(gridExtra)
library(parallel) # can do this on WSL2 (if using Windows)

# Load --------------------------------------------------------------------

load(file = "data/intermediate_data/diversification_analyses/env_data_list.RData")
trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

is.ultrametric(trees[[1]])
fos_tree <- phytools::force.ultrametric(trees[[1]],"extend")
plotTree(fos_tree)
is.ultrametric(fos_tree)

# Fit models ---------------------------------------------------------------

# Crown age
tot_time <- max(node.age(fos_tree)$ages)
cond = "crown"

## Fraction of diversity represented. In our tree we have 
## 48 Described species of Anilios + 2 more grypus undescribed + 1 more ligatus - A. splendidus
## So our "true" species at the moment is 50 and we have 41 tips in our tree

length(fos_tree$tip.label)
frac = 38/50

# Explanation of the models we'll fit:
# BCST
#   + standard Yule model, speciation rate is constant and no extinction
#	BCST_DCST
#		+ standard birth-death model, with constant speciation and extinction rates
#	BEnvVar_EXPO
#		+ speciation rate is correlated exponentially to the rate of environmental variation, with no extinction (mu=0)
#	BEnvVarDCST_EXPO
#		+ speciation rate is correlated exponentially to the rate of environmental variation, with constant extinction
#	BCSTDEnvVar_EXPO
#		+ speciation rate is constant, and extinction rate is correlated exponentially to the rate of environmental variation
#	BEnvVarDEnvVar
#		+ speciation and extinction rates are correlated exponentially to the rate of environmental variation

#####################################################
################ Reference models ###################
#####################################################

# BCST (YULE) Speciation rate constant through time with no extinction (mu = 0)
f.lamb <- function(t,y){y[1]}
f.mu <- function(t,y){0}
lamb_par <- c(0.1) 
mu_par <- c(0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=T

tree_BCST <- fit_bd(phylo = fos_tree, tot_time = tot_time, 
                    f.lamb = f.lamb, 
                    f.mu, lamb_par, 
                    mu_par = mu_par, f = frac, 
                    cst.lamb=cst.lamb, cst.mu=cst.mu, expo.lamb=expo.lamb, 
                    expo.mu=expo.mu, fix.mu=fix.mu, dt = 1e-3,
                    cond="crown")

tree_BCST$model <- "BCST"

#	BCST_DCST (constant birth-death (Yule))
f.lamb <- function(t,y){y[1]}
f.mu <- function(t,y){y[1]}
lamb_par <- c(0.1)
mu_par <- c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

tree_BCSTDCST <-  fit_bd(fos_tree, tot_time, f.lamb, f.mu, lamb_par, mu_par, f = frac, 
                         cst.lamb=cst.lamb, cst.mu=cst.mu, expo.lamb=expo.lamb, expo.mu=expo.mu, fix.mu=fix.mu, dt = 1e-3,
                         cond="crown")

tree_BCSTDCST$model <- "BCST_DCST"


# FUNCTION ENV_FITR -------------------------------------------------------

env_fitr <- function(i){
  
  #####################################################
  ### Environment dependent (exponential variation) ###
  #####################################################
  
  #	BEnvVar_EXPO
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){0}
  lamb_par<-c(0.1, 0)
  mu_par<-c()
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  BEnvVar_expo <- fit_env(fos_tree,env_data = i, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                          cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                          fix.mu=fix.mu,cond=cond)
  
  ########################
  ###	BEnvVarDCST_EXPO ###
  ########################
  
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){y[1]}
  lamb_par<-c(abs(BEnvVar_expo$lamb_par[1]), BEnvVar_expo$lamb_bar[2]) # from lamb_par[1] and lamb_par[2] of previous model output
  mu_par<-c(0.01)
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
  
  BEnvVarDCST <- fit_env(fos_tree,env_data = i, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac, 
                         cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  
  
  ########################
  ### BCSTDEnvVar_EXPO ###
  ########################
  
  f.lamb<-function(t,x,y){y[1]}
  f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
  lamb_par<-c(tree_BCSTDCST[1])
  mu_par<-c(0.01, 0)
  cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  BCSTDenvVar <- fit_env(fos_tree,env_data = i, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                         cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                         fix.mu=fix.mu,cond=cond)
  
  ###########################
  ### BEnvVarDEnvVar ###
  ###########################
  
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
  lamb_par<-c(abs(BEnvVarDCST$lamb_par[1]), BEnvVarDCST$lamb_par[2])
  mu_par<-c(0.01, 0)
  cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  BEnvVarDEnvVar <- fit_env(fos_tree,env_data = i, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                            cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                            fix.mu=fix.mu,cond=cond)
  
  
  return(list(tree_BCST, tree_BCSTDCST, BEnvVar_expo, BEnvVarDCST, BCSTDenvVar, BEnvVarDEnvVar))
}

result_envfitr <- mclapply(X = env_data_list, FUN = env_fitr, mc.cores = 4)
