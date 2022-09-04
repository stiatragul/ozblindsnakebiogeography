# 02_diversification_analyses.R
# Putter Tiatragul
# June 2022

# The purpose of this script is to fit diversification models with RPANDA
#### Pre-reqs
# source("code/utility/00_diversification_env_data_prep.R")
# source("code/utility/00_diversification_tree_data_prep.R")
####

## IMPORTANT ## 
# These analyses can be slow, especially when running multiple models. 
# From initial testing, using R on WSL2 is faster than R on Windows even when using the same source files on Windows. 
# So when running analyses, I'm using Rstudio server launched from WSL2. 

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
#	BEnvVarDEnvVar_EXPO
#		+ speciation and extinction rates are correlated exponentially to the rate of environmental variation
#	BEnvVar_LIN
#		+ + speciation rate is correlated linearly to the rate of environmental variation, with no extinction (mu=0)
#	BEnvVarDCST_LIN
#		+ speciation rate is correlated linearly to the rate of environmental variation, with constant extinction
#	BCSTDEnvVar_LIN
#		+ speciation rate is constant, and extinction rate is correlated linearly to the rate of environmental variation
#	BEnvVarDEnvVar_LIN
#		+ speciation and extinction rates are correlated linearly to the rate of environmental variation

## Helpful function for tabling results from fitting RPANDA
tab_func <- function(i, .mu, gen_name = NULL){
  
  data.frame(model = paste(names(env_data_list), gen_name),
             aicc = sapply(i, '[[', 'aicc'),
             logL = sapply(i, '[[', 'LH'),
             Lamb0 = sapply(1:length(i),function(x){i[[x]]$lamb_par[1]}),
             alpha = sapply(1:length(i),function(x){i[[x]]$lamb_par[2]})) %>% 
    dplyr::mutate(mu0 = ifelse(.mu == TRUE, sapply(1:length(i),function(x){i[[x]]$mu_par[1]}), NA),
                  beta = ifelse(.mu == TRUE, sapply(1:length(i),function(x){i[[x]]$mu_par[2]}), NA),
                  AICw = phytools::aic.w(as.numeric(aicc))) %>% arrange(aicc)
}


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

#####################################################
### Environment dependent (exponential variation) ###
#####################################################

#	BEnvVar_EXPO
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.1, 0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

# Parallel on 4 cores
# nclust <- parallel::detectCores()
num_cores <- detectCores()

system.time(
  tree_BEnvVar_EXPO_res <- mclapply(env_data_list, 
                                    FUN=function(i){fit_env(fos_tree,env_data = i, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                                                            cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                                                            fix.mu=fix.mu,cond=cond)},
                                    mc.cores = 8)
)

# tab_func(tree_BEnvVar_EXPO_res, .mu = FALSE, gen_name = "BEnvVar_expo")
# plot_fit_env(tree_BEnvVar_EXPO_res$arid_str, env_data_list$arid_str, tot_time)

########################
###	BEnvVarDCST_EXPO ###
########################

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(0.01),0.001) # from lamb_par[1] and lamb_par[2] of previous model output
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

system.time(
  
  tree_BEnvVarDCST_EXPO_res <- mclapply(env_data_list, 
                                        FUN=function(i){fit_env(fos_tree, env_data = i, 
                                                                tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac, 
                                                                cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,
                                                                expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)},
                                        mc.cores = 8)
  
)

# tab_func(tree_BEnvVarDCST_EXPO_res, .mu = TRUE, gen_name = "BEnvVarDCST_expo")
# plot_fit_env(tree_BEnvVarDCST_EXPO_res$mean_sco, env_data_list[1], tot_time)

########################
### BCSTDEnvVar_EXPO ###
########################

f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(0.1)
mu_par<-c(0.01, 0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

system.time(
  tree_BCSTDenvVar_EXPO_res <- mclapply(env_data_list, 
                                        FUN=function(i){fit_env(fos_tree,env_data = i, 
                                                                tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                                                                cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                                                                fix.mu=fix.mu,cond=cond)},
                                        mc.cores = 8)
)

# tab_func(tree_BCSTDenvVar_EXPO_res, .mu = TRUE, gen_name = "BCSTDenvVar_expo")
# plot_fit_env(tree_BCSTDenvVar_EXPO_res$mean_sco, env_data_list[1], tot_time)

###########################
### BEnvVarDEnvVar_EXPO ###
###########################

f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(0.01), 0)
mu_par<-c(0.01, 0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

system.time(
  tree_BEnvVarDEnvVar_EXPO_res <- mclapply(env_data_list, 
                                           FUN=function(i){fit_env(fos_tree,env_data = i, 
                                                                   tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,
                                                                   cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,
                                                                   fix.mu=fix.mu,cond=cond)},
                                           mc.cores = 8)
)


# plot_fit_env(BEnvVarDEnvVar_expo_tmp, env_data_list$mean_sco, tot_time)



# Tabulate the results ----------------------------------------------------

rbind(tab_func(tree_BCSTDenvVar_EXPO_res, .mu = TRUE, gen_name = "BCSTDenvVar_expo"),
      tab_func(tree_BEnvVarDCST_EXPO_res, .mu = TRUE, gen_name = "BEnvVarDCST_expo"),
      tab_func(tree_BEnvVarDEnvVar_EXPO_res, .mu = TRUE, gen_name = "BEnvVarDEnvVar_expo"),
      tab_func(tree_BEnvVar_EXPO_res, .mu = FALSE, gen_name = "BEnvVar_expo")) %>% 
  arrange(desc(aicc))

plot_fit_env(tree_BEnvVar_EXPO_res$arid_val, env_data_list$arid_val, tot_time)

# save.image(file='output/fit_env_output.RData')
# load("output/fit_env_output.RData")


# Plot results ------------------------------------------------------------

plot_fit_env(tree_BEnvVar_EXPO_res$mean_sco, env_data_list[1], tot_time)


# Individual models
# BEnvVarDEnvVar_expo_tmp <- fit_env(fos_tree,env_data = env_data_list$mean_sco, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
# BEnvVarDEnvVar_expo_tmp_straume <- fit_env(fos_tree,env_data = env_data_list$mean_str, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
# BEnvVarDEnvVar_expo_glo <- fit_env(fos_tree,env_data = env_data_list$global_temp, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
# BEnvVarDEnvVar_expo_arid <- fit_env(fos_tree,env_data = env_data_list$arid_sco, tot_time,f.lamb,f.mu,lamb_par,mu_par,f=frac,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)




################# ARCHIVE


# time dependent ----------------------------------------------------------


# # Exponential
# f.lamb_exp <- function(t,y){y[1]*exp(y[2]*t)}
# f.mu_exp <- function(t,y){y[1]*exp(y[2]*t)}
# 
# # Linear
# f.lamb_lin <- function(t,y){abs(y[1]+y[2]*t)}
# f.mu_lin <- function(t,y){y[1]}
# 
# # Constant
# f.lamb_cst <- function(t,y){y[1]}
# f.mu_cst <- function(t,y){y[1]}
# 
# # No extinction
# f.mu_null <- function(t,y){0}
# 
# # initial speciation (lambda) and extinction (mu)
# lamb_par_init_exp <- c(0.6,0.001)
# mu_par_init_exp <- c(0.005, 0.001)
# 
# lamb_par_init_lin <- c(0.6, 0.001)
# mu_par_init_lin <- c(0.005, 0.001)
# 
# mu_par_init_cst <- c(0.005)
# 
# # No extinction
# mu_par_init_null <- c()
# 
# 
# # Total time
# tot_time <- max(node.age(fos_tree)$ages)
# 
# # Plot Fitted speciation rate
# ## Here we tweak out the function `RPANDA::plot_fit_bd` to one plotting device.
# plot_fit_bd_res <- function(fit_bd_res, tot_time){
#   
#   fit_bd_res <- fit_bd_res
#   t <- seq(0, tot_time, length.out = 100)
#   r <- function(t) {fit_bd_res$f.lamb(t) - fit_bd_res$f.mu(t)}
#   
#   par(mfrow = c(2,2))
#   plot(-t, fit_bd_res$f.lamb(t), type = "l", xlab = "Time", ylab = "Speciation rate", main = "Fitted speciation rate")
#   plot(-t, fit_bd_res$f.mu(t), type = "l", xlab = "Time", ylab = "Extinction rate", main = "Fitted extinction rate")
#   plot(-t, r(t), type = "l", xlab = "Time", ylab = "Net diversification rate", main = "Fitted net diversification rate")
#   
# }
# 
# # dev.off()
# # plot_fit_bd_res(result_lin, tot_time)
# 
# 
# # Diversity through time
# plot_dtt(result_lin, tot_time, N0 = 48)


