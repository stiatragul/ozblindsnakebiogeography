#02_table_fit_env_results.R
# September 2022

# This script will organise the results from the parallel fit and the individual fit together in one table

# From 01_diversification_analyses_parallel.R I successfully fitted several environment-dependent birth-death models to
# min_sco ;arid_sco; mean_str; min_str;  arid_str;  min_val; max_val; arid_val. 
# The following had to be fitted individually using initial lamb_par 0.1 and/or 0. 
# mean_val; # max_str; mean_sco; max_sco; global temp


# libraries ---------------------------------------------------------------

library(RPANDA); library(dplyr); library(phytools)

# Load results ------------------------------------------------------------
# from parallel fit and individual fit
load('output/20220823_multiRPANDA_fit.Rdata')
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_global_temp.Rdata'); global_temp_res <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_max_sco.Rdata'); max_sco_res <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_max_str.Rdata'); max_str_res <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_mean_sco.Rdata'); mean_sco_res <- Anilios_res

# Data --------------------------------------------------------------------
load(file = "data/intermediate_data/diversification_analyses/env_data_list.RData")
trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

# Setting up analyses -----------------------------------------------------
fos_tree <- phytools::force.ultrametric(trees[[1]],"extend")
fos_tree$edge.length <- fos_tree$edge.length * 100
phylo<-fos_tree # If a posterior trees distribution is used, then it would be "posteriors[[i]]"
tot_time<-max(node.age(phylo)$ages)
f<-Ntip(phylo)/50 # As I found in Reptile Database
cond="crown"

# Fit models ---------------------------------------------------------------

# Crown age
tot_time <- max(node.age(fos_tree)$ages)
cond = "crown"

# Limit year to max total time
subset_env_age <- function(df) {
  df <- df %>% filter(time_mya < 26)
}

env_data_list <- lapply(env_data_list, subset_env_age)

## Fraction of diversity represented. In our tree we have 
## 48 Described species of Anilios + 2 more grypus undescribed + 1 more ligatus - A. splendidus
## So our "true" species at the moment is 50 and we have 41 tips in our tree

length(fos_tree$tip.label)
frac = 38/50

# BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(0.1)
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDCST<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDCST)

# BCST Yule (constant Birth (Yule))
print("BCST Yule")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(0.1)
mu_par<-c(0)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BCSTYule<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTYule)


# Organise table ----------------------------------------------------------

tablr <- function(env_res, name_it){
  results<-matrix(NA,10,9)
  colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta", "AICcWt")
  #Models
  results[,1]<-c("BEnvVar_EXPO","BEnvVarDCST_EXPO","BCSTDEnvVar_EXPO","BEnvVarDEnvVar_EXPO",
                 "BEnvVar_LIN","BEnvVarDCST_LIN","BCSTDEnvVar_LIN","BEnvVarDEnvVar_LIN", "BCST_DCST", "BCST_Yule")
  # BCST DCST (constant Birth-death)
  print("BCST DCST")
  f.lamb<-function(x,y){y[1]}
  f.mu<-function(x,y){y[1]}
  lamb_par<-c(0.1)
  mu_par<-c(0.01)
  cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
  treei_BCSTDCST<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  # print(treei_BCSTDCST)
  # BCST Yule (constant Birth (Yule))
  f.lamb<-function(x,y){y[1]}
  f.mu<-function(x,y){y[1]}
  lamb_par<-c(0.1)
  mu_par<-c(0)
  cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  treei_BCSTYule<-fit_bd(phylo,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  # print(treei_BCSTYule)
  #Parameters
  results[1,2]<-2; results[2,2]<-3; results[3,2]<-3
  results[4,2]<-4; results[5,2]<-2; results[6,2]<-3
  results[7,2]<-3; results[8,2]<-4; results[9,2]<-2
  results[10,2]<-1
  #logL
  results[1,3]<-env_res$BEnvVar_EXPO$LH; results[2,3]<-env_res$BEnvVarDCST_EXPO$LH
  results[3,3]<-env_res$BCSTDEnvVar_EXPO$LH; results[4,3]<-env_res$BEnvVarDEnvVar_EXPO$LH
  results[5,3]<-env_res$BEnvVar_LIN$LH; results[6,3]<-env_res$BEnvVarDCST_LIN$LH
  results[7,3]<-env_res$BCSTDEnvVar_LIN$LH; results[8,3]<-env_res$BEnvVarDEnvVar_LIN$LH
  results[9,3]<-env_res$BCSTDCST$LH; 
  results[10,3]<-treei_BCSTYule$LH
  #AICc
  results[1,4]<-env_res$BEnvVar_EXPO$aicc; results[2,4]<-env_res$BEnvVarDCST_EXPO$aicc
  results[3,4]<-env_res$BCSTDEnvVar_EXPO$aicc; results[4,4]<-env_res$BEnvVarDEnvVar_EXPO$aicc
  results[5,4]<-env_res$BEnvVar_LIN$aicc; results[6,4]<-env_res$BEnvVarDCST_LIN$aicc
  results[7,4]<-env_res$BCSTDEnvVar_LIN$aicc; results[8,4]<-env_res$BEnvVarDEnvVar_LIN$aicc
  results[9,4]<-env_res$BCSTDCST$aicc; 
  results[10,4]<-treei_BCSTYule$aicc
  #Lambda0
  results[1,5]<-abs(env_res$BEnvVar_EXPO$lamb_par[1]); results[2,5]<-abs(env_res$BEnvVarDCST_EXPO$lamb_par[1])
  results[3,5]<-abs(env_res$BCSTDEnvVar_EXPO$lamb_par[1]); results[4,5]<-abs(env_res$BEnvVarDEnvVar_EXPO$lamb_par[1])
  results[5,5]<-abs(env_res$BEnvVar_LIN$lamb_par[1]); results[6,5]<-abs(env_res$BEnvVarDCST_LIN$lamb_par[1])
  results[7,5]<-abs(env_res$BCSTDEnvVar_LIN$lamb_par[1]); results[8,5]<-abs(env_res$BEnvVarDEnvVar_LIN$lamb_par[1])
  results[9,5]<-abs(env_res$BCSTDCST$lamb_par[1]); 
  results[10,5]<-abs(treei_BCSTYule$lamb_par[1])
  #Alpha Temp
  results[1,6]<-env_res$BEnvVar_EXPO$lamb_par[2]; results[2,6]<-env_res$BEnvVarDCST_EXPO$lamb_par[2]
  results[4,6]<-env_res$BEnvVarDEnvVar_EXPO$lamb_par[2]; results[5,6]<-env_res$BEnvVar_LIN$lamb_par[2]
  results[6,6]<-env_res$BEnvVarDCST_LIN$lamb_par[2]; results[8,6]<-env_res$BEnvVarDEnvVar_LIN$lamb_par[2]	
  #Mu0
  results[2,7]<-abs(env_res$BEnvVarDCST_EXPO$mu_par[1]); results[3,7]<-abs(env_res$BCSTDEnvVar_EXPO$mu_par[1])
  results[4,7]<-abs(env_res$BEnvVarDEnvVar_EXPO$mu_par[1]); results[6,7]<-abs(env_res$BEnvVarDCST_LIN$mu_par[1])
  results[7,7]<-abs(env_res$BCSTDEnvVar_LIN$mu_par[1]); results[8,7]<-abs(env_res$BEnvVarDEnvVar_LIN$mu_par[1])
  results[9,6]<-env_res$BCSTDCST$mu_par[1]
  #Beta Temp
  results[3,8]<-env_res$BCSTDEnvVar_EXPO$mu_par[2]; results[4,8]<-env_res$BEnvVarDEnvVar_EXPO$mu_par[2]
  results[7,8]<-env_res$BCSTDEnvVar_LIN$mu_par[2]; results[8,8]<-env_res$BEnvVarDEnvVar_LIN$mu_par[2]
  #AICc Weights
  results[,9] <- aic.w(as.numeric(results[,4])); #ALL RESULTS
  final_Anilios <- as.data.frame(results)
  final_Anilios[, 2:9] <- lapply(final_Anilios[, 2:9], as.numeric)
  final_Anilios[, 2:4] <- lapply(final_Anilios[, 2:4], round, 2)
  final_Anilios$AICcWt <- round(final_Anilios$AICcWt, 2)
  final_Anilios[5:8] <- lapply(final_Anilios[5:8], formatC, format = "e", digits = 2)
  final_Anilios$env <- name_it
  return(final_Anilios)
}

# Extract results in data frame
tab1 <- tablr(result_fit$min_sco[[1]], "min_sco")
tab2 <- tablr(result_fit$arid_sco[[1]], "arid_sco")
tab3 <- tablr(result_fit$mean_str[[1]], "mean_str")
tab4 <- tablr(result_fit$min_str[[1]], "min_str")
tab5 <- tablr(result_fit$arid_str[[1]], "arid_str")
tab6 <- tablr(result_fit$min_val[[1]], "min_val")
tab7 <- tablr(result_fit$max_val[[1]], "max_val")
tab8 <- tablr(result_fit$arid_val[[1]], "arid_val")
tab9 <- tablr(global_temp_res[[1]], "global_temp")
tab10 <- tablr(mean_sco_res[[1]], "mean_sco")
tab11 <- tablr(max_sco_res[[1]], "max_sco")
tab12 <- tablr(max_str_res[[1]], "max_str")


# Combine the table
rbind(tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8,
      tab9, tab10, tab11, tab12) %>% 
  dplyr::arrange(AICc)
