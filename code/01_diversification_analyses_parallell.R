# 01_diversification_parallell_fit.R

#Required R packages
library(picante)
library(pspline)
library(RPANDA)
library(phytools); library(gridExtra)
library(parallel) # can do this on WSL2 (if using Windows)

# R codes and functions which need to be sourced
# source("/YOUR_PATH/fit_bd.R")
# source("/YOUR_PATH/fit_env_bd.R")
# source("/YOUR_PATH/likelihood_bd.R")
# source("/YOUR_PATH/Phi.R")
# source("/YOUR_PATH/Psi.R")
# source("/YOUR_PATH/integrate.R")

# Explanation of the models we'll fit:
#	BCST_DCST
#		+ standard birth-death model, with constant speciation and extinction rates
#	BCST_Yule
#		+ standard Yule model, speciation rate is constant and no extinction
#	BAltiVar_EXPO
#		+ speciation rate is correlated exponentially to the rate of altitude variation, with no extinction (mu=0)
#	BAltiVarDCST_EXPO
#		+ speciation rate is correlated exponentially to the rate of altitude variation, with constant extinction
#	BCSTDAltiVar_EXPO
#		+ speciation rate is constant, and extinction rate is correlated exponentially to the rate of altitude variation
#	BAltiVarDAltiVar_EXPO
#		+ speciation and extinction rates are correlated exponentially to the rate of altitude variation
#	BAltiVar_LIN
#		+ + speciation rate is correlated linearly to the rate of altitude variation, with no extinction (mu=0)
#	BAltiVarDCST_LIN
#		+ speciation rate is correlated linearly to the rate of altitude variation, with constant extinction
#	BCSTDAltiVar_LIN
#		+ speciation rate is constant, and extinction rate is correlated linearly to the rate of altitude variation
#	BAltiVarDAltiVar_LIN
#		+ speciation and extinction rates are correlated linearly to the rate of altitude variation

# Load --------------------------------------------------------------------

load(file = "data/intermediate_data/diversification_analyses/env_data_list.RData")
trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

is.ultrametric(trees[[1]])
fos_tree <- phytools::force.ultrametric(trees[[1]],"extend")
fos_tree$edge.length <- fos_tree$edge.length * 100
plot(fos_tree); axisPhylo()
is.ultrametric(fos_tree)

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

#################
## Fit the models ##
#################

env_fitter <- function(env_data_list){
# 
  env_data <- as.data.frame(env_data_list)

  Anilios <- fos_tree

  nb.trees <- 1 # This can be useful to run over a posterior trees distribution: select the number of trees you want to run
  #posteriors <- sample(Anilios, nb.trees) # Uncomment this if you want to run the models over a posterior distribution

  final_Anilios<-NULL
  Anilios_res<-NULL

  for (i in 1: nb.trees) # This can be useful to run over a posterior trees distribution
#     
  env_data <- env_data_list$mean_sco
  
    print(i)
  phylo<-Anilios # If a posterior trees distribution is used, then it would be "posteriors[[i]]"
  tot_time<-max(node.age(phylo)$ages)
  f<-Ntip(phylo)/50 # As I found in Reptile Database
  cond="crown"
  
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
  
  #####################################################
  ###### Env Dependence (exponential variation) ######
  #####################################################
  print(i)
  
  # BEnvVar EXPO
  print("BEnvVar EXPO")
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){0}
  lamb_par<-c(0.1,0)
  mu_par<-c()
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  treei_BEnvVar_EXPO<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVar_EXPO)
  
  
  # BEnvVar DCST EXPO
  print("BEnvVar DCST EXPO")
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){y[1]}
  lamb_par<-c(abs(treei_BEnvVar_EXPO$lamb_par[1]),treei_BEnvVar_EXPO$lamb_par[2])
  mu_par<-c(0.01)
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BEnvVarDCST_EXPO<-fit_env(phylo,as.data.frame(env_data),tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVarDCST_EXPO)
  
  
  # BCST DEnvVar EXPO
  print("BCST DEnvVar EXPO")
  f.lamb<-function(t,x,y){y[1]}
  f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
  lamb_par<-c(treei_BCSTDCST$lamb_par[1])
  mu_par<-c(0.01,0)
  cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BCSTDEnvVar_EXPO<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BCSTDEnvVar_EXPO)
  
  
  # BEnvVar DEnvVar EXPO
  print("BEnvVar DEnvVar EXPO")
  f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
  f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
  lamb_par<-c(abs(treei_BEnvVarDCST_EXPO$lamb_par[1]),treei_BEnvVarDCST_EXPO$lamb_par[2])
  mu_par<-c(0.01,0)
  cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BEnvVarDEnvVar_EXPO<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVarDEnvVar_EXPO)
  
  ################################################
  ###### Env Dependence (linear variation) ######
  ################################################
  print(i)
  
  # BEnvVar LIN
  print("BEnvVar LIN")
  f.lamb<-function(t,x,y){y[1]+y[2]*x}
  f.mu<-function(t,x,y){0}
  lamb_par<-c(abs(treei_BEnvVar_EXPO$lamb_par[1]),0)
  mu_par<-c()
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
  
  treei_BEnvVar_LIN<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVar_LIN)
  
  
  # BEnvVar DCST LIN
  print("BEnvVar DCST LIN")
  f.lamb<-function(t,x,y){y[1]+y[2]*x}
  f.mu<-function(t,x,y){y[1]}
  lamb_par<-c(abs(treei_BEnvVar_LIN$lamb_par[1]),treei_BEnvVar_LIN$lamb_par[2])
  mu_par<-c(0.01)
  cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BEnvVarDCST_LIN<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVarDCST_LIN)
  
  
  # BCST DEnvVar LIN
  print("BCST DEnvVar LIN")
  f.lamb<-function(t,x,y){y[1]}
  f.mu<-function(t,x,y){y[1]+y[2]*x}
  lamb_par<-c(treei_BCSTDCST$lamb_par[1])
  mu_par<-c(0.02,0)
  cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BCSTDEnvVar_LIN<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BCSTDEnvVar_LIN)
  
  
  # BEnvVar DEnvVar LIN
  print("BEnvVar DEnvVar LIN")
  f.lamb<-function(t,x,y){y[1]+y[2]*x}
  f.mu<-function(t,x,y){y[1]+y[2]*x}
  lamb_par<-c(abs(treei_BEnvVarDCST_LIN$lamb_par[1]),treei_BEnvVarDCST_LIN$lamb_par[2])
  mu_par<-c(0.02,0)
  cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
  
  treei_BEnvVarDEnvVar_LIN<-fit_env(phylo,env_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
  print(treei_BEnvVarDEnvVar_LIN)
  
  
  ############# RESULTS ###########################################
  
  results<-matrix(NA,10,9)
  colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta", "AICcWt")
  
  #Models
  results[,1]<-c("BEnvVar_EXPO","BEnvVarDCST_EXPO","BCSTDEnvVar_EXPO","BEnvVarDEnvVar_EXPO",
                 "BEnvVar_LIN","BEnvVarDCST_LIN","BCSTDEnvVar_LIN","BEnvVarDEnvVar_LIN", "BCST_DCST", "BCST_Yule")
  
  #Parameters
  results[1,2]<-2
  results[2,2]<-3
  results[3,2]<-3
  results[4,2]<-4
  results[5,2]<-2
  results[6,2]<-3
  results[7,2]<-3
  results[8,2]<-4
  results[9,2]<-2
  results[10,2]<-1
  
  #logL
  results[1,3]<-treei_BEnvVar_EXPO$LH
  results[2,3]<-treei_BEnvVarDCST_EXPO$LH
  results[3,3]<-treei_BCSTDEnvVar_EXPO$LH
  results[4,3]<-treei_BEnvVarDEnvVar_EXPO$LH
  results[5,3]<-treei_BEnvVar_LIN$LH
  results[6,3]<-treei_BEnvVarDCST_LIN$LH
  results[7,3]<-treei_BCSTDEnvVar_LIN$LH
  results[8,3]<-treei_BEnvVarDEnvVar_LIN$LH
  results[9,3]<-treei_BCSTDCST$LH
  results[10,3]<-treei_BCSTYule$LH
  
  #AICc
  results[1,4]<-treei_BEnvVar_EXPO$aicc
  results[2,4]<-treei_BEnvVarDCST_EXPO$aicc
  results[3,4]<-treei_BCSTDEnvVar_EXPO$aicc
  results[4,4]<-treei_BEnvVarDEnvVar_EXPO$aicc
  results[5,4]<-treei_BEnvVar_LIN$aicc
  results[6,4]<-treei_BEnvVarDCST_LIN$aicc
  results[7,4]<-treei_BCSTDEnvVar_LIN$aicc
  results[8,4]<-treei_BEnvVarDEnvVar_LIN$aicc
  results[9,4]<-treei_BCSTDCST$aicc
  results[10,4]<-treei_BCSTYule$aicc
  
  #Lambda0
  results[1,5]<-abs(treei_BEnvVar_EXPO$lamb_par[1])
  results[2,5]<-abs(treei_BEnvVarDCST_EXPO$lamb_par[1])
  results[3,5]<-abs(treei_BCSTDEnvVar_EXPO$lamb_par[1])
  results[4,5]<-abs(treei_BEnvVarDEnvVar_EXPO$lamb_par[1])
  results[5,5]<-abs(treei_BEnvVar_LIN$lamb_par[1])
  results[6,5]<-abs(treei_BEnvVarDCST_LIN$lamb_par[1])
  results[7,5]<-abs(treei_BCSTDEnvVar_LIN$lamb_par[1])
  results[8,5]<-abs(treei_BEnvVarDEnvVar_LIN$lamb_par[1])
  results[9,5]<-abs(treei_BCSTDCST$lamb_par[1])
  results[10,5]<-abs(treei_BCSTYule$lamb_par[1])
  
  #Alpha Temp
  results[1,6]<-treei_BEnvVar_EXPO$lamb_par[2]
  results[2,6]<-treei_BEnvVarDCST_EXPO$lamb_par[2]
  results[4,6]<-treei_BEnvVarDEnvVar_EXPO$lamb_par[2]
  results[5,6]<-treei_BEnvVar_LIN$lamb_par[2]
  results[6,6]<-treei_BEnvVarDCST_LIN$lamb_par[2]
  results[8,6]<-treei_BEnvVarDEnvVar_LIN$lamb_par[2]	
  
  #Mu0
  results[2,7]<-abs(treei_BEnvVarDCST_EXPO$mu_par[1])
  results[3,7]<-abs(treei_BCSTDEnvVar_EXPO$mu_par[1])
  results[4,7]<-abs(treei_BEnvVarDEnvVar_EXPO$mu_par[1])
  results[6,7]<-abs(treei_BEnvVarDCST_LIN$mu_par[1])
  results[7,7]<-abs(treei_BCSTDEnvVar_LIN$mu_par[1])
  results[8,7]<-abs(treei_BEnvVarDEnvVar_LIN$mu_par[1])
  results[9,6]<-treei_BCSTDCST$mu_par[1]
  
  #Beta Temp
  results[3,8]<-treei_BCSTDEnvVar_EXPO$mu_par[2]
  results[4,8]<-treei_BEnvVarDEnvVar_EXPO$mu_par[2]
  results[7,8]<-treei_BCSTDEnvVar_LIN$mu_par[2]
  results[8,8]<-treei_BEnvVarDEnvVar_LIN$mu_par[2]
  
  #AICc Weights
  results[,9] <- aic.w(as.numeric(results[,4]))
  #ALL RESULTS	
  final_Anilios[[i]]<-results
  
  resi<-list("Clade_age"=tot_time,"Taxon_sampling"=Ntip(phylo),"Sampling_fraction"=f,
             "BEnvVar_EXPO"=treei_BEnvVar_EXPO,"BEnvVarDCST_EXPO"=treei_BEnvVarDCST_EXPO,"BCSTDEnvVar_EXPO"=treei_BCSTDEnvVar_EXPO,"BEnvVarDEnvVar_EXPO"=treei_BEnvVarDEnvVar_EXPO,
             "BEnvVar_LIN"=treei_BEnvVar_LIN,"BEnvVarDCST_LIN"=treei_BEnvVarDCST_LIN,"BCSTDEnvVar_LIN"=treei_BCSTDEnvVar_LIN,"BEnvVarDEnvVar_LIN"=treei_BEnvVarDEnvVar_LIN,
             "BCSTDCST"=treei_BCSTDCST, "BCST_Yule"=treei_BCSTYule)
  
  Anilios_res<-c(Anilios_res,list(resi))
  return(final_Anilios, Anilios_res)
}

result_fit <- mclapply(X = env_data_list, FUN=function(i)env_fitter(i), mc.cores = 6)

# save.image(file = 'output/20220823_multiRPANDA_fit.Rdata')

write.table(final_Anilios,file="/YOUR_PATH/Results_Anilios_EnvDep.txt", quote=FALSE,sep="\t",row.names=FALSE)
save(final_Anilios,file="/YOUR_PATH/final_AniliosTrees_EnvDep.RData")
save(Anilios_res,file="/YOUR_PATH/Anilios_EnvDep.RData")


