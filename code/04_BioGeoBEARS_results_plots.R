# 04_BioGeoBEARS_results_plots.R
# August 2022

# Test best model fit from BioGeoBEARS analyses

# load results from analysis script
# source("code/03_BioGeoBEARS_analyses_parallel.R")


# Libraries ---------------------------------------------------------------

library(GenSA);library(FD)      
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)
# library(parallel); library(snow)


# Tree --------------------------------------------------------------------

tree.c <- ape::read.tree('data/intermediate_data/bears/blindsnake_b.tre')


# Saved results -----------------------------------------------------------
## These saved results are from source('code/03_BioGeoBEARS_analyses_parallel.R)

bio_geo_bears_files_list <- list.files(path = "./data/intermediate_data/bears/",
                                       pattern = "*.Rdata", full.names = T)

# Load individual saved files
lapply(bio_geo_bears_files_list, load, .GlobalEnv)


# Results table -----------------------------------------------------------
# Extract parameters and make results table to compare AICc scores

restable = NULL
teststable = NULL

fin.dec = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.decj = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.decx = extract_params_from_BioGeoBEARS_results_object(results_object=resDECx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.decjx = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.diva = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.divaj = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.divax = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.divajx = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.bayes = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.bayesj = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.bayesx = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.bayesjx = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

restable <- rbind(fin.dec,fin.decj,fin.decx,fin.decjx,fin.diva,fin.divaj,fin.divax,fin.divajx,fin.bayes,fin.bayesj,fin.bayesx,fin.bayesjx)
row.names(restable) = c("fin.dec","fin.decj","fin.decx","fin.decjx","fin.diva","fin.divaj","fin.divax","fin.divajx","fin.bayes","fin.bayesj","fin.bayesx","fin.bayesjx")
restable = BioGeoBEARS::put_jcol_after_ecol(restable)

restable$AICc <- calc_AICc_vals(LnL_vals=restable$LnL,nparam_vals=restable$numparams,samplesize=length(tree.c$tip.label))

# Save this result table
# save(restable, file="data/intermediate_data/bears/restable.Rdata")

# Load result table
# load("data/intermediate_data/bears/restable.Rdata")

aicc.vals <- restable$AICc
names(aicc.vals) <- rownames(restable)

(aicc.pref <- names(sort(aicc.vals))[1])

# BioGeoBEARS::lrttest likelihood ratio test
lrt.fin <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                   LnL_2= restable$LnL[which(rownames(restable)=="fin.decj")],
                   numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                   numparams2= restable$numparams[which(rownames(restable)=="fin.decj")],
                   returnwhat = "all")

lrt.fin2 <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                    LnL_2= restable$LnL[which(rownames(restable)=="fin.decx")],
                    numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                    numparams2= restable$numparams[which(rownames(restable)=="fin.decx")],
                    returnwhat = "all")

lrt.fin3 <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                    LnL_2= restable$LnL[which(rownames(restable)=="fin.dec")],
                    numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                    numparams2= restable$numparams[which(rownames(restable)=="fin.dec")],
                    returnwhat = "all")

lrt.fin
lrt.fin2
lrt.fin3
save(lrt.fin,file="data/intermediate_data/bears/LRT_res1.Rdata")
save(lrt.fin2,file="data/intermediate_data/bears/LRT_res2.Rdata")
save(lrt.fin3,file="data/intermediate_data/bears/LRT_res3.Rdata")
# load("LRT_res1.Rdata")
# load("LRT_res2.Rdata")
# load("LRT_res3.Rdata")



# Plots -------------------------------------------------------------------

# Choose best fitting model to plot

# load("VaranusDEC_JX.Rdata")

# Modify to match best fitting model
analysis_titletxt ="BEST FITTING MODEL" 

# Setup
results_object = resjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
tr=tree.c

# States
resstates = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j","x"), 
                                     plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, 
                                     splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                     include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j","x"), plotwhat="pie", 
                         label.offset=0.45, tipcex=0.02, statecex=0.4, splitcex=0.4, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, plotlegend = T,
                         legend_cex = 0.5)


############################
#####Stochastic mapping#####
############################

# BSM = Biogeographic Stochastic Mapping

# 
# clado_events_tables = NULL
# ana_events_tables = NULL
# lnum = 0
# 
# BSM_inputs_fn = "BSM_inputs_file.Rdata"
# runInputsSlow = TRUE
# 
# if (runInputsSlow){
#   stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=resjx)
#   save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
# } else {
#   # Loads to "stochastic_mapping_inputs_list"
#   load(BSM_inputs_fn)
# } # END if (runInputsSlow)
# names(stochastic_mapping_inputs_list)
# set.seed(seed=as.numeric(Sys.time()))
# 
# runBSMslow = TRUE
# 
# if (runBSMslow == TRUE){
#   # Saves to: RES_clado_events_tables.Rdata
#   # Saves to: RES_ana_events_tables.Rdata
#   BSM_output = runBSM(resjx, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, 
#                       maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, 
#                       save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
#   
#   RES_clado_events_tables = BSM_output$RES_clado_events_tables
#   RES_ana_events_tables = BSM_output$RES_ana_events_tables
# } else {
#   # Load previously saved...
#   
#   # Loads to: RES_clado_events_tables
#   load(file="RES_clado_events_tables.Rdata")
#   # Loads to: RES_ana_events_tables
#   load(file="RES_ana_events_tables.Rdata")
#   BSM_output = NULL
#   BSM_output$RES_clado_events_tables = RES_clado_events_tables
#   BSM_output$RES_ana_events_tables = RES_ana_events_tables
# } # END if (runBSMslow == TRUE)
# 
# clado_events_tables = BSM_output$RES_clado_events_tables
# ana_events_tables = BSM_output$RES_ana_events_tables
# 
# #sample(1:50,1)
# 
# randomly.selected.map <- 7
# 
# sel.clado <- clado_events_tables[[randomly.selected.map]]
# sel.ana <- ana_events_tables[[randomly.selected.map]]
# 
# save(sel.clado,file="clado_events_table.Rdata")
# save(sel.ana,file="ana_events_table.Rdata")
