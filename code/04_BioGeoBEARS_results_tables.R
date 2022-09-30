# 04_BioGeoBEARS_results_table.R
# Sept 2022

# Format best table and likelihood ratio tests
# load results from analysis script
# source("code/03_BioGeoBEARS_analyses_parallel.R") or load from saved intermediate data same thing.

# Libraries ---------------------------------------------------------------

library(GenSA);library(FD)      
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)
library(dplyr)
# library(parallel); library(snow)

# Tree --------------------------------------------------------------------
tree.c <- ape::read.tree('data/intermediate_data/bears/blindsnake_b.tre')
tree.c$edge.length <- tree.c$edge.length * 100

# Saved results -----------------------------------------------------------
load("data/intermediate_data/bears/bears_model_fit.Rdata")

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

# combine these objects into a table
restable <- rbind(fin.dec,fin.decj,fin.decx,fin.decjx,fin.diva,fin.divaj,fin.divax,fin.divajx,fin.bayes,fin.bayesj,fin.bayesx,fin.bayesjx)

# Change row names to help identify what model
row.names(restable) = c("DEC","DEC+J","DEC+X","DEC+J+X",
                        "DIVALIKE","DIVALIKE+J","DIVALIKE+X","DIVALIKE+J+X",
                        "BAYAREALIKE","BAYAREALIKE+J","BAYAREALIKE+X","BAYAREALIKE+J+X")
restable = BioGeoBEARS::put_jcol_after_ecol(restable)

# Save this result table
# write.csv(restable, file="output/supp_biogeobears_table.csv")
# save(restable, file="data/intermediate_data/bears/restable.Rdata")

# Load result table
# load("data/intermediate_data/bears/restable.Rdata")
restable2 = restable

# With AICs:
AICtable = calc_AIC_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams)
restable = cbind(restable, AICtable)
restable_AIC_rellike = AkaikeWeights_on_summary_table(restable=restable, colname_to_use="AIC")
restable_AIC_rellike = put_jcol_after_ecol(restable_AIC_rellike)
restable_AIC_rellike

# With AICcs -- factors in sample size
samplesize = length(tr$tip.label)
AICtable = calc_AICc_column(LnL_vals=restable$LnL, nparam_vals=restable$numparams, samplesize=samplesize)
restable2 = cbind(restable2, AICtable)
restable_AICc_rellike = AkaikeWeights_on_summary_table(restable=restable2, colname_to_use="AICc")
restable_AICc_rellike = put_jcol_after_ecol(restable_AICc_rellike)
restable_AICc_rellike

## VIEW TABLE ##

sorted_resttable <- restable_AICc_rellike %>% dplyr::arrange(AICc) 

## d = parameter for the rate of anagenetic dispersal (range-expansion)
## e = parameter for the rate of anagenetic extinction (range-contraction)
## j = parameter for the rate of cladogenetic long-distance dispersal (founder speciation)

write.csv(sorted_resttable, file="output/supp_biogeobears_table_weights.csv")

# Likelihood ratio tests (LRT) --------------------------------------------

# BioGeoBEARS::lrttest likelihood ratio test to compare top models
# Top five

# Is DEC+j+x significantly better than DECj
lrt.fin_base <- BioGeoBEARS::lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                                LnL_2= restable$LnL[which(rownames(restable)=="fin.dec")],
                                numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                                numparams2= restable$numparams[which(rownames(restable)=="fin.dec")],
                                returnwhat = "all")

lrt.fin_j <- BioGeoBEARS::lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                                LnL_2= restable$LnL[which(rownames(restable)=="fin.decj")],
                                numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                                numparams2= restable$numparams[which(rownames(restable)=="fin.dec")],
                                returnwhat = "all")

lrt.fin_x <- BioGeoBEARS::lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                                LnL_2= restable$LnL[which(rownames(restable)=="fin.decx")],
                                numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                                numparams2= restable$numparams[which(rownames(restable)=="fin.dec")],
                                returnwhat = "all")



# Is DEC+j+x significantly better than DIVALIKE+j+x
lrt.fin2 <- BioGeoBEARS::lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                                 LnL_2= restable$LnL[which(rownames(restable)=="fin.divajx")],
                                 numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                                 numparams2= restable$numparams[which(rownames(restable)=="fin.divajx")],
                                 returnwhat = "all")

# Is DEC+j+x significantly better than BAYAREALIKE+j+x
lrt.fin3 <- BioGeoBEARS::lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decjx")],
                                 LnL_2= restable$LnL[which(rownames(restable)=="fin.bayesjx")],
                                 numparams1= restable$numparams[which(rownames(restable)=="fin.decjx")],
                                 numparams2= restable$numparams[which(rownames(restable)=="fin.bayesjx")],
                                 returnwhat = "all")

# Fit in the DEC+j+x is a statistically significant improvement over DIVALIKE+j+x and BAYAREALIKE+j+x.
lrt.fin
lrt.fin2
lrt.fin3

save(lrt.fin,file="data/intermediate_data/bears/LRT_res1.Rdata")
save(lrt.fin2,file="data/intermediate_data/bears/LRT_res2.Rdata")
save(lrt.fin3,file="data/intermediate_data/bears/LRT_res3.Rdata")
# load("LRT_res1.Rdata")
# load("LRT_res2.Rdata")
# load("LRT_res3.Rdata")


#######################################################
# Plot ancestral states - DEC+j+x
#######################################################
resjx$inputs$states_list

colour_lists <- data.frame(biome = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                           colours = c("#dcc674", #0
                                       "#1f4733", #1
                                       "#e03f28", #2
                                       "#4f937f", #3
                                       "#179d48", #4 
                                       "#5265ae", #5
                                       "#d2e6ae", #6
                                       "#92d1bd", #7
                                       "#eaf400"))  #8 

# mix colours in meyerweb.com color blender

colors_list_for_states  = c("gray87","#dcc674","#1f4733","#e03f28","#4f937f",
                            "#dd9f47","#5265ae","#d2e6ae",
                            "#92d1bd","#d2e6ae", 
                            "#758151", #A and B
                            "#DE8951", #A and C
                            "#8FAA7A", #A and D
                            "#DDB15B", #A and E
                            "#848899", #A and F
                            "#88432D", #B and C
                            "#1B763E", #B and D
                            "#1B763E", #B and E
                            "#3B5776", #B and F
                            "#916D57", #C and D
                            "#935471", #C and F
                            "#9E6550", #D and C
                            "#517A99", #D and F
                            "#gray87", #E and A
                            "#DFEE4F", #G and I
                            "#BAE167", #H and I
                            rep("gray87", 46-26))


# Modify to match best fitting model
analysis_titletxt ="DEC+j+x results"

# Setup
results_object = resjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
tr=tree.c

dev.off()
# States
resstates = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45,
                                     tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                     include_null_range=TRUE, tr=tr,
                                     colors_list_for_states=colors_list_for_states,
                                     tipranges=tipranges)

dev.off()

# pdf('output/pie_ancestral.pdf')
# Pie chart
BioGeoBEARS::plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j","x"), plotwhat="pie", 
                                      label.offset=0.45, tipcex=0.02, statecex=0.4, splitcex=0.4, titlecex=0.8, plotsplits=F, 
                                      cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, plotlegend = T,
                                      # colors_list_for_states=colors_list_for_states,
                                      legend_cex = 0.5)
tiplabels(tr$tip.label, adj = c(0,0.5), bg = NULL, col = NULL, frame = "none")

dev.off()


