# 04_BioGeoBEARS_results_bsm.R
# August 2022
# Code adapted from Carlos Pavon Vazquez

# Stochastic mapping

# load results from analysis script

# Libraries ---------------------------------------------------------------

library(GenSA);library(FD)      
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)
# library(dplyr)
library(parallel); library(snow)

# Tree --------------------------------------------------------------------

tree.c <- ape::read.tree('data/intermediate_data/bears/blindsnake_b.tre')
tree.c$edge.length <- tree.c$edge.length * 100

# Saved results -----------------------------------------------------------

load("data/intermediate_data/bears/bears_model_fit.Rdata")

############################
#####Stochastic mapping#####
############################

load("data/intermediate_data/bears/AniliosDEC_JX_full.Rdata")
resjx
tr <- tree.c
# BSM = Biogeographic Stochastic Mapping for the best fitting model

model_name = "DECjx"
res = resjx

pdffn = paste0("output/","Anilios_", model_name, "_v1.pdf")
pdf(pdffn, width=6, height=6)

analysis_titletxt = paste0(model_name, " on Anilios")

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j", "x"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, 
                                plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j", "x"), plotwhat="pie", label.offset=0.45, 
                         tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                         include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it

#######################################################
# Stochastic mapping on DEC+J+x
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = "data/intermediate_data/bears/BSM_inputs_file.Rdata"
runInputsSlow = FALSE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...set runslow == FALSE
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
# max_range_size = 5

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0("output/", model_name, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, width=6, height=6)

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j", "x"), 
                         label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE,
                         colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, 
                              colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), 
                              root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = stratified

# Loop through the maps and plot to PDF
pdffn = paste0("output/", model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=11.33, height=7.5)

nummaps_goal = 50
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j", "x"),
           label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, 
           include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times


clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables,
                                     areanames, actual_names, timeperiod = c(0,0.1))

count_ana_clado_events(clado_events_tables, ana_events_tables,
                       areanames, actual_names)

counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables,
                                     areanames, actual_names, timeperiod = c(5,10))



#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0('data/intermediate_data/bears/',tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)


# SAVE

save.image(file = 'data/intermediate_data/bears/BSB_50.Rdata')

