# 00_anilios_BioGeoBEARS.R

# base code from Octavio

# Libraries ---------------------------------------------------------------

library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
# library(parallel)
library(snow)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

tree.c <- ape::read.tree('data/tree/anilios_newick_st.tre')

plot(tree.c)

trfn <- "data/tree/anilios_newick_st.tre"
geogfn <- "data/bears_txt/geofile.txt"
# adj <- "data/bears_txt/area_adj.txt"

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 6

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

coln <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

#                     A   B   C   D   E   F   G   H   I
adjacency <- matrix(c(1,	1,	1,	1,	1,	1,	1,	0,	1,  # A
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # B
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # C
                      1,	1,	1,	1,	0,	1,	0,	0,	0,  # D
                      1,	1,	0,	0,	1,	1,	0,	0,	1,  # E
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # F
                      1,	0,	0,	0,	0,	0,	1,	0,	1,  # G
                      0,	0,	0,	0,	0,	0,	0,	1,	1,  # H
                      0,	0,	0,	0,	0,	0,	1,	1,	1), # I
                      ncol=9, nrow=9, byrow=T)


maxAreas = 6
states_list <- rcpp_areas_list_to_states_list(areas=coln, maxareas=maxAreas, include_null_range=TRUE)

# Set up list of TRUEs
states_to_keep <- rep(TRUE, times=length(states_list))

# Delete states that are not adjacent in the adjacency matrix
# Skip the first state, if it's null
if ( is.na(states_list[[1]]) || states_list[[1]] == "_" ){
  start_i = 2
} else {
  start_i = 1
}

for (i in start_i:length(states_list)){
  areas_in_this_state_range_0based = states_list[[i]]
  areas_in_this_state_range_1based = areas_in_this_state_range_0based + 1
  adjacent_if_1s = adjacency[areas_in_this_state_range_1based, areas_in_this_state_range_1based]
  
  print(adjacent_if_1s)
  
  if ( length(adjacent_if_1s) == sum(adjacent_if_1s) ){
    states_to_keep[i] = TRUE
  } else {
    states_to_keep[i] = FALSE
  }
}

# Modify the states_list
states_list <- states_list[states_to_keep]
# https://groups.google.com/d/msg/biogeobears/zkHP_YpA9kY/8c3vR-gIqP0J
statenames <- areas_list_to_states_list_new(areas=coln, maxareas=maxAreas, include_null_range=TRUE, split_ABC=FALSE)
statenames <- statenames[states_to_keep]
unlist(statenames)


############
#####DEC####
############

# NOTE: 'DEC' is identical with Lagrange DEC model.

DEC <- define_BioGeoBEARS_run()
DEC$trfn = trfn
DEC$geogfn = geogfn
# DEC$areas_adjacency_fn = adj
DEC$max_range_size = max_range_size
DEC$min_branchlength = 0.000001; DEC$include_null_range = TRUE
DEC$on_NaN_error = -1e50; DEC$speedup = TRUE; DEC$use_optimx = "GenSA"
DEC$num_cores_to_use = 4

DEC = readfiles_BioGeoBEARS_run(DEC)

DEC$return_condlikes_table = TRUE
DEC$calc_TTL_loglike_from_condlikes_table = TRUE
DEC$calc_ancprobs = TRUE    # get ancestral states from optim run
DEC$states_list = states_list # Need to define manually (using code from Octavio)

check_BioGeoBEARS_run(DEC)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "AniliosDEC.Rdata"
if (runslow)
{
  res = bears_optim_run(DEC)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}


###############
#####DEC+J#####
###############
# D = dispersal
# E = Extinction
# C = cladogenesis
# +
# J is for 'jump'

DECj = define_BioGeoBEARS_run()
DECj$trfn = trfn
DECj$geogfn = geogfn
DECj$max_range_size = max_range_size
DECj$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DECj$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

DECj$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
DECj$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
DECj$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
DECj$num_cores_to_use = 4
DECj$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
DECj$states_list = states_list # Need to define manually

### For areas allowed and time periods ###
# DECj$areas_allowed_fn = area.all
# DECj$timesfn = tiempos

DECj = readfiles_BioGeoBEARS_run(DECj)
# DECj = section_the_tree(inputs=DECj, make_master_table=TRUE, plot_pieces=FALSE)

DECj$return_condlikes_table = TRUE
DECj$calc_TTL_loglike_from_condlikes_table = TRUE
DECj$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
DECj$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DECj$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DECj$BioGeoBEARS_model_object@params_table["e","init"] = estart
DECj$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
DECj$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DECj$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DECj$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(DECj)

resfnj = "AniliosDEC_J.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  resj = bears_optim_run(DECj)
  resj    
  
  save(resj, file=resfnj)
  
  resDECj = resj
} else {
  # Loads to "res"
  load(resfnj)
  resDECj = resj
}

#######################################################
# PDF plots
#######################################################
pdffn = "Anilios_DEC.pdf"
pdf(pdffn, width=8, height=11)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Anilios"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE,
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Anilios"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text",
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", 
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it



######################################
### DIVALIKE AND DIVALIKE+J ANALYSIS ##
#######################################
# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
#
# DIVALIKE is a likelihood interpretation of parsimony
# DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.

#######################################################
# Run DIVALIKE
#######################################################
DIVA = define_BioGeoBEARS_run()
DIVA$trfn = trfn
DIVA$geogfn = geogfn
DIVA$max_range_size = max_range_size
DIVA$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVA$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

DIVA$on_NaN_error = -1e50    
DIVA$speedup = TRUE         
DIVA$use_optimx = TRUE    
DIVA$num_cores_to_use = 4
DIVA$force_sparse = FALSE   
DIVA$states_list = states_list # Need to define manually

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
DIVA = readfiles_BioGeoBEARS_run(DIVA)

DIVA$return_condlikes_table = TRUE
DIVA$calc_TTL_loglike_from_condlikes_table = TRUE
DIVA$calc_ancprobs = TRUE    # get ancestral states from optim run


# Set up DIVALIKE model
# Remove subset-sympatry
DIVA$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVA$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVA$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVA$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVA$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVA$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVA$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
DIVA$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVA$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVA$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

check_BioGeoBEARS_run(DIVA)

runslow = TRUE
resfn = "AniliosDIVA.Rdata"
if (runslow)
{
  res = bears_optim_run(DIVA)
  res    
  
  save(res, file=resfn)
  resDIVALIKE = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKE = res
}


#######################################################
# Run DIVALIKE + J
#######################################################
DIVAj = define_BioGeoBEARS_run()
DIVAj$trfn = trfn
DIVAj$geogfn = geogfn
DIVAj$max_range_size = max_range_size
DIVAj$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVAj$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

DIVAj$on_NaN_error = -1e50    
DIVAj$speedup = TRUE         
DIVAj$use_optimx = TRUE    
DIVAj$num_cores_to_use = 4
DIVAj$force_sparse = FALSE   
DIVAj$states_list = states_list # Need to define manually

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
DIVAj = readfiles_BioGeoBEARS_run(DIVAj)

DIVAj$return_condlikes_table = TRUE
DIVAj$calc_TTL_loglike_from_condlikes_table = TRUE
DIVAj$calc_ancprobs = TRUE    # get ancestral states from optim run


# Set up DIVALIKE+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001
# Input starting values for d, e
DIVAj$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DIVAj$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DIVAj$BioGeoBEARS_model_object@params_table["e","init"] = estart
DIVAj$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Remove subset-sympatry
DIVAj$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVAj$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVAj$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVAj$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVAj$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVAj$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVAj$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

# Allow classic, widespread vicariance; all events equiprobable
DIVAj$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVAj$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVAj$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add jump dispersal/founder-event speciation
DIVAj$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DIVAj$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DIVAj$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under DIVALIKE+J, the max of "j" should be 2, not 3 (as is default in DEC+J)
DIVAj$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
DIVAj$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(DIVAj)

resfn = "AniliosDIVA_J.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(DIVAj)
  res    
  
  save(res, file=resfn)
  
  resDIVALIKEj = res
} else {
  # Loads to "res"
  load(resfn)
  resDIVALIKEj = res
}

pdffn = "Anilios_DIVA_DIVAj.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DIVALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE on Anilios"

# Setup
results_object = resDIVALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

#######################################################
# Plot ancestral states - DIVALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS DIVALIKE+J on Anilios"

# Setup
results_object = resDIVALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#####################
#####BAYAREALIKE#####
#####################

# Likelihood interpretation of Bayesian BayArea. Test importance of cladogenesis.

BAYES = define_BioGeoBEARS_run()
BAYES$trfn = trfn
BAYES$geogfn = geogfn
BAYES$max_range_size = max_range_size
BAYES$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYES$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# BAYES$areas_allowed_fn = area.all
# BAYES$timesfn = tiempos

# Speed options and multicore processing if desired
BAYES$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BAYES$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BAYES$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BAYES$num_cores_to_use = 4
BAYES$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BAYES$states_list = states_list 

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BAYES = readfiles_BioGeoBEARS_run(BAYES)
# BAYES = section_the_tree(inputs=BAYES, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BAYES$return_condlikes_table = TRUE
BAYES$calc_TTL_loglike_from_condlikes_table = TRUE
BAYES$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BAYES$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BAYES$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BAYES$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BAYES$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BAYES$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BAYES$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BAYES$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BAYES$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BAYES$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BAYES$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BAYES$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BAYES$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Check the inputs
check_BioGeoBEARS_run(BAYES)

runslow = TRUE
resfnbayes = "AniliosBAYES.Rdata"
if (runslow)
{
  resbayes = bears_optim_run(BAYES)
  resbayes    
  
  save(resbayes, file=resfnbayes)
  resBAYAREALIKE = resbayes
} else {
  # Loads to "res"
  load(resfnbayes)
  resBAYAREALIKE = resbayes
}

#######################
#####BAYAREALIKE+J#####
#######################

BAYESj = define_BioGeoBEARS_run()
BAYESj$trfn = trfn
BAYESj$geogfn = geogfn
BAYESj$max_range_size = max_range_size
BAYESj$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYESj$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

# BAYESj$areas_allowed_fn = area.all
# BAYESj$timesfn = tiempos

# Speed options and multicore processing if desired
BAYESj$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BAYESj$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BAYESj$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BAYESj$num_cores_to_use = 4
BAYESj$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
BAYES$states_list = states_list # Need to define manually

BAYESj = readfiles_BioGeoBEARS_run(BAYESj)
BAYESj = section_the_tree(inputs=BAYESj, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BAYESj$return_condlikes_table = TRUE
BAYESj$calc_TTL_loglike_from_condlikes_table = TRUE
BAYESj$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resBAYAREALIKE$outputs@params_table["d","est"]
estart = resBAYAREALIKE$outputs@params_table["e","est"]
jstart = 0.0001

# Set up BAYAREALIKE model
# No subset sympatry
BAYESj$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BAYESj$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BAYESj$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BAYESj$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BAYESj$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BAYESj$BioGeoBEARS_model_object@params_table["e","init"] = estart
BAYESj$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No vicariance
BAYESj$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BAYESj$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BAYESj$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BAYESj$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BAYESj$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BAYESj$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BAYESj$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BAYESj$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BAYESj$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BAYESj$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BAYESj$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BAYESj$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BAYESj$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

BAYESj$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BAYESj$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BAYESj$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BAYESj$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BAYESj$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BAYESj$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Check the inputs
check_BioGeoBEARS_run(BAYESj)

runslow = TRUE
resfnbayesj = "AniliosBAYES_J.Rdata"
if (runslow)
{
  resbayesj = bears_optim_run(BAYESj)
  resbayesj    
  
  save(resbayesj, file=resfnbayesj)
  resBAYAREALIKEj = resbayesj
} else {
  # Loads to "res"
  load(resfnbayesj)
  resBAYAREALIKEj = resbayesj
}

pdffn = "Anilios_BAYAREA.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - BayAREALIKE
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE on Anilios"

# Setup
results_object = resBAYAREALIKE
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", 
                                label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                                cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie",
                         label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, 
                         cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

#######################################################
# Plot ancestral states - BAYAREALIKE+J
#######################################################
analysis_titletxt ="BioGeoBEARS BAYAREALIKE+J on Anilios"

# Setup
results_object = resBAYAREALIKEj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tree.c, tipranges=tipranges)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)







##############
#####TEST#####
##############

restable = NULL
teststable = NULL

## Result table ##

fin.dec = extract_params_from_BioGeoBEARS_results_object(results_object=resDEC, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
fin.decj = extract_params_from_BioGeoBEARS_results_object(results_object=resDECj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# fin.decx = extract_params_from_BioGeoBEARS_results_object(results_object=resDECx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# fin.decjx = extract_params_from_BioGeoBEARS_results_object(results_object=resDECjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.diva = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
fin.divaj = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# fin.divax = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# fin.divajx = extract_params_from_BioGeoBEARS_results_object(results_object=resDIVALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
fin.bayes = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKE, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
fin.bayesj = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEj, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
# fin.bayesx = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)
# fin.bayesjx = extract_params_from_BioGeoBEARS_results_object(results_object=resBAYAREALIKEjx, returnwhat="table", addl_params=c("j","x"), paramsstr_digits=4)

restable <- rbind(fin.dec,fin.decj, fin.diva,fin.divaj, fin.bayes,fin.bayesj)
row.names(restable) = c("fin.dec","fin.decj", "fin.diva","fin.divaj","fin.bayes","fin.bayesj")
restable = put_jcol_after_ecol(restable)

arb <- ape::read.tree('data/tree/anilios_newick_st.tre')
restable$AICc <- calc_AICc_vals(LnL_vals=restable$LnL,nparam_vals=restable$numparams,samplesize=length(arb$tip.label))

save(restable, file="restable.Rdata")
# load("restable.Rdata")

## Test Table ##
aicc.vals <- restable$AICc
names(aicc.vals) <- rownames(restable)

(aicc.pref <- names(sort(aicc.vals))[1])

lrt.fin <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.decj")],
                   LnL_2= restable$LnL[which(rownames(restable)=="fin.dec")],
                   numparams1= restable$numparams[which(rownames(restable)=="fin.decj")],
                   numparams2= restable$numparams[which(rownames(restable)=="fin.dec")],
                   returnwhat = "all")

lrt.fin2 <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.divaj")],
                    LnL_2= restable$LnL[which(rownames(restable)=="fin.diva")],
                    numparams1= restable$numparams[which(rownames(restable)=="fin.divaj")],
                    numparams2= restable$numparams[which(rownames(restable)=="fin.diva")],
                    returnwhat = "all")

lrt.fin3 <- lrttest(LnL_1= restable$LnL[which(rownames(restable)=="fin.bayesj")],
                    LnL_2= restable$LnL[which(rownames(restable)=="fin.bayes")],
                    numparams1= restable$numparams[which(rownames(restable)=="fin.bayesj")],
                    numparams2= restable$numparams[which(rownames(restable)=="fin.bayes")],
                    returnwhat = "all")

result_table <- rbind(lrt.fin, lrt.fin2, lrt.fin3)

result_table$alt <- c("DECj", "DIVAj", "BAYAREAj")
result_table$null <- c("DEC", "DIVA", "BAYAREA")

result_table$LnLalt <- as.numeric(result_table$LnLalt) %>% round(., digits = 2)
result_table$LnLnull <- as.numeric(result_table$LnLnull) %>% round(., digits = 2)
result_table$Dstatistic <- as.numeric(result_table$Dstatistic) %>% round(., digit = 3)
result_table$pval <- as.numeric(result_table$pval) %>% round(., digit = 4)

result_table

save(result_table, file = "result_table_prelim.Rdata")

### More details 
lrt_tester <- function(res_alt, res_null, ip_1, ip_2){
  lnl_alt = get_LnL_from_BioGeoBEARS_results_object(res_alt)
  lnl_null = get_LnL_from_BioGeoBEARS_results_object(res_null)
  f_numparams1 = ip_1;  f_numparams2 = ip_2
  alt_name <- gsub(x = deparse(substitute(res_alt)), pattern = "res", replacement = "")
  null_name <- gsub(x = deparse(substitute(res_null)), pattern = "res", replacement = "")
  
  stats = AICstats_2models(lnl_alt, lnl_null, f_numparams1, f_numparams2)
  
  # Null model
  f_res2 = extract_params_from_BioGeoBEARS_results_object(results_object=res_null, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  # Alternative model
  f_res1 = extract_params_from_BioGeoBEARS_results_object(results_object=res_alt, returnwhat="table", addl_params=c("j"), paramsstr_digits=4)
  
  stats$null <- null_name; stats$alt <- alt_name
  
  stats[, 3:9] <- as.numeric(stats[, 3:9])
  
  return(stats)
}


lrt_DEC <- lrt_tester(res_alt =  resDECj, res_null = resDEC, 
                      ip_1 = restable$numparams[which(rownames(restable)=="fin.decj")],
                      ip_2 = restable$numparams[which(rownames(restable)=="fin.dec")]
)

lrt_DIVA <- lrt_tester(res_alt =  resDIVALIKEj, res_null = resDIVALIKE, 
                       ip_1 = restable$numparams[which(rownames(restable)=="fin.divaj")],
                       ip_2 = restable$numparams[which(rownames(restable)=="fin.diva")]
)

lrt_BAYAREA <- lrt_tester(res_alt =  resBAYAREALIKEj, res_null = resBAYAREALIKE, 
                          ip_1 = restable$numparams[which(rownames(restable)=="fin.bayesj")],
                          ip_2 = restable$numparams[which(rownames(restable)=="fin.bayes")]
)

rbind(lrt_DEC, lrt_DIVA, lrt_BAYAREA)
