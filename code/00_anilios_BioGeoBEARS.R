# 00_anilios_BioGeoBEARS.R


# Libraries ---------------------------------------------------------------

library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
# library(parallel)
library(snow)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

tree.c <- ape::read.tree('data/tree/anilios_newick.tre')

trfn <- "data/tree/anilios_newick.tre"
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
states_list = rcpp_areas_list_to_states_list(areas=coln, maxareas=maxAreas, include_null_range=TRUE)

# Set up list of TRUEs
states_to_keep = rep(TRUE, times=length(states_list))

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
states_list = states_list[states_to_keep]
# https://groups.google.com/d/msg/biogeobears/zkHP_YpA9kY/8c3vR-gIqP0J
statenames = areas_list_to_states_list_new(areas=coln, maxareas=maxAreas, include_null_range=TRUE, split_ABC=FALSE)
statenames = statenames[states_to_keep]
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
DEC$states_list = states_list # Need to define manually

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
# DECj$areas_allowed_fn = area.all
# DECj$timesfn = tiempos
DECj$states_list = states_list # Need to define manually

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

