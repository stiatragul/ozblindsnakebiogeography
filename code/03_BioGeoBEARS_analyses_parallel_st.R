# 03_anilios_BioGeoBEARS_parallel_st.R
# August 2022
# THIS one runs with tree using skewT and skew normal priors only

# Script to fit different biogeography models with BioGeoBEARS and save results.
# base code from BioGeoBEARS website, Octavio Jimenez Robles and Carlos Pavon Vazquez

# Notes:
# +J is for jump dispersal; + X is to account for distance between biome
# Models: DEC, BAYAREALIKE, DIVALIKE
# Fit models, models+j, models+x, models+j+x

# Libraries ---------------------------------------------------------------

library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(parallel); library(snow)
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)

# Check tree
tree.c <- ape::read.tree('data/intermediate_data/bears/blindsnake_st.tre')
plot(tree.c)

# Load inputs -------------------------------------------------------------

# Tree file
trfn <- "data/intermediate_data/bears/blindsnake_st.tre"
# Distance matrix file 
dist.mat <- "data/bears_txt/biome_distance.txt"
# Geo file
geogfn <- "data/bears_txt/geofile.txt"
# adj <- "data/bears_txt/area_adj.txt"



# Specifications ----------------------------------------------------------

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 5

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

coln <- c("A", "B", "C", "D", "E", "F", "G", "H", "I")

# Manually setting up area-adjacency using Octavio's code
#                     A   B   C   D   E   F   G   H   I
adjacency <- matrix(c(1,	1,	1,	1,	1,	1,	1,	1,	1,  # A
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # B
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # C
                      1,	1,	1,	1,	0,	1,	0,	0,	0,  # D
                      1,	1,	0,	0,	1,	1,	0,	0,	1,  # E
                      1,	1,	1,	1,	1,	1,	0,	0,	0,  # F
                      1,	0,	0,	0,	0,	0,	1,	0,	1,  # G
                      0,	0,	0,	0,	0,	0,	0,	1,	1,  # H
                      0,	0,	0,	0,	0,	0,	1,	1,	1), # I
                    ncol=9, nrow=9, byrow=T)

maxAreas = 5
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

# Modify the states_list to reduce analyses time
states_list <- states_list[states_to_keep]
# https://groups.google.com/d/msg/biogeobears/zkHP_YpA9kY/8c3vR-gIqP0J
statenames <- areas_list_to_states_list_new(areas=coln, maxareas=maxAreas, include_null_range=TRUE, split_ABC=FALSE)
statenames <- statenames[states_to_keep]
unlist(statenames)

# Number of cores if run parallel
number_cores <- detectCores()/2


# Fitting models ----------------------------------------------------------


############
#####DEC####
############

####

# NOTE: 'DEC' is identical with Lagrange DEC model.

DEC <- define_BioGeoBEARS_run()
DEC$trfn = trfn
DEC$geogfn = geogfn
DEC$max_range_size = max_range_size
DEC$min_branchlength = 0.000001; DEC$include_null_range = TRUE
DEC$on_NaN_error = -1e50; DEC$speedup = TRUE; DEC$use_optimx = "GenSA"
DEC$num_cores_to_use = number_cores

DEC = readfiles_BioGeoBEARS_run(DEC)

DEC$return_condlikes_table = TRUE
DEC$calc_TTL_loglike_from_condlikes_table = TRUE
DEC$calc_ancprobs = TRUE    # get ancestral states from optim run
DEC$states_list = states_list # Need to define manually (using code from Octavio)

check_BioGeoBEARS_run(DEC)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = "data/intermediate_data/bears/AniliosDEC_st.Rdata"
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
# D = dispersal; E = Extinction; C = cladogenesis; + J is for 'jump'


#  Regular set up here
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
DECj$states_list = states_list # Need to define manually. This replaces areas_allowed_fn

### For time periods ###
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

resfnj = "data/intermediate_data/bears/AniliosDEC_J_st.Rdata"
runslow = TRUE
if (runslow)
{
  resj = bears_optim_run(DECj)
  resj    
  save(resj, file=resfnj)
  resDECj = resj
} else {
  # Loads to "res"
  load(resfnj)
  resDECj = resj
}




###############
#####DEC+X#####
###############
# Distance matrix

DECx <- define_BioGeoBEARS_run()
DECx$trfn = trfn
DECx$geogfn = geogfn
DECx$max_range_size = max_range_size
DECx$min_branchlength = 0.000001
DECx$include_null_range = TRUE
DECx$on_NaN_error = -1e50
#DECx$speedup = FALSE  
#DECx$use_optimx = FALSE
DECx$num_cores_to_use = number_cores
DECx$states_list = states_list
DECx$distsfn = dist.mat
# DECx$timesfn = tiempos

DECx = readfiles_BioGeoBEARS_run(DECx)
# DECx = section_the_tree(inputs=DECx, make_master_table=TRUE, plot_pieces=FALSE)

DECx$return_condlikes_table = TRUE
DECx$calc_TTL_loglike_from_condlikes_table = TRUE
DECx$calc_ancprobs = TRUE    # get ancestral states from optim run

# Add distance-dependent dispersal (+x)
DECx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
DECx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
DECx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
DECx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
DECx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(DECx)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfnx = "data/intermediate_data/bears/AniliosDEC_X_st.Rdata"
if (runslow)
{
  resx = bears_optim_run(DECx)
  resx    
  
  save(resx, file=resfnx)
  resDECx = resx
} else {
  # Loads to "res"
  load(resfnx)
  resDECx = resx
}




#################
#####DEC+J+X#####
#################
# Distance matrix and jump



DECjx = define_BioGeoBEARS_run()
DECjx$trfn = trfn
DECjx$geogfn = geogfn
DECjx$max_range_size = max_range_size
DECjx$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DECjx$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

DECjx$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
#DECjx$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
#DECjx$use_optimx = FALSE    # if FALSE, use optim() instead of optimx()
DECjx$num_cores_to_use = number_cores
DECjx$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale
DECjx$states_list = states_list
DECjx$distsfn = dist.mat
# DECjx$timesfn = tiempos

DECjx = readfiles_BioGeoBEARS_run(DECjx)
# DECjx = section_the_tree(inputs=DECjx, make_master_table=TRUE, plot_pieces=FALSE)

DECjx$return_condlikes_table = TRUE
DECjx$calc_TTL_loglike_from_condlikes_table = TRUE
DECjx$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDECx$outputs@params_table["d","est"]
estart = resDECx$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
DECjx$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DECjx$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DECjx$BioGeoBEARS_model_object@params_table["e","init"] = estart
DECjx$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
DECjx$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DECjx$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DECjx$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Add distance-dependent dispersal (+x)
DECjx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
DECjx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
DECjx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
DECjx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
DECjx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(DECjx)

resfnjx = "data/intermediate_data/bears/AniliosDEC_JX_st.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  resjx = bears_optim_run(DECjx)
  resjx    
  
  save(resjx, file=resfnjx)
  
  resDECjx = resjx
} else {
  # Loads to "res"
  load(resfnjx)
  resDECjx = resjx
}



#######################################
### DIVALIKE AND DIVALIKE+J ANALYSIS ##
#######################################

# NOTE: The BioGeoBEARS "DIVALIKE" model is not identical with 
# Ronquist (1997)'s parsimony DIVA. It is a likelihood
# interpretation of DIVA, constructed by modelling DIVA's
# processes the way DEC does, but only allowing the 
# processes DIVA allows (widespread vicariance: yes; subset
# sympatry: no; see Ronquist & Sanmartin 2011, Figure 4).
# DIVALIKE is a likelihood interpretation of parsimony
# DIVA, and it is "like DIVA" -- similar to, but not
# identical to, parsimony DIVA.



##################
#####DIVALIKE#####
##################


DIVA = define_BioGeoBEARS_run()
DIVA$trfn = trfn
DIVA$geogfn = geogfn
DIVA$max_range_size = max_range_size
DIVA$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVA$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
DIVA$states_list = states_list
# DIVA$areas_allowed_fn = area.all
# DIVA$timesfn = tiempos

# Speed options and multicore processing if desired
DIVA$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
DIVA$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
DIVA$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
DIVA$num_cores_to_use = number_cores
DIVA$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

DIVA = readfiles_BioGeoBEARS_run(DIVA)
# DIVA = section_the_tree(inputs=DIVA, make_master_table=TRUE, plot_pieces=FALSE)

DIVA$return_condlikes_table = TRUE
DIVA$calc_TTL_loglike_from_condlikes_table = TRUE
DIVA$calc_ancprobs = TRUE    # get ancestral states from optim run

DIVA$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVA$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVA$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVA$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVA$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVA$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVA$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

DIVA$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVA$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVA$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

check_BioGeoBEARS_run(DIVA)

runslow = TRUE
resfndiva = "data/intermediate_data/bears/AniliosDIVA_st.Rdata"
if (runslow)
{
  resdiva = bears_optim_run(DIVA)
  resdiva    
  
  save(resdiva, file=resfndiva)
  resDIVALIKE = resdiva
} else {
  # Loads to "res"
  load(resfndiva)
  resDIVALIKE = resdiva
}


####################
#####DIVALIKE+J#####
####################

DIVAj = define_BioGeoBEARS_run()
DIVAj$trfn = trfn
DIVAj$geogfn = geogfn
DIVAj$max_range_size = max_range_size
DIVAj$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVAj$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
DIVAj$states_list = states_list
# DIVAj$areas_allowed_fn = area.all
# DIVAj$timesfn = tiempos

# Speed options and multicore processing if desired
DIVAj$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
DIVAj$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
DIVAj$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
DIVAj$num_cores_to_use = 4
DIVAj$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

DIVAj = readfiles_BioGeoBEARS_run(DIVAj)
# DIVAj = section_the_tree(inputs=DIVAj, make_master_table=TRUE, plot_pieces=FALSE)

DIVAj$return_condlikes_table = TRUE
DIVAj$calc_TTL_loglike_from_condlikes_table = TRUE
DIVAj$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resDIVALIKE$outputs@params_table["d","est"]
estart = resDIVALIKE$outputs@params_table["e","est"]
jstart = 0.0001

DIVAj$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DIVAj$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DIVAj$BioGeoBEARS_model_object@params_table["e","init"] = estart
DIVAj$BioGeoBEARS_model_object@params_table["e","est"] = estart

DIVAj$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVAj$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVAj$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVAj$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVAj$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVAj$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVAj$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

DIVAj$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVAj$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVAj$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

DIVAj$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DIVAj$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DIVAj$BioGeoBEARS_model_object@params_table["j","est"] = jstart

DIVAj$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
DIVAj$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

check_BioGeoBEARS_run(DIVAj)

runslow = TRUE
resfndivaj = "AniliosDIVA_J_st.Rdata"
if (runslow)
{
  resdivaj = bears_optim_run(DIVAj)
  resdivaj    
  
  save(resdivaj, file=resfndivaj)
  resDIVALIKEj = resdivaj
} else {
  # Loads to "res"
  load(resfndivaj)
  resDIVALIKEj = resdivaj
}



####################
#####DIVALIKE+X#####
####################
# Includes distance matrix

DIVAx = define_BioGeoBEARS_run()
DIVAx$trfn = trfn
DIVAx$geogfn = geogfn
DIVAx$max_range_size = max_range_size
DIVAx$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVAx$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
DIVAx$states_list = states_list
DIVAx$distsfn = dist.mat
# DIVAx$timesfn = tiempos

# Speed options and multicore processing if desired
DIVAx$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
#DIVAx$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
#DIVAx$use_optimx = FALSE    # if FALSE, use optim() instead of optimx()
DIVAx$num_cores_to_use = 4
DIVAx$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

DIVAx = readfiles_BioGeoBEARS_run(DIVAx)
# DIVAx = section_the_tree(inputs=DIVAx, make_master_table=TRUE, plot_pieces=FALSE)

DIVAx$return_condlikes_table = TRUE
DIVAx$calc_TTL_loglike_from_condlikes_table = TRUE
DIVAx$calc_ancprobs = TRUE    # get ancestral states from optim run

DIVAx$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVAx$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVAx$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVAx$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVAx$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVAx$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVAx$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

DIVAx$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVAx$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVAx$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

# Add distance-dependent dispersal (+x)
DIVAx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
DIVAx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
DIVAx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
DIVAx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
DIVAx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(DIVAx)

runslow = TRUE
resfndivax = "AniliosDIVA_X_st.Rdata"
if (runslow)
{
  resdivax = bears_optim_run(DIVAx)
  resdivax    
  
  save(resdivax, file=resfndivax)
  resDIVALIKEx = resdivax
} else {
  # Loads to "res"
  load(resfndivax)
  resDIVALIKEx = resdivax
}



######################
#####DIVALIKE+J+X#####
######################
# Distance matrix (x) and jump dispersal


DIVAjx = define_BioGeoBEARS_run()
DIVAjx$trfn = trfn
DIVAjx$geogfn = geogfn
DIVAjx$max_range_size = max_range_size
DIVAjx$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
DIVAjx$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
DIVAjx$states_list = states_list
DIVAjx$distsfn = dist.mat
# DIVAjx$timesfn = tiempos

# Speed options and multicore processing if desired
DIVAjx$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
#DIVAjx$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
#DIVAjx$use_optimx = FALSE    # if FALSE, use optim() instead of optimx()
DIVAjx$num_cores_to_use = number_cores
DIVAjx$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

DIVAjx = readfiles_BioGeoBEARS_run(DIVAjx)
# DIVAjx = section_the_tree(inputs=DIVAjx, make_master_table=TRUE, plot_pieces=FALSE)

DIVAjx$return_condlikes_table = TRUE
DIVAjx$calc_TTL_loglike_from_condlikes_table = TRUE
DIVAjx$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resDIVALIKEx$outputs@params_table["d","est"]
estart = resDIVALIKEx$outputs@params_table["e","est"]
jstart = 0.0001

DIVAjx$BioGeoBEARS_model_object@params_table["d","init"] = dstart
DIVAjx$BioGeoBEARS_model_object@params_table["d","est"] = dstart
DIVAjx$BioGeoBEARS_model_object@params_table["e","init"] = estart
DIVAjx$BioGeoBEARS_model_object@params_table["e","est"] = estart

DIVAjx$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
DIVAjx$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
DIVAjx$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

DIVAjx$BioGeoBEARS_model_object@params_table["ysv","type"] = "2-j"
DIVAjx$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/2"
DIVAjx$BioGeoBEARS_model_object@params_table["y","type"] = "ysv*1/2"
DIVAjx$BioGeoBEARS_model_object@params_table["v","type"] = "ysv*1/2"

DIVAjx$BioGeoBEARS_model_object@params_table["mx01v","type"] = "fixed"
DIVAjx$BioGeoBEARS_model_object@params_table["mx01v","init"] = 0.5
DIVAjx$BioGeoBEARS_model_object@params_table["mx01v","est"] = 0.5

DIVAjx$BioGeoBEARS_model_object@params_table["j","type"] = "free"
DIVAjx$BioGeoBEARS_model_object@params_table["j","init"] = jstart
DIVAjx$BioGeoBEARS_model_object@params_table["j","est"] = jstart

DIVAjx$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
DIVAjx$BioGeoBEARS_model_object@params_table["j","max"] = 1.99999

# Add distance-dependent dispersal (+x)
DIVAjx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
DIVAjx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
DIVAjx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
DIVAjx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
DIVAjx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

check_BioGeoBEARS_run(DIVAjx)

runslow = TRUE
resfndivajx = "data/intermediate_data/bears/AniliosDIVA_JX_st.Rdata"
if (runslow)
{
  resdivajx = bears_optim_run(DIVAjx)
  resdivajx    
  
  save(resdivajx, file=resfndivajx)
  resDIVALIKEjx = resdivajx
} else {
  # Loads to "res"
  load(resfndivajx)
  resDIVALIKEjx = resdivajx
}



#####################
#####BAYAREALIKE#####
#####################
# Assumption: nothig in particular happens at cladogenesis

BAYES = define_BioGeoBEARS_run()
BAYES$trfn = trfn
BAYES$geogfn = geogfn
BAYES$max_range_size = max_range_size
BAYES$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYES$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BAYES$states_list = states_list
# BAYES$timesfn = tiempos

# Speed options and multicore processing if desired
BAYES$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BAYES$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BAYES$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BAYES$num_cores_to_use = number_cores
BAYES$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

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
resfnbayes = "data/intermediate_data/bears/AniliosBAYES_st.Rdata"
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
# + jump parameter



BAYESj = define_BioGeoBEARS_run()
BAYESj$trfn = trfn
BAYESj$geogfn = geogfn
BAYESj$max_range_size = max_range_size
BAYESj$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYESj$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BAYESj$states_list = states_list
# BAYESj$timesfn = tiempos

# Speed options and multicore processing if desired
BAYESj$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BAYESj$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BAYESj$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BAYESj$num_cores_to_use = number_cores
BAYESj$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

BAYESj = readfiles_BioGeoBEARS_run(BAYESj)
# BAYESj = section_the_tree(inputs=BAYESj, make_master_table=TRUE, plot_pieces=FALSE)

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
resfnbayesj = "data/intermediate_data/bears/AniliosBAYES_J_st.Rdata"
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




#######################
#####BAYAREALIKE+X#####
#######################
# +x Distance matrix 



BAYESx = define_BioGeoBEARS_run()
BAYESx$trfn = trfn
BAYESx$geogfn = geogfn
BAYESx$max_range_size = max_range_size
BAYESx$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYESx$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.

BAYESx$states_list = states_list
BAYESx$distsfn = dist.mat
# BAYESx$timesfn = tiempos

# Speed options and multicore processing if desired
BAYESx$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
#BAYESx$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
#BAYESx$use_optimx = FALSE    # if FALSE, use optim() instead of optimx()
BAYESx$num_cores_to_use = number_cores
BAYESx$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BAYESx = readfiles_BioGeoBEARS_run(BAYESx)
# BAYESx = section_the_tree(inputs=BAYESx, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BAYESx$return_condlikes_table = TRUE
BAYESx$calc_TTL_loglike_from_condlikes_table = TRUE
BAYESx$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up BAYAREALIKE model
# No subset sympatry
BAYESx$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BAYESx$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BAYESx$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

# No vicariance
BAYESx$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BAYESx$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BAYESx$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# Adjust linkage between parameters
BAYESx$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BAYESx$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BAYESx$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BAYESx$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BAYESx$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BAYESx$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

# Add distance-dependent dispersal (+x)
BAYESx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BAYESx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BAYESx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BAYESx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BAYESx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

# Check the inputs
check_BioGeoBEARS_run(BAYESx)

runslow = TRUE
resfnbayesx = "data/intermediate_data/bears/AniliosBAYES_X_st.Rdata"
if (runslow)
{
  resbayesx = bears_optim_run(BAYESx)
  resbayesx   
  
  save(resbayesx, file=resfnbayesx)
  resBAYAREALIKEx = resbayesx
} else {
  # Loads to "res"
  load(resfnbayesx)
  resBAYAREALIKEx = resbayesx
}




#########################
#####BAYAREALIKE+J+X#####
#########################
# + jump and distance



BAYESjx = define_BioGeoBEARS_run()
BAYESjx$trfn = trfn
BAYESjx$geogfn = geogfn
BAYESjx$max_range_size = max_range_size
BAYESjx$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BAYESjx$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
BAYESjx$states_list = states_list
BAYESjx$distsfn = dist.mat
# BAYESjx$timesfn = tiempos

# Speed options and multicore processing if desired
BAYESjx$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
#BAYESjx$speedup = FALSE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
#BAYESjx$use_optimx = FALSE    # if FALSE, use optim() instead of optimx()
BAYESjx$num_cores_to_use = number_cores
BAYESjx$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

BAYESjx = readfiles_BioGeoBEARS_run(BAYESjx)
# BAYESjx = section_the_tree(inputs=BAYESjx, make_master_table=TRUE, plot_pieces=FALSE)

# Good default settings to get ancestral states
BAYESjx$return_condlikes_table = TRUE
BAYESjx$calc_TTL_loglike_from_condlikes_table = TRUE
BAYESjx$calc_ancprobs = TRUE    # get ancestral states from optim run

dstart = resBAYAREALIKEx$outputs@params_table["d","est"]
estart = resBAYAREALIKEx$outputs@params_table["e","est"]
jstart = 0.0001

# Set up BAYAREALIKE model
# No subset sympatry
BAYESjx$BioGeoBEARS_model_object@params_table["s","type"] = "fixed"
BAYESjx$BioGeoBEARS_model_object@params_table["s","init"] = 0.0
BAYESjx$BioGeoBEARS_model_object@params_table["s","est"] = 0.0

BAYESjx$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BAYESjx$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BAYESjx$BioGeoBEARS_model_object@params_table["e","init"] = estart
BAYESjx$BioGeoBEARS_model_object@params_table["e","est"] = estart

# No vicariance
BAYESjx$BioGeoBEARS_model_object@params_table["v","type"] = "fixed"
BAYESjx$BioGeoBEARS_model_object@params_table["v","init"] = 0.0
BAYESjx$BioGeoBEARS_model_object@params_table["v","est"] = 0.0

# *DO* allow jump dispersal/founder-event speciation (set the starting value close to 0)
BAYESjx$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BAYESjx$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BAYESjx$BioGeoBEARS_model_object@params_table["j","est"] = jstart

# Under BAYAREALIKE+J, the max of "j" should be 1, not 3 (as is default in DEC+J) or 2 (as in DIVALIKE+J)
BAYESjx$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Adjust linkage between parameters
BAYESjx$BioGeoBEARS_model_object@params_table["ysv","type"] = "1-j"
BAYESjx$BioGeoBEARS_model_object@params_table["ys","type"] = "ysv*1/1"
BAYESjx$BioGeoBEARS_model_object@params_table["y","type"] = "1-j"

# Only sympatric/range-copying (y) events allowed, and with 
# exact copying (both descendants always the same size as the ancestor)
BAYESjx$BioGeoBEARS_model_object@params_table["mx01y","type"] = "fixed"
BAYESjx$BioGeoBEARS_model_object@params_table["mx01y","init"] = 0.9999
BAYESjx$BioGeoBEARS_model_object@params_table["mx01y","est"] = 0.9999

BAYESjx$BioGeoBEARS_model_object@params_table["d","min"] = 0.0000001
BAYESjx$BioGeoBEARS_model_object@params_table["d","max"] = 4.9999999

BAYESjx$BioGeoBEARS_model_object@params_table["e","min"] = 0.0000001
BAYESjx$BioGeoBEARS_model_object@params_table["e","max"] = 4.9999999

BAYESjx$BioGeoBEARS_model_object@params_table["j","min"] = 0.00001
BAYESjx$BioGeoBEARS_model_object@params_table["j","max"] = 0.99999

# Add distance-dependent dispersal (+x)
BAYESjx$BioGeoBEARS_model_object@params_table["x","type"] = "free"
BAYESjx$BioGeoBEARS_model_object@params_table["x","init"] = 0.0
BAYESjx$BioGeoBEARS_model_object@params_table["x","est"] = -0.00001
BAYESjx$BioGeoBEARS_model_object@params_table["x","min"] = -5.0
BAYESjx$BioGeoBEARS_model_object@params_table["x","max"] = 0.0

# Check the inputs
check_BioGeoBEARS_run(BAYESjx)

runslow = TRUE
resfnbayesjx = "data/intermediate_data/bears/AniliosBAYES_JX_st.Rdata"
if (runslow)
{
  resbayesjx = bears_optim_run(BAYESjx)
  resbayesjx    
  
  save(resbayesjx, file=resfnbayesjx)
  resBAYAREALIKEjx = resbayesjx
} else {
  # Loads to "res"
  load(resfnbayesjx)
  resBAYAREALIKEjx = resbayesjx
}

save.image(file = "data/intermediate_data/bears/bears_model_fit_st.Rdata")
