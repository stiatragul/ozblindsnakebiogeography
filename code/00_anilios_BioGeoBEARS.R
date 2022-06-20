# 00_anilios_BioGeoBEARS.R


# Libraries ---------------------------------------------------------------

library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
# library(parallel)
library(snow)
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)

