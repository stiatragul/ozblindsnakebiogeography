# 02_diversification_plots.R
# Putter
# Plotting diversification models


# libraries ---------------------------------------------------------------

library(RPANDA); library(dplyr)
library(phytools); library(gridExtra)

# Load results from 01_diversification_analyses ---------------------------
# These models take a lot of time to run so I ran it on WSL2 and saved the image file. 

load("output/20220818_model_output_env.RData")

source("code/02_table_fit_env_results.R")

# Best fitting model
result_fit$arid_sco[[1]]$BEnvVarDCST_EXPO
result_fit$min_sco[[1]]$BCSTDEnvVar_EXPO


RPANDA::plot_fit_env(result_fit$arid_sco[[1]]$BEnvVar_EXPO, env_data = env_data_list$arid_sco, tot_time = tot_time)
RPANDA::plot_fit_env(result_fit$min_sco[[1]]$BCSTDEnvVar_EXPO, env_data = env_data_list$min_sco, tot_time = tot_time)
