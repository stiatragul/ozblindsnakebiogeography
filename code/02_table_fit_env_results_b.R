# 02_table_fit_env_results_b.R
# September 2022

# This script will organise the results from the parallel fit and the individual
# fit together in one table from fitting on B tree (with cauchy priors).

# From 01_fit_*.R. I successfully fitted several environment-dependent birth-death models.
# The following had to be fitted individually using initial lamb_par 0.1 and/or 0. 
# min_sco ;arid_sco; mean_str; min_str;  arid_str;  min_val; max_val; arid_val. 
# mean_val; # max_str; mean_sco; max_sco; global temp

# libraries ---------------------------------------------------------------

library(RPANDA); library(dplyr); library(phytools);
library(stringr)

# Load results ------------------------------------------------------------

# Get full path to load data the txt files with the results
file_n <- list.files(path = "data/intermediate_data/diversification_analyses/", pattern = "_b.txt", full.names = T)

# get the names of the model
file_n_m <- stringr::str_remove(list.files(path = "data/intermediate_data/diversification_analyses/",
                                           pattern = "_b.txt"), pattern = "Results_Anilios_")
mod_names <- stringr::str_remove(file_n_m, pattern = ".txt")

# List of models with their names
filelist <- lapply(file_n, read.delim)

# Environmental data in list form
load('data/intermediate_data/diversification_analyses/env_data_list.RData')

# Put environment information in the table as a column
for(i in 1:length(filelist)){
  filelist[[i]]$env <- mod_names[i]
}

# Combine tables into one big table ----------------------------------------

res_table <- do.call(rbind, filelist)

# Round the values
options("scipen"=2, "digits"=3)
res_table[3:4] <- lapply(res_table[3:4], round, 2)

# Table to present in paper
b_table <- res_table %>% arrange(desc(Lambda))
b_table[, 2:4] <- lapply(b_table[, 2:4], round, 2)
b_table$AICcWt <- round(b_table$AICcWt, 2)
b_table[5:8] <- lapply(b_table[5:8], formatC, format = "e", digits = 2)

b_table$Lambda <- as.numeric(b_table$Lambda)

# write.csv(b_table, file = 'output/prelim_RPANDA_fit.csv')



