# 02_diversification_plots.R
# Putter
# Plotting diversification models


# libraries ---------------------------------------------------------------

library(RPANDA); library(dplyr)
library(phytools); library(gridExtra)

# Load results from 01_diversification_analyses ---------------------------
# These models take a lot of time to run so I ran it on WSL2 and saved the image file. 

load("output/20220818_model_output_env.RData")


# Sort result table -------------------------------------------------------

result_table <- rbind(tab_func(tree_BCSTDenvVar_EXPO_res, .mu = TRUE, gen_name = "BCSTDenvVar_expo"),
      tab_func(tree_BEnvVarDCST_EXPO_res, .mu = TRUE, gen_name = "BEnvVarDCST_expo"),
      tab_func(tree_BEnvVarDEnvVar_EXPO_res, .mu = TRUE, gen_name = "BEnvVarDEnvVar_expo"),
      tab_func(tree_BEnvVar_EXPO_res, .mu = FALSE, gen_name = "BEnvVar_expo")) %>% 
  arrange(aicc)

head(result_table)

result_table %>% 
  dplyr::mutate_if(is.numeric, round, digits =2)


tree_BEnvVar_EXPO_res$arid_sco

RPANDA::plot_fit_env(tree_BEnvVar_EXPO_res$arid_str, env_data = env_data_list$arid_str, tot_time = tot_time)
