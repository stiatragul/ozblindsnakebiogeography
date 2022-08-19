# 00_diversification_env_data_prep.R
# Putter Tiatragul
# Started on 27 May 2022

# The purpose of this script is to prepare the data required for diversification analyses.
# --- Environment data -- global paleo climate, aridity, and precipitation (Alex Skeels)
# Temperature variation based on combination of paleo digital elevation models (DEMs) and paleo climate bands. 
# Contemporary mean values for Australian Koppen bands [@scoteseAtlas2021] projected back in time by Alex Skeels (quoted).
# Reconstructions based different global DEMs by (1) Scotese @scotesePALEOMAP2018, and (2) 
# Straume @straumeGlobal2020 is reconstructed by Alex Skeels using `gen3sis` @hagenGen3sis2020.


# libraries ---------------------------------------------------------------
library(tidyverse); library(RPANDA); library(stringr)

# Climate data ------------------------------------------------------------

# Global temp from RPANDA (isotopes)
data("InfTemp")

# Australia temp and aridity from Alex Skeels
# --- aridity index (coded as p = precipitation because mislabeled from original data set
# --- higher values means more arid)
climate_data_sco <- read.csv("data/paleo_env/Australia_climate_data_45Mya_Scotese.csv")
climate_data_str <- read.csv("data/paleo_env/Australia_climate_data_45Mya_Straume.csv")
climate_data_val <- read.csv("data/paleo_env/Australia_climate_data_45Mya_Valdes_160K.csv") %>% dplyr::select(-prec_mean, -prec_median, -prec_sd)

# Multiply p_mean by 100
climate_data_sco$p_mean <- climate_data_sco$p_mean * 100
climate_data_str$p_mean <- climate_data_str$p_mean * 100
climate_data_val$p_mean <- climate_data_val$p_mean * 100

# Temporal resolution in RPANDA::InfTemp is greater than reconstructions Alex Skeels' models 
# so need to make same length by sub sampling RPANDA temp to similar timestep

# choose one random Age within a % range of Age from Scotese temperature data. 
climate_data_sco$time_mya[2]; climate_data_sco$time_mya[3]

#  get timesteps from RPANDA's dataset that already match Alex's
panda_true <- InfTemp %>% 
  dplyr::filter(Age < 46) %>% 
  dplyr::mutate(choose = ifelse(Age %in% c(climate_data_sco$time_mya), "TRUE", "FALSE")) %>% 
  dplyr::filter(choose == "TRUE") %>% dplyr::distinct(Age, .keep_all = TRUE)

# get time that did not match 
time_match <- climate_data_sco %>% 
  dplyr::mutate(need = ifelse(time_mya %in% c(panda_true$Age), "NO", "YES")) %>% 
  dplyr::filter(need == "YES") %>% 
  dplyr::select(time_mya)

# Function to sample temperature for a given Age value
TempSampler <- function(.age_value, .tempdata, .fraction){
  
  .age_value <- as.numeric(.age_value)
  p <- .fraction
  tmp <- InfTemp[which(.age_value - .age_value*p <= .tempdata$Age & .tempdata$Age <= .age_value + .age_value*p), ]
  number <- sample(tmp$Temperature, size = 1)
  return(number)
  
}

# Apply TempSampler to List of age values that still need temperature values
randomised_temp_value <- lapply(X = as.list(time_match$time_mya), 
                                FUN = TempSampler, .tempdata = InfTemp, .fraction = 0.05)

# Append the vector to the time_match dataframe
time_match$Temperature <- unlist(randomised_temp_value)
names(time_match) <- c("Age", "Temperature")

# Subset InfTemp with same time steps as Scotese and Straume
subset_InfTemp <- rbind(time_match, panda_true[1:2]) %>% sort()
names(subset_InfTemp) <- c("time_mya", "t_mean")

# Make a list of environmental data we are interested in testing in the analysis 
env_data_list <- list(mean_sco = climate_data_sco[, c("time_mya", "t_mean")],
                      min_sco = climate_data_sco[, c("time_mya", "t_min")],
                      max_sco = climate_data_sco[, c("time_mya", "t_max")],
                      arid_sco = climate_data_sco[, c("time_mya", "p_mean")],
                      # arid_prop_sco = climate_data_sco[, c("time_mya", "prop_arid")],
                      mean_str = climate_data_str[, c("time_mya", "t_mean")],
                      min_str = climate_data_str[, c("time_mya", "t_min")],
                      max_str = climate_data_str[, c("time_mya", "t_max")],
                      arid_str = climate_data_str[, c("time_mya", "p_mean")],
                      # arid_prop_str = climate_data_str[, c("time_mya", "prop_arid")],
                      mean_val = climate_data_val[, c("time_mya", "t_mean")],
                      min_val = climate_data_val[, c("time_mya", "t_min")],
                      max_val = climate_data_val[, c("time_mya", "t_max")],
                      arid_val = climate_data_val[, c("time_mya", "p_mean")],
                      # arid_prop_val = climate_data_val[, c("time_mya", "prop_arid")],
                      global_temp = subset_InfTemp
)

save(env_data_list, file = "data/intermediate_data/diversification_analyses/env_data_list.RData")
