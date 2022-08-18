# 00_diversification_tree_data_prep.R
# Putter Tiatragul
# Started on 27 May 2022

# The purpose of this script is to prepare the data required for diversification analyses.
# --- Tree data -- SqCL exon-capture tree (two versions depends on prior)
# --- Environment data -- global paleo climate, aridity, and precipitation (Alex Skeels)

# libraries ---------------------------------------------------------------

library(tidyverse); library(RPANDA)
library(phytools); library(ape)
library(gridExtra); library(ggplot2)

# phylogeny ---------------------------------------------------------------
# Time calibrated tree from SqCL probe and mcmctree analysis with fossils

# Full tree with fossils
fos_tree_full <- ape::read.nexus("data/tree/20220420mcmctree.tre")
phy_full <- ape::read.nexus("data/tree/anilios_v3_st.tre")
phy_full_b <- ape::read.nexus("data/tree/anilios_v3_b.tre")

# Rename cfwaitii
phy_full$tip.label[which(phy_full$tip.label == "Ramphotyphlops_cfwaitii_R51244")] <- "Anilios_centralis_SAMA_R51244"
phy_full_b$tip.label[which(phy_full_b$tip.label == "Ramphotyphlops_cfwaitii_R51244")] <- "Anilios_centralis_SAMA_R51244"

# Rename polygrammicus to Sundatyphlops
phy_full$tip.label[which(phy_full$tip.label == "Anilios_polygrammicus_R98715")] <- "Sundatyphlops_polygrammicus_R98715"


# Subset species tree to Acutotyphlops & Anilios ------------------------------------------

# Additional species to drop
# A. splendidus, we decided to synonymise this species with A. pinguis
# A. centralis R51244 is redundant
# only need one unguirostris at the moment. Will only keep 180002
to_drop <- c("Anilios_splendidus_R119900", "Anilios_centralis_R51244",
             "Anilios_unguirostris_R115861", "Anilios_unguirostris_R115861", "Anilios_unguirostris_R21669", "Anilios_unguirostris_R54430")

plot(phy_full);nodelabels()

phy_small <- phy_full %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

phy_small_b <- phy_full_b %>% 
  ape::extract.clade(phy = ., 159) %>% 
  ape::drop.tip(phy = ., tip = to_drop)

# Rename tips (for grypus, ligatus we have three potential species)
phy_small$tip.label[which(phy_small$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297', 'Anilios_ligatus_R31019', 'Anilios_grypus_R55272'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_ligatusE_R31019', 'Anilios_grypusET_R55272')


# Just keep the species name (without rego)
fos_tree$tip.label <- stringr::str_replace(string = fos_tree$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")


# Drop additional taxa' ---------------------------------------------------



to_drop <- c("Anilios_splendidus_R119900", "Anilios_centralis_R51244",
             "Anilios_unguirostris_R115861", "Anilios_unguirostris_R115861", "Anilios_unguirostris_R21669", "Anilios_unguirostris_R54430")

phy_small <- 

# tips we want for blindsnakes
blindsnakes_samples <- c('Anilios_affinis_J93635','Anilios_ammodytes_R158097','Anilios_aspina_J91822','Anilios_australis_R115859',
                         'Anilios_bicolor_R55342','Anilios_bituberculatus_R63990','Anilios_broomi_J87805','Anilios_centralis_R14631',
                         'Anilios_diversus_R112027','Anilios_endoterus_R102725','Anilios_ganei_R165000', 
                         'Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272',
                         'Anilios_guentheri_R108431','Anilios_hamatus_R136277','Anilios_howi_R146381','Anilios_kimberleyensis_R164213', 
                         'Anilios_leptosoma_R119241','Anilios_leucoproctus_J87547', 
                         'Anilios_ligatus_R31019','Anilios_longissimus_R120049','Anilios_margaretae_R163269',
                         'Anilios_nigrescens_R31022','Anilios_obtusifrons_R146400','Anilios_pilbarensis_R108813',
                         'Anilios_pinguis_R146995', 'Anilios_proximus_R132486', 
                         'Anilios_silvia_J46128','Anilios_splendidus_R119900','Anilios_systenos_R114894', 
                         'Anilios_torresianus_J47640','Anilios_tovelli_R21063','Anilios_troglodytes_R146048', 
                         'Anilios_unguirostris_R21669','Anilios_waitii_R113302','Anilios_wiedii_J59020', 
                         'Anilios_yirrikalae_J85695','Ramphotyphlops_cfwaitii_R51244', 
                         'Anilios_ligatus_R19109'
)

# Subset the tree
fos_tree <- ape::keep.tip(phy = fos_tree_full, tip = blindsnakes_samples)

plotTree(fos_tree)
# Rename tips (for grypus we have three potential species)
fos_tree$tip.label[which(fos_tree$tip.label %in% c('Anilios_grypus_R108596','Anilios_grypus_R157297','Anilios_grypus_R55272', 'Anilios_ligatus_R31019'))] <- c('Anilios_grypusW_R108596', 'Anilios_grypusNW_R157297', 'Anilios_grypusET_R55272', 'Anilios_ligatusE_R31019')

# Just keep the species name (without rego)
fos_tree$tip.label <- stringr::str_replace(string = fos_tree$tip.label, pattern = regex("_([A-z][0-9]+)"), replacement = "")

# Climate data ------------------------------------------------------------

# Global temp from RPANDA (isotopes)
data("InfTemp")

# Australia temp and aridity from Alex Skeels
# --- aridity index (coded as p = precipitation because mislabeled from original data set
# --- higher values means more arid)
climate_data_sco <- read.csv("data/Australia_climate_data_45Mya_Scotese.csv")
climate_data_str <- read.csv("data/Australia_climate_data_45Mya_Straume.csv")
climate_data_val <- read.csv("data/Australia_climate_data_45Mya_Valdes_160K.csv") %>% dplyr::select(-prec_mean, -prec_median, -prec_sd)

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
                      # min_sco = climate_data_sco[, c("time_mya", "t_min")],
                      arid_sco = climate_data_sco[, c("time_mya", "p_mean")],
                      # arid_prop_sco = climate_data_sco[, c("time_mya", "prop_arid")],
                      mean_str = climate_data_str[, c("time_mya", "t_mean")],
                      # min_str = climate_data_str[, c("time_mya", "t_min")],
                      arid_str = climate_data_str[, c("time_mya", "p_mean")],
                      # arid_prop_str = climate_data_str[, c("time_mya", "prop_arid")],
                      mean_val = climate_data_val[, c("time_mya", "t_mean")],
                      # min_val = climate_data_val[, c("time_mya", "t_min")],
                      arid_val = climate_data_val[, c("time_mya", "p_mean")],
                      # arid_prop_val = climate_data_val[, c("time_mya", "prop_arid")],
                      global_temp = subset_InfTemp
)

