 # 02_diversification_plots_b.R
# Putter Tiatragul
# Plotting diversification _b models

# libraries ---------------------------------------------------------------

library(RPANDA); library(dplyr)
library(phytools); library(gridExtra)
library(deeptime) # for deeptime axis
library(ggplot2); library(viridis)
library(patchwork)

# source("code/01_env_data_plots.R")
# Load results from 01_diversification_analyses ---------------------------
source("code/02_table_fit_env_results_b.R") 
source("code/utility/plot_fit_env_new.R") # load custom functions and some data

b_table$AICc <- as.numeric(b_table$AICc)

# Check best fitting model by checking lowest AICc values
b_table %>% 
  dplyr::arrange(AICc) %>% 
  head()

# Plot best fitting -------------------------------------------------------

# Best fitting models based on lowest AICc
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_arid_sco_b.Rdata'); arid_sco <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_arid_str_b.Rdata'); arid_str <- Anilios_res

# arid_sco[[1]]$BEnvVarDCST_EXPO

fit_env_arid <- arid_sco[[1]]$BEnvVarDCST_EXPO

# RPANDA Base R plot
plot_fit_env_new(fit_env_arid, 
                 env_data = env_data_list$arid_sco, tot_time = tot_time,
                 name_env = 'Aridity')


# GGPLOT
# fit_env_df_make is a custom function 
arid_sco_df <- fit_env_df_make(env_data = env_data_list$arid_sco, 
                               tot_time = tot_time,
                              fit.env = fit_env_arid)

# Custom theme
theme_1 <- function (base_size = 11, base_family = "", base_line_size = base_size/22, base_rect_size = base_size/22) {
  theme_grey(base_size = base_size,
             base_family = base_family, base_line_size = base_line_size,
             base_rect_size = base_rect_size
             ) %+replace%
    theme(panel.background = element_rect(fill = "white", colour = NA),      
          panel.grid = element_line(colour = "grey92"),
          strip.background = element_rect(fill = "grey85", colour = "grey20"),
          legend.key = element_rect(fill = "white", colour = NA),
          axis.title = element_text(size = 16),
          axis.text = element_text(size = 14),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.line = element_line(),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          complete = TRUE
    )
}
  
p2 <- arid_sco_df %>% 
  ggplot(data = ., aes(x = age)) +
  geom_line(aes(y = speciation), size = 1) +
  deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  labs(y = "Speciation (Event per Myrs)") +  
  scale_x_reverse("Age (Mya)", limits = c(26, 0), breaks = seq(25, 0, by = -5)) + 
  theme_1()

# Env data fitted
p3 <- arid_sco_df %>% 
  ggplot(data = ., aes(x = environment_data)) +
  geom_point(aes(y = speciation), size = 1.3) +
  # deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  labs(y = "Speciation (Event per Myrs)", x = "Aridity") +  
  theme_1()

# Environment -------------------------------------------------------------

# Plot environment
p4 <- env_data_list$arid_sco %>% 
  ggplot(data = ., aes(x = time_mya)) +
  geom_line(aes(y = p_mean), size = 1, colour = "orange") +
  # geom_point(aes(y = vals)) + 
  deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  # scale_colour_manual(values = viridis(4), name = "Data source", labels = c("RPANDA Global", "Scotese Australia", "Straume Australia", "Valdes Australia")) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 65, 10)) + labs(y = "Aridity index") +
  scale_x_reverse("Age (Mya)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_1() 


# Lineage through time
lineage_time <- as.data.frame(ltt.plot.coords(fos_tree))

lineage_time$time <- lineage_time$time * -1

p5 <- lineage_time %>% 
  ggplot(data = ., aes(x = time)) +
  geom_step(aes(y = N), size = 1) +
  scale_x_reverse("Age (Mya)", limits = c(26, 0), breaks = seq(25, 0, by = -5)) + theme_1() +
  deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE)  +
  scale_y_continuous(limits = c(0, 40), breaks = seq(0, 40, 5)) + labs(y = "Species")

# Print to PDF ------------------------------------------------------------

(p4 | p2) /
(p5 | p3)

pdf(file = 'output/RPANDA_results_arid_sco.pdf', width = 11.33, height = 8.5)
(p4 | p2) /
  (p5 | p3)
dev.off()



pdf(file = "output/RPANDA_arid_sco.pdf")
arid_sco_plot
dev.off()





# Weird values ------------------------------------------------------------
# These models estimated lambda values that are very high > 100, which does not
# make sense. So plot to check.

load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_mean_sco_b.Rdata'); mean_sco <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_max_sco_b.Rdata'); max_sco <- Anilios_res
load('data/intermediate_data/diversification_analyses/Anilios_EnvDep_min_str_b.Rdata'); min_str <- Anilios_res

dev.off()
plot_fit_env_new(mean_sco[[1]]$BEnvVar_EXPO,
                 env_data = env_data_list$arid_sco, tot_time = tot_time, 
                 name_env = "Mean temp Scotese")

dev.off()
plot_fit_env_new(max_sco[[1]]$BEnvVar_EXPO,
                 env_data = env_data_list$arid_sco, tot_time = tot_time,
                 name_env = "Max temp Scotese")