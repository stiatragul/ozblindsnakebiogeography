# name: 01_env_data_plots.R

# libraries ---------------------------------------------------------------

library(deeptime) # for deeptime axis
library(ggplot2); library(viridis)

# load data
load(file = "data/intermediate_data/diversification_analyses/env_data_list.RData")

# paleoclimate plots ------------------------------------------------------

str(env_data_list)

# Combine all data into one dataframe to plot using GGPLOT2. Need to make dataframe
# Make new column for each data set to have name of data as a column

for(i in 1:length(env_data_list)){
  env_data_list[[i]]$source <- names(env_data_list)[i]
}

# New column names so we can rbind them AGE, VAL(UES), SOURCE
env_data_list <- lapply(env_data_list, setNames, nm = c("age", "vals", "source"))

# Merge as df
env_df <- do.call(rbind, env_data_list)
rownames(env_df) <- NULL

env_df

# Function to plot
env_plottr <- function(.df){
  .df %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
    geom_line(aes(y = vals), size = 1) +
    # geom_point(aes(y = vals)) + 
    deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
    # scale_colour_manual(values = viridis(4), name = "Data source", labels = c("RPANDA Global", "Scotese Australia", "Straume Australia", "Valdes Australia")) +
    scale_y_continuous(limits = c(-5, 30), breaks = seq(-5, 30, 5)) + labs(y = "Mean temperature (°C)") +
    scale_x_reverse("Time (Ma)", limits = c(25, 0), breaks = seq(25, 0, by = -5)) + theme_bw() 
}
# env_plottr <- function(.df){
#   .df %>% 
#   ggplot(data = ., aes(x = age, colour = source)) +
#     geom_line(aes(y = vals), size = 1) +
#     geom_point(aes(y = vals)) + 
#     deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
#     # scale_colour_manual(values = viridis(4), name = "Data source", labels = c("RPANDA Global", "Scotese Australia", "Straume Australia", "Valdes Australia")) +
#     scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, 5)) + labs(y = "Mean temperature (°C)") +
#     scale_x_reverse("Age (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() 
# }

# pdf(file = "output/paleo_env_data_plots.pdf", height = 10, width = 10)

env_df[which(env_df$source %in% c("mean_sco", "max_sco", "min_sco", "global_temp")),] %>% env_plottr()

env_df[which(env_df$source %in% c("mean_str", "max_str", "min_str")),] %>% env_plottr()

env_df[which(env_df$source %in% c("mean_val", "max_val", "min_val")),] %>% env_plottr()

dev.off()  

# temp_plots <- env_df[which(env_df$source %in% c("mean_sco", "mean_str", "mean_val", "global_temp")),] %>% env_plottr()


temp_plots <- env_df[which(env_df$source %in% c("mean_sco", "mean_str", "mean_val", "global_temp")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  # geom_point(aes(y = vals)) + 
  deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  scale_colour_manual(values = c("#F0027F", "#386CB0", "#FDC086", "#666666"), name = "Data source", labels = c("RPANDA Global", "Scotese and Wright 2018", "Straume et al., 2020", "Valdes et al., 2021")) +
  # scale_color_brewer(type='qual',palette='Set1', name = "Data source", labels = c("RPANDA Global", "Scotese and Wright 2018", "Straume et al. 2019", "Valdes et al. 2021")) +
  scale_y_continuous(limits = c(-5, 30), breaks = seq(-5, 30, 5)) + labs(y = "Mean temperature (°C)") +
  scale_x_reverse("Time (Ma)", limits = c(25, 0), breaks = seq(25, 0, by = -5)) + theme_bw() +
  theme(panel.border = element_blank(), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  theme(legend.position = "top")


arid_plots <-env_df[which(env_df$source %in% c("arid_sco", "arid_str", "arid_val")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  scale_x_reverse("Time (Ma)", limits = c(25, 0), breaks = seq(25, 0, by = -5)) + theme_bw() +
  deeptime::coord_geo(xlim = c(26, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  labs(y = "Aridity index") +
  # scale_color_brewer(type='qual',palette='Set1') +
  # scale_color_brewer(type='qual',palette='Set1', name = "Data source", labels = c("Scotese and Wright 2018", "Straume et al. 2019", "Valdes et al. 2021")) +
  scale_colour_manual(values = c("#386CB0", "#FDC086", "#666666"), name = "Data source", labels = c("Scotese and Wright 2018", "Straume et al., 2020", "Valdes et al., 2021")) +
  theme(panel.border = element_blank(), axis.line.x = element_line(color = "black", size = 0.5), axis.line.y = element_line(color = "black", size = 0.5)) +
  theme(legend.position = "none")

library(patchwork)

# pdf(file = "output/paleo_data_compare_plots.pdf", height = 11.33, width = 8.5)
pdf(file = "output/paleo_data_compare_plots.pdf", height = 7, width = 4.33)
temp_plots /
  arid_plots
dev.off()



p_sco <- env_df[which(env_df$source %in% c("mean_sco", "max_sco", "min_sco")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  geom_point(aes(y = vals)) +
  deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  scale_colour_manual(values = viridis(3), name = "Scotese & Wright, 2018", labels = c("Max. temp", "Meann. temp", "Min. temp")) +
  scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, 5)) + labs(y = "Temperature (°C)") +
  scale_x_reverse("Time (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() +
  theme(legend.position = "top")


p_str <- env_df[which(env_df$source %in% c("mean_str", "max_str", "min_str")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  geom_point(aes(y = vals)) +
  deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  scale_colour_manual(values = viridis(3), name = "Straume et al., 2020", labels = c("Max. temp", "Meann. temp", "Min. temp")) +
  scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, 5)) + labs(y = "Temperature (°C)") +
  scale_x_reverse("Time (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() +
  theme(legend.position = "top")

p_val <- env_df[which(env_df$source %in% c("mean_val", "max_val", "min_val")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  geom_point(aes(y = vals)) +
  deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  scale_colour_manual(values = viridis(3), name = "Valdes et al., 2021", labels = c("Max. temp", "Meann. temp", "Min. temp")) +
  scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, 5)) + labs(y = "Temperature (°C)") +
  scale_x_reverse("Time (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() +
  theme(legend.position = "top")


pdf(file = "output/supp_paleo_temp.pdf", height = 7.02, width = 16.86)
p_sco + p_str + p_val
dev.off()



