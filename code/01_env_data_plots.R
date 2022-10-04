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
    deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
    # scale_colour_manual(values = viridis(4), name = "Data source", labels = c("RPANDA Global", "Scotese Australia", "Straume Australia", "Valdes Australia")) +
    scale_y_continuous(limits = c(-5, 40), breaks = seq(-5, 40, 5)) + labs(y = "Mean temperature (°C)") +
    scale_x_reverse("Age (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() 
}


# pdf(file = "output/paleo_env_data_plots.pdf", height = 10, width = 10)

env_df[which(env_df$source %in% c("mean_sco", "max_sco", "min_sco", "global_temp")),] %>% env_plottr()

env_df[which(env_df$source %in% c("mean_str", "max_str", "min_str")),] %>% env_plottr()

env_df[which(env_df$source %in% c("mean_val", "max_val", "min_val")),] %>% env_plottr()

env_df[which(env_df$source %in% c("arid_sco", "arid_str", "arid_val")),] %>% 
  ggplot(data = ., aes(x = age, colour = source)) +
  geom_line(aes(y = vals), size = 1) +
  scale_x_reverse("Age (Ma)", limits = c(45, 0), breaks = seq(45, 0, by = -5)) + theme_bw() +
  deeptime::coord_geo(xlim = c(46, 0), dat = c("epochs"), pos = 'bottom', expand = TRUE) +
  labs(y = "Aridity index")

dev.off()  
