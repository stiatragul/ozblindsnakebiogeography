# 04_BioGeoBEARS_bsm_plots.R
# Nov 2022
# Plotting figures from 50 BioGeoBEARS biostochastic maps


# Libraries ---------------------------------------------------------------
library(dplyr)
library(stringr)
library(patchwork)
library(corrplot)

# Load data ---------------------------------------------------------------
load('data/intermediate_data/bears/BSB_50.Rdata')

# EXTRACT BSM OUTPUT ------------------------------------------------------

# For DEC+j+x 
## CLADOGENETIC Events -- sympatry (subset) + vicariance (range splitting)
## ANAGETIC Events -- Founder (jump j parameter) + extinction

# Range-expansion dispersal (all observed 'd' dispersals): controlled by the parameter 'd', which is freely estimated in DEC and related models.
                                                         # If the d parameter is 0, these counts are zero, but otherwise there should be a few events.

# Anagenetic dispersal (mean of all observed anagenetic 'a' or 'd' dispersals): since a = 0, it's just d
# Cladogenetic dispersal (mean of all observed jump 'j' dispersals): founder event 'j '

# ALL dispersal (mean of all observed anagenetic 'a', 'd' dispersals, PLUS cladogenetic founder/jump dispersal):
summary_counts_BSMs = counts_list$summary_counts_BSMs

# Summary of BSM
conditional_format_table(summary_counts_BSMs)

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0("output/", model_name, "_histograms_of_event_counts.pdf"))

# Count individual stuff --------------------------------------------------

### Sum up dispersal event to different locations ###

# Change row and column names to be more informative (from and to)
dimnames(counts_list$anagenetic_dispersals_counts_cube) <- list(paste0('from', LETTERS[1:9]), paste0('to', LETTERS[1:9]), 1:50)

# Find sum of each column (to get gross number of dispersal TO each biome per stochastic map)
sum_disp_area <- rowSums(counts_list$anagenetic_dispersals_counts_cube, dims = 2)

counts_list$anagenetic_dispersals_counts_cube

rownames(sum_disp_area) <- LETTERS[1:9]
sum_disp_area_df <- data.frame(sum_disp_area)
sum_disp_area_df$from <- rownames(sum_disp_area_df)

outward_disp_df <- sum_disp_area_df %>% 
  tidyr::pivot_longer(!from, names_to = "to", values_to = "count") %>% 
  dplyr::mutate(state_from = ifelse(from == "C", "arid", "non_arid")) %>% 
  dplyr::mutate(state_to = ifelse(to == "toC", "arid", "non_arid")) %>% 
  dplyr::mutate(subset_from = ifelse(from %in% c("G", "H", "I"), "island", from)) %>% 
  dplyr::mutate(subset_to = ifelse(to %in% c("toG", "toH", "toI"), "tozisland", to))

# x-labels
biome_labs <- c("Trop Grass", "Temp Forest", "Arid", "Mediterranean", "Trop Forest", "Temp Grass", "Sunda", "N Guinea", "Bismarck Arc")
biome_values <- unique(outward_disp_df$to)
biome_colours <- c("#dcc674", #0
                   "#1f4733", #1
                   "#e03f28", #2
                   "#4f937f", #3
                   "#179d48", #4 
                   "#5265ae", #5
                   "#d2e6ae", #6
                   "#92d1bd", #7
                   "#eaf400")

# Mean transitions
mean_in_disp <- counts_list$all_dispersals_counts_fromto_means 

# change row names
mean_in_disp$to <- paste0("to", rownames(mean_in_disp))


mean_in_disp_area <- mean_in_disp %>% 
  tidyr::pivot_longer(!to, names_to = "from", values_to = "count") %>% 
  dplyr::mutate(state_from = ifelse(from == "C", "arid", "non_arid")) %>% 
  dplyr::mutate(state_to = ifelse(to == "toC", "arid", "non_arid")) %>% 
  dplyr::mutate(subset_from = ifelse(from %in% c("G", "H", "I"), "island", from)) %>% 
  dplyr::mutate(subset_to = ifelse(to %in% c("toG", "toH", "toI"), "tozisland", to))

# outwards dispersals 
outwards_disp <- ggplot(data = inward_disp_area, aes(fill= from, y = count, x = to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_x_discrete(labels = biome_labs) + 
  scale_fill_manual(labels = biome_labs, values = biome_colours) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none")

# outwards dispersals geo state
outwards_disp_geo <- ggplot(data = inward_disp_area, aes(fill= state_from, y = count, x = state_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_fill_manual(labels = c('arid', 'non-arid'), values = c(biome_colours[3], "grey")) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "bottom")

# outwards dispersals grp
outwards_disp_grp <- ggplot(data = inward_disp_area, aes(fill= subset_from, y = count, x = subset_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_x_discrete(labels = c(biome_labs[1:6], 'islands')) +
  scale_fill_manual(values = c(biome_colours[1:6], 'grey')) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none")







# RAW FREQUENCY FROM 50 BSM -----------------------------------------------

# Get average score instead
outward_disp_df$count <- outward_disp_df$count/50

# plot outward dispersal from - > to
outward_dispersal <- ggplot(data = outward_disp_df, aes(fill= to, y = count, x = from)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Outward dispersals", y = "frequency") + 
  scale_x_discrete(labels = biome_labs) + 
  scale_fill_manual(labels = biome_labs, values = biome_colours) +
  guides(fill = guide_legend(title = "Destination/From")) +
  theme(legend.position = "top")

# Outward dispersal by geographic state
outward_dispersal_geo <- ggplot(data = outward_disp_df, aes(fill= state_to, y = count, x = state_from)) + 
  geom_bar(position="stack", stat="identity") + 
  labs(x = "Outward dispersals", y = "frequency") + 
  scale_fill_manual(labels = c('arid', 'non-arid'), values = c(biome_colours[3], "grey")) +
  guides(fill = guide_legend(title = "Destination")) +
  theme(legend.position = "top")

# Island grouped
outward_dispersal_grp <- ggplot(data = outward_disp_df, aes(fill= subset_to, y = count, x = subset_from)) + 
  geom_bar(position="stack", stat="identity") + 
  labs(x = "Outward dispersals", y = "frequency") + 
  scale_x_discrete(labels = c(biome_labs[1:6], 'islands')) + 
  scale_fill_manual(values = c(biome_colours[1:6], 'grey')) +
  guides(fill = guide_legend(title = "Destination")) +
  theme(legend.position = "top")

outward_disp_df %>% group_by(from) %>% 
  dplyr::summarise(mean = mean(count))

outward_disp_df %>% group_by(from) %>% 
  dplyr::summarise(total = sum(count))

# # Inward dispersal

# # transpose the matrix

counts_list$all_dispersals_counts_fromto_means
inward_disp_area_df <- t(counts_list$all_dispersals_counts_fromto_means) %>% as.data.frame()

# # t_ana_disp_list <- lapply(ana_disp_list, FUN = t)
# # # Find the sum of BSM from transposed 
# inward_disp_area_df <- Reduce('+', t_ana_disp_list) %>% data.frame() 
inward_disp_area_df$to <- paste0("to", rownames(inward_disp_area_df))

inward_disp_area <- inward_disp_area_df %>% 
  tidyr::pivot_longer(!to, names_to = "from", values_to = "count") %>% 
  dplyr::mutate(state_from = ifelse(from == "C", "arid", "non_arid")) %>% 
  dplyr::mutate(state_to = ifelse(to == "toC", "arid", "non_arid")) %>% 
  dplyr::mutate(subset_from = ifelse(from %in% c("G", "H", "I"), "island", from)) %>% 
  dplyr::mutate(subset_to = ifelse(to %in% c("toG", "toH", "toI"), "tozisland", to))

# inwards dispersals 
inwards_dispersal <- ggplot(data = inward_disp_area, aes(fill= from, y = count, x = to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_x_discrete(labels = biome_labs) + 
  scale_fill_manual(labels = biome_labs, values = biome_colours) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none")

# inwards dispersals geo state
inwards_dispersal_geo <- ggplot(data = inward_disp_area, aes(fill= state_from, y = count, x = state_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_fill_manual(labels = c('arid', 'non-arid'), values = c(biome_colours[3], "grey")) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "bottom")

# inwards dispersals grp
inwards_dispersal_grp <- ggplot(data = inward_disp_area, aes(fill= subset_from, y = count, x = subset_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inward dispersal", y = "frequency") + 
  scale_x_discrete(labels = c(biome_labs[1:6], 'islands')) +
  scale_fill_manual(values = c(biome_colours[1:6], 'grey')) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none")

### stitch together with patchwork package

# Save PDF to be included as inset in BioGeoBEARS plot
pdf(file = 'output/biogeobears_dispersal_freq.pdf', width = 11.5, height = 5.7)
outward_dispersal_grp +inwards_dispersal_grp
dev.off()

# Save PDF for supplementary
pdf(file = 'output/supp_biogeobears_dispersal.pdf', width = 11.33, height = 8.5)
outward_dispersal + outward_dispersal_geo + scale_y_continuous(limits = c(0, 1400)) +
  inwards_dispersal + inwards_dispersal_geo +scale_y_continuous(limits = c(0, 1400)) +
  patchwork::plot_layout(widths = c(2,1))
dev.off()

# Heat map of transitions -------------------------------------------------

# Combine island transitions into  one
means_dispersal <- as.matrix(counts_list$all_dispersals_counts_fromto_means)
means_dispersal_df <- counts_list$all_dispersals_counts_fromto_means

means_disp_df <- means_dispersal_df[1:6, 1:6]
means_disp_df$island <- c(0.12,0,0,0,0,0)
means_disp_df <- rbind(means_disp_df, c(0.04,0,0,0,0,0,0.1))

rownames(means_disp_df) <- c(LETTERS[1:6], 'island')

.M <- as.matrix(means_disp_df)

dev.off()

pdf(file = 'output/supp_biogeobears_biome_transitions.pdf')
p_disp_island <- corrplot(.M, is.corr = FALSE, addCoef.col = 'black', col = COL2('RdBu', 100),type = 'full') 
p_disp_all <- corrplot(means_dispersal, is.corr = FALSE, col = COL2('RdBu', 100), type = 'full')
dev.off()

# StDEV
counts_list$all_dispersals_counts_fromto_sds %>% tibble()




dev.off()



## STATs ##
# Check individual data
clado_events_tables[2][[1]] %>% 
  dplyr::filter(stringr::str_detect(clado_event_type, pattern = regex("[A-z]+"))) %>% 
  dplyr::select(clado_event_type, node, ord_ndname, clado_event_txt, clado_dispersal_to, everything()) %>% 
  tibble()

# Dataframes for each stochastic map labeled

clado_events_tables = BSMs_w_sourceAreas$clado_events_tables # resets

clado_events_tables[1][[1]] %>% class()

for(i in 1:length(clado_events_tables)){
  
  clado_events_tables[i][[1]]$stochastic_no <- paste0("s", i) 
  
}

# merge the n number of dataframes together

sm_df <- do.call(rbind, clado_events_tables)

# Founder events
founder_events_df <- sm_df %>%
  # filter for clado_events we care about
  dplyr::filter(str_detect(clado_event_type, pattern = regex("(j)"))) %>% 
  # select columns that we want 
  dplyr::select(clado_event_type, node, ord_ndname, clado_event_txt, clado_dispersal_to, stochastic_no, everything()) %>% 
  dplyr::mutate(from = str_extract(clado_event_txt, pattern = "[A-Z]")) %>% 
  dplyr::mutate(transitions = paste0(from, "-", clado_dispersal_to)) %>% 
  dplyr::group_by(node, transitions) %>% 
  dplyr::summarise(n = n()) %>% 
  dplyr::mutate(freq = n / sum(n)) 

founder_events_df$node <- as.character(founder_events_df$node)

# Get rid of those with one transition

founder_events_df


plot(tr)
nodelabels()

obj <- founder_events_df %>% 
  dplyr::select(node, transitions, freq) %>% 
  tidyr::pivot_wider(names_from = transitions, values_from = freq) 

# replace NA with 0
obj[is.na(obj)] <- 0

obj[-1]

pie_matrix <- as.matrix(obj[-1])
rownames(pie_matrix) <- obj$node

pie_matrix

dim(pie_matrix)[2]

library(viridis)
cols<-setNames(viridis(dim(pie_matrix)[2]), colnames(pie_matrix))

dev.off()
plot(tr)
nodelabels(pie=pie_matrix, piecol=cols,cex=0.5)
legend(x="bottomright", legend=colnames(pie_matrix), 
       pch=21,pt.cex=2,pt.bg=cols,cex=1)

pie_df <- founder_events_df %>% 
  dplyr::select(node, transitions, freq) %>% 
  tidyr::pivot_wider(names_from = transitions, values_from = freq) %>% 
  tidyr::pivot_longer(!node, names_to = "props") %>% 
  mutate(ypos = cumsum(value)- 0.5*value)

pie_df

library(ggplot2)

pie_charts_order_proportion <- ggplot(data=pie_df, aes(x=" ", y=value, group=node, colour=props, fill=props)) +
  geom_bar(width = 1, stat = "identity", colour = "white") +
  # scale_fill_manual("Proportions", values=c("#648FFF", "#FFB000")) +
  geom_text(aes(y = ypos, label = props), position = position_stack(vjust = 0.5), color = "white")  +
  coord_polar("y") + 
  facet_wrap(~ node) +theme_void() 

pie_charts_order_proportion 
