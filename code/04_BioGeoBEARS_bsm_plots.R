# 04_BioGeoBEARS_bsm_plots.R
# Nov 2022
# Plotting figures from 50 BioGeoBEARS biostochastic maps


# Libraries ---------------------------------------------------------------
library(dplyr)
library(stringr)
library(patchwork)
library(corrplot)
library(GenSA);library(FD)      
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)
library(ggplot2)

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
# hist_event_counts(counts_list, pdffn=paste0("output/", model_name, "_histograms_of_event_counts.pdf"))

counts_list$all_dispersals_counts_fromto_means
counts_list$all_dispersals_counts_fromto_sds


## Largest outwards dispersal
rowSums(counts_list$all_dispersals_counts_fromto_means)
rowSums(counts_list$all_dispersals_counts_fromto_sds)

## Largest inwards dispersals

rowSums(t(counts_list$all_dispersals_counts_fromto_means))
rowSums(t(counts_list$all_dispersals_counts_fromto_sds))


# Mean transitions --------------------------------------------------------
###
### OUTWARDS DISPERSAL ###
####
mean_out_disp <- counts_list$all_dispersals_counts_fromto_means 
sd_out_disp <- counts_list$all_dispersals_counts_fromto_sds

### INWARDS DISPERSAL ###
# Opposite of outwards dispersal so just need to transpose matrix t()

# Transpose matrix for inwards transition
mean_in_disp <- as.data.frame(t(counts_list$all_dispersals_counts_fromto_means))
sd_in_disp <- t(counts_list$all_dispersals_counts_fromto_sds)


# change row names for both inwards and outwards df
mean_out_disp$from <- paste0("from", rownames(mean_in_disp))
mean_in_disp$to <- paste0("to", rownames(mean_in_disp))

# Pivot longer and prep for ggplot 
mean_out_disp_area <- mean_out_disp %>%
  tidyr::pivot_longer(!from, names_to = "to", values_to = "count") %>% 
  dplyr::mutate(state_from = ifelse(from == "fromC", "arid", "non_arid")) %>% 
  dplyr::mutate(state_to = ifelse(to == "C", "arid", "non_arid")) %>% 
  dplyr::mutate(subset_from = ifelse(from %in% c("fromG", "fromH", "fromI"), "island", from)) %>% 
  dplyr::mutate(subset_to = ifelse(to %in% c("G", "H", "I"), "tozisland", to))  

mean_out_disp_prop <- mean_out_disp_area %>% 
  dplyr::group_by(subset_from) %>% 
  dplyr::summarise(Freq = sum(count)) %>% 
  dplyr::mutate(Prop = 100*Freq/sum(Freq))

mean_in_disp_area <- mean_in_disp %>% 
  tidyr::pivot_longer(!to, names_to = "from", values_to = "count") %>% 
  dplyr::mutate(state_from = ifelse(from == "C", "arid", "non_arid")) %>% 
  dplyr::mutate(state_to = ifelse(to == "toC", "arid", "non_arid")) %>% 
  dplyr::mutate(subset_from = ifelse(from %in% c("G", "H", "I"), "island", from)) %>% 
  dplyr::mutate(subset_to = ifelse(to %in% c("toG", "toH", "toI"), "tozisland", to))

mean_in_disp_prop <- mean_in_disp_area %>% 
  dplyr::group_by(to) %>% 
  dplyr::summarise(Freq = sum(count)) %>% 
  dplyr::mutate(Prop = 100*Freq/sum(Freq))
  

# Plotting ----------------------------------------------------------------

# x-labels
biome_labs <- c("Trop Grass", "Temp Forest", "Arid", "Mediterranean", "Trop Forest", "Temp Grass", "Sunda", "N Guinea", "Bismarck Arc")
biome_values <- unique(mean_in_disp_prop$to)
biome_colours <- c("#dcc674", #0
                   "#1f4733", #1
                   "#e03f28", #2
                   "#4f937f", #3
                   "#179d48", #4 
                   "#5265ae", #5
                   "#d2e6ae", #6
                   "#92d1bd", #7
                   "#eaf400")

### OUTWARDS TRANSITION PLOT

# outwards dispersals 
outwards_disp <- ggplot(data = mean_out_disp_area, aes(fill= to, y = count, x = from)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Outwards dispersal", y = "Mean outwards transition frequency") + 
  scale_x_discrete(labels = biome_labs) + 
  scale_fill_manual(labels = biome_labs, values = biome_colours) +
  scale_y_continuous(limits = c(0, 17)) +
  guides(fill = guide_legend(title = "From")) 
  # theme(legend.position = "none")

# outwards dispersals geo state
outwards_disp_geo <- ggplot(data = mean_in_disp_area, aes(fill= state_to, y = count, x = state_from)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Outwards dispersal", y = "Mean outwards transition frequency") + 
  scale_fill_manual(labels = c('arid', 'non-arid'), values = c(biome_colours[3], "grey")) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "bottom")

# outwards dispersals grp
outwards_disp_grp <- ggplot(data = mean_in_disp_area, aes(y = count, fill= subset_to, x = subset_from)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Outwards dispersal", y = "Mean outwards transition frequency") + 
  scale_x_discrete(labels = c(biome_labs[1:6], 'islands')) +
  scale_fill_manual(values = c(biome_colours[1:6], 'grey')) +
  scale_y_continuous(limits = c(0, 17)) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none") +
  # geom_text(data = mean_out_disp_prop,
  #           aes(y = Freq + 1, x = subset_from,
  #               label = sprintf('%.2f%%', Prop)))

  geom_text(aes(label = count), position = position_stack(vjust = 0.5))

### Inward transitions

# inwards dispersals 
inwards_disp <- ggplot(data = mean_in_disp_area, aes(fill= from, y = count, x = to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inwards dispersal", y = "Mean inwards transition frequency") + 
  scale_x_discrete(labels = biome_labs) + 
  scale_fill_manual(labels = biome_labs, values = biome_colours) +
  scale_y_continuous(limits = c(0, 17)) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none")

# geom_text(data = df2,
#           aes(y = Freq + 1,
#               label = sprintf('%.2f%%', Prop)))# inwards dispersals geo state
inwards_disp_geo <- ggplot(data = mean_in_disp_area, aes(fill= state_from, y = count, x = state_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inwards dispersal", y = "Mean inwards transition frequency") + 
  scale_fill_manual(labels = c('arid', 'non-arid'), values = c(biome_colours[3], "grey")) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "bottom")

# inwards dispersals grp
inwards_disp_grp <- ggplot(data = mean_in_disp_area, aes(fill= subset_from, y = count, x = subset_to)) + 
  geom_bar(position="stack", stat="identity") +
  labs(x = "Inwards dispersal", y = "Mean inwards transition frequency") + 
  scale_x_discrete(labels = c(biome_labs[1:6], 'islands')) +
  scale_fill_manual(values = c(biome_colours[1:6], 'grey')) +
  scale_y_continuous(limits = c(0, 17)) +
  guides(fill = guide_legend(title = "From")) +
  theme(legend.position = "none") 
  # geom_text(aes(label = count), position = position_stack(vjust = 0.5))



# SAVE TO PDF -------------------------------------------------------------

### stitch together with patchwork package

# Save PDF to be included as inset in BioGeoBEARS plot
# pdf(file = 'output/biogeobears_dispersal_freq.pdf', width = 11.5, height = 5.7)
pdf(file = 'output/biogeobears_dispersal_mean.pdf', width = 11.5, height = 5.7)
outwards_disp_grp +inwards_disp_grp
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
