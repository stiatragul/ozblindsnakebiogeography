# 06_CLaDS_plot_tips.R

# Cladogenetic diversification in rate shifts 
# Estimates lineage specific diversification rates on a phylogeny. 
# Accounts for diverse sources of variation in diversification rates that occur
# during the evolutionary history of clades. 



# libraries ---------------------------------------------------------------
library(ape); library(phytools); library(RPANDA)
library(phyloch); data(strat2012); library(treeio)
library(patchwork)

# Load data ---------------------------------------------------------------
# Load ClaDS output from julia script code/06_ClaDS2.jl
load('data/intermediate_data/ClaDS/clads_output.Rdata')
geohisse_df <- read.csv('data/intermediate_data/geohisse/arid_nonarid_both_states.csv')

# Plotting ----------------------------------------------------------------
# Tree
s_tree <- CladsOutput$tree

# Shorten tip labels
s_tree$tip.label <- gsub(pattern = "Anilios", replacement = "A.", s_tree$tip.label)

source('code/utility/plot_ClaDS_phylo_colour.R') # custom plot function
source('code/utility/squish_trans.R') # squish continuous axis

# Plot branch-specific speciation rate
dev.off()
pdf(file = 'output/supp_ClaDS_magma.pdf', width = 8.5, height = 11.33)
plot_ClaDS_phylo_colour(s_tree, rates = CladsOutput$lambdai_map, 
                 show.tip.label = T, log = T)
axisGeo(GTS = strat2012, unit = "epoch"); axisPhylo()
dev.off()

pdf(file = 'output/fig_ClaDS_colour.pdf', width = 11.33, height = 8.5)
plot_ClaDS_phylo(s_tree, rates = CladsOutput$lambdai_map, 
                 show.tip.label = T, log = T)
axisGeo(GTS = strat2012, unit = "epoch"); axisPhylo()
dev.off()

# Plot tip specific diversification rate
pdf(file = 'output/fig_ClaDS_tip_div.pdf', width = 11.33, height = 8.5)
plot_ClaDS_phylo_colour(s_tree, rates = CladsOutput$lambdatip_map, 
                                        show.tip.label = T, log = T)
axisGeo(GTS = strat2012, unit = "epoch"); axisPhylo()
dev.off()

# Get tip rates -----------------------------------------------------------

# Make new phylo object
s_tree_rates <-  CladsOutput$tree

# add vector of speciation rates onto edge lengths of phylo (according to RPANDA doc, they are in the same order)
s_tree_rates$edge.length <- CladsOutput$lambdai_map

plot(s_tree_rates, show.tip.label = T, use.edge.length = F)
nodelabels()
# Write tree with ClaDS rates as edge lengths
# write.tree(s_tree_rates, file = "data/intermediate_data/ClaDS/final_tree_with_CLaDs_rates_as_edge_lengths.tree")

# Correlate with biome ----------------------------------------------------
 
s_tree_df <- s_tree_rates %>% 
  treeio::as_tibble() %>% 
  dplyr::filter(!is.na(label)) 

branch_specific_rates <- geohisse_df %>% left_join(s_tree_df, by = c("taxon" = "label"))

branch_specific_rates$state <- as.character(branch_specific_rates$state)

# Non-arid species
# rates2 <- tip_rates %>% dplyr::filter(state == 2)

p1 <- branch_specific_rates %>% 
  dplyr::filter(state != 2) %>% 
ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1, alpha = 0.4) +
  scale_x_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(labels = c("Widespread", "Arid", "Non-arid"), 
                    values = c("blue", "red")) +
  labs(x = "Branch-specific rates") +
  theme_classic()

p2 <- branch_specific_rates %>% 
  dplyr::filter(state == 2) %>% 
  ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1, alpha = 0.4) +
  scale_fill_manual(labels = c("Non-arid"), 
                    values = c("purple")) +
  labs(x = "Branch-specific rates") +
  theme_classic()


library(ggbreak)

p_all_bsr <- branch_specific_rates %>% 
  ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1 , alpha = 0.7) +
  scale_x_continuous(limits = c(-0.05, 6.35), expand = c(0, 0)) +
  scale_x_break(c(0.45, 2.29)) +
  scale_x_break(c(2.4, 4.68)) +
  scale_x_break(c(4.78, 5.7)) +
  scale_x_break(c(5.71, 6.27)) +
  scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  expand_limits(x = 0, y = 0) +
  scale_fill_manual(labels = c("Widespread", "Arid", "Non-arid"),
                    values = c("#5B5F97", "#FF6B6C", "#FFD400")) +
  labs(x = "Branch-specific rates", y = "Density") +
  theme_classic()

p_all_bsr

branch_specific_rates %>% 
  ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1.1, alpha = 0.4) +
  scale_x_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(labels = c("Widespread", "Arid", "Non-arid"), 
                    values = c("blue", "red", "purple")) +
  # scale_x_continuous(limits = c(0, 0.18), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Branch-specific rates") +
  theme_classic()

pdf(file = 'output/supp_ClaDS_rates.pdf', height = 10, width = 8.5)
p_all_bsr 
dev.off()



# Boxplot -----------------------------------------------------------------

branch_specific_rates %>% 
  ggplot(aes(x = state, y = branch.length, group = state, fill = state)) +
  # geom_density(adjust = 1.1, alpha = 0.4) +
  geom_boxplot() +
  scale_y_continuous(limits = c(0, 0.5)) +
  scale_fill_manual(labels = c("Widespread", "Arid", "Non-arid"), 
                    values = c("blue", "red", "purple")) +
  scale_x_discrete(limits = c("1", "2", "0")) + 
  labs(x = "Geographic states", y = "Branch-specific speciation rate (per Myr)") +
  theme_classic()


