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

tip_rates <- geohisse_df %>% left_join(s_tree_df, by = c("taxon" = "label"))

tip_rates$state <- as.character(tip_rates$state)

# Non-arid species
# rates2 <- tip_rates %>% dplyr::filter(state == 2)


p1 <- tip_rates %>% 
  dplyr::filter(state != 2) %>% 
ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1, alpha = 0.4) +
  scale_x_continuous(limits = c(0, 0.2)) +
  scale_fill_manual(labels = c("Widespread", "Arid", "Non-arid"), 
                    values = c("blue", "red")) +
  theme_classic()

p2 <- tip_rates %>% 
  dplyr::filter(state == 2) %>% 
  ggplot(aes(x = branch.length, group = state, fill = state)) +
  geom_density(adjust = 1, alpha = 0.4) +
  scale_fill_manual(labels = c("Non-arid"), 
                    values = c("purple")) +
  theme_classic()

pdf(file = 'output/supp_ClaDS_rates.pdf', height = 9, width = 8.5)
p1 /
  p2
dev.off()
