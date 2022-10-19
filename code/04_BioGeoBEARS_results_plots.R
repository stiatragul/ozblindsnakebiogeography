# 04_BioGeoBEARS_results_plots.R
# August 2022
# Code adapted from Carlos Pavon Vazquez

# Test best model fit from BioGeoBEARS analyses

# load results from analysis script
# source("code/03_BioGeoBEARS_analyses_parallel.R") or load from saved intermediate data same thing.


# Libraries ---------------------------------------------------------------

library(GenSA);library(FD)      
library(rexpokit);library(cladoRcpp); library(BioGeoBEARS)
# library(dplyr)
library(parallel); library(snow)


# Tree --------------------------------------------------------------------

tree.c <- ape::read.tree('data/intermediate_data/bears/blindsnake_b.tre')
tree.c$edge.length <- tree.c$edge.length * 100

# Saved results -----------------------------------------------------------

load("data/intermediate_data/bears/bears_model_fit.Rdata")


#######################################################
# Plot ancestral states - DEC+j+x
#######################################################
resjx$inputs$states_list

colour_lists <- data.frame(biome = c("A", "B", "C", "D", "E", "F", "G", "H", "I"),
                           colours = c("#dcc674", #0
                                       "#1f4733", #1
                                       "#e03f28", #2
                                       "#4f937f", #3
                                       "#179d48", #4 
                                       "#5265ae", #5
                                       "#d2e6ae", #6
                                       "#92d1bd", #7
                                       "#eaf400"))  #8 

# mix colours in meyerweb.com color blender

colors_list_for_states  = c("gray87","#dcc674","#1f4733","#e03f28","#4f937f",
                            "#dd9f47","#5265ae","#d2e6ae",
                            "#92d1bd","#d2e6ae", 
                            "#758151", #A and B
                            "#DE8951", #A and C
                            "#8FAA7A", #A and D
                            "#DDB15B", #A and E
                            "#848899", #A and F
                            "#88432D", #B and C
                            "#1B763E", #B and D
                            "#1B763E", #B and E
                            "#3B5776", #B and F
                            "#916D57", #C and D
                            "#935471", #C and F
                            "#9E6550", #D and C
                            "#517A99", #D and F
                            "#gray87", #E and A
                            "#DFEE4F", #G and I
                            "#BAE167", #H and I
                            rep("gray87", 46-26))


# Modify to match best fitting model
analysis_titletxt ="DEC+j+x results"

# Setup
results_object = resjx
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
tr=tree.c

dev.off()
# States

# pdf('output/bgb_letter_ancestral.pdf', width = 11.5, height = 8.5)
resstates = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j","x"), plotwhat="text", label.offset=0.45,
                                     tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, 
                                     include_null_range=TRUE, tr=tr,
                                     colors_list_for_states=colors_list_for_states,
                                     tipranges=tipranges)
tiplabels(tr$tip.label, adj = c(0,0.5), bg = NULL, col = NULL, frame = "none")

dev.off()



# pdf('output/pie_ancestral.pdf')
# Pie chart
BioGeoBEARS::plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j","x"), plotwhat="pie", 
                                      label.offset=0.45, tipcex=0.02, statecex=0.4, splitcex=0.4, titlecex=0.8, plotsplits=F, 
                                      cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges, plotlegend = T,
                                      # colors_list_for_states=colors_list_for_states,
                                      legend_cex = 0.5)
tiplabels(tr$tip.label, adj = c(0,0.5), bg = NULL, col = NULL, frame = "none")

dev.off()


# Results write up --------------------------------------------------------

plot(tr)
nodelabels()

