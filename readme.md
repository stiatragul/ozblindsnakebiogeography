# Data set for "Paleoenvironmental models for Australia and the impact of aridification on blind snake diversification"

---

These data sets were used to perform analyses included in the research paper "Paleoenvironmental models for Australia and the impact of aridification on blind snake diversification." 

Main aims for the project:

1. Estimate the historical biogeography of Australian blind snakes using 'BioGeoBEARS.'
1. Fit birth-death models and estimate diversification rates under different paleoenvironmental conditions using 'RPANDA.'
1. Compare diversification rates between arid-adapted and mesic-adapted lineages using 'ClaDS' and 'hisse.'

## Data and file structure

**/tree/** -- this folder contains the subset of phylogenies for *Anilios*. 

  - anilios_newick_b.tre - tree with *A. splendidus*
  - anilios_newick_st.tre - tree with *A. splendidus*
  - subset_anilios_newick_b.tre - tree without *A. splendidus*
  - subset_anilios_newick_st.tre - tree without *A. splendidus*

**/bears_txt/** -- this folder contains .txt files necessary for fitting biogeographical models in' BioGeoBEARS' 

  - geofile.txt - geographic range for each species
  - biome_distance.txt - modify distance for +x analysis

**/intermediate_data/bears/**
  - blindsnake_b.tre - phylogeny for fitting BioGeoBEARS models. 

**/paleo_env** -- this folder contains the reconstructed paleoenvironment data from three data sources as described in the paper. 

  - Australia_climate_data_40Mya_Valdes.csv - reconstructed from Valdes et al., (2021)
  - Australia_climate_data_45Mya_Scotese.csv - reconstructed from Scotese & Wright, (2018)
  - Australia_climate_data_45Mya_Straume.csv - reconstructed from Straume et al., (2020)

**/intermediate_data/diversification_analyses/** -- file for 'RPANDA' analyses
  - blindsnake.trees - contains two versions of *Anilios* trees. 

**/2021_ALA_blindsnake_occurence_data/** -- occurrence data downloaded from Atlas of the Living Australia (ALA)

  - Edited_ALA_2021_blindsnake_distribution.csv - occurrence data downloaded from ALA but edited by ST
  - citation.csv - citation information

**/intermediate_data/geohisse/**

  - arid_nonarid_both_states.csv - list of lineages and their geographic state. 0 = widespread, 1 = arid, and 2 = mesic. Designation of geographic states were based on distribution data and literature review. 
  - 2022_species_list_arid_nonarid_widespread.csv - includes information about which species are in the phylogeny to calculate fraction.

## Other information

Zenodo repository for "Data set: Australia's hidden radiation - phylogenomic analysis reveals rapid Miocene radiation of blind snakes"
  *  *Anilios* phylogeny can be found at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7155340.svg)](https://doi.org/10.5281/zenodo.7155340)

Occurrence data
  * Atlas of the Living Australia data - [DOI](https://doi.org/10.26197/ala.d92678b1-ad2d-437b-9457-9f52737ba003)

<!-- Github -->

## Code/Software

All scripts can be run using open source software.

  * R is required to run R scripts (.R).
  * Julia is required to run Julia scripts (.jl).

**/Code**

  - 00_*.R - scripts are used for preparing data for analyses
  - 01_fit_*.R - scripts were used to fit various RPANDA models. Note difference in initial lambda parameters for some models.
  - 01_env_data_plots.R - script to plot different paleoenvironmental data
  - 02_diversification_plots_b.R - script for plotting results
  - 02_diversification_plots_b.R - script for plotting results
  - 02_table_fit_env_results_b.R - script to summarise model fit
  - 03_BioGeoBEARS_analyses_parallel.R - script for fitting multiple BioGeoBEARS models 
  - 04_BioGeoBEARS_results_bsm.R - Conduct Biogeographic Stochastic Mapping for the best fitting model
  - 04_BioGeoBEARS_results_plots.R - plotting Biogeographic Stochastic Mapping from the best fitting model
  - 04_BioGeoBEARS_results_tables.R - Biogeographic Stochastic Mapping result table from the best fitting model
  - 05_GeoHiSSE_fit.R - fit GeoSSE and GeoHiSSE using 'hisse' package. 
  - 06_ClaDS2.jl - Julia script to estimate branch-specific rate under ClaDS
  - 06_ClaDS_plot_tips.R - plot results from ClaDS

## Contact

Should you have questions about these analysis scripts, please do not hesitate to contact Sarin Tiatragul (contact information can be found in the paper) or on Github (https://github.com/stiatragul)
