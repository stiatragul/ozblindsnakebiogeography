# 06_ClaDS2.jl
# Species-specific diversification rates -- Implementation of ClaDS2 on Julia. 

# Remember to activate environment and install required packages
# Pkg.add("PANDA") # not Pandas because that's the python pandas. 

using PANDA

# Subset tree that includes ligatus and grypus OTUs
my_tree = load_tree("./data/tree/subset_anilios_newick_b.tre")
my_tree_st = load_tree("./data/tree/subset_anilios_newick_st.tre")


# Run ClaDS2 
output = infer_ClaDS(my_tree, f = 0.76, print_state = 100)
output_st = infer_ClaDS(my_tree_st, f = 0.76, print_state = 100)

# Save ClaDS run to disk.
using JLD2, FileIO
# @save "C:/Users/ST/Documents/repo/blindsnakebiogeography/data/intermediate_data/ClaDS/julia_ClaDS.jld2" output

# Load results
@load "./data/intermediate_data/ClaDS/julia_ClaDS.jld2" output

# Save as R object
save_ClaDS_in_R(output, "./data/intermediate_data/ClaDS/clads_output.Rdata")
save_ClaDS_in_R(output_st, "./data/intermediate_data/ClaDS/clads_output_st.Rdata")

# Results from the output

# Can check tip specific diversification rate
tip_rate(output, "Anilios_affinis")
tip_rate(output_st, "Anilios_silvia")


# Plotting
plot_CladsOutput(output, options = "lwd = 2")
plot_CladsOutput(output_st, options = "lwd = 2")

# Diversity through time
plot_CladsOutput(output, method = "DTT")

# Rate through time
plot_CladsOutput(output, method = "RTT")

# Marginal posterior densities
plot_CladsOutput(output, method = "chain")


## sigma = heterogeneity parameters
plot_CladsOutput(output, method = "density", id_par = "sigma")

## trend parameter alpha
plot_CladsOutput(output, method = "density", id_par = "alpha")

## turnover rate epsilon
plot_CladsOutput(output, method = "density", id_par = "epsilon")

## mean change in rate (m = alpha x epsilon ^ sigma^2 / 2)
plot_CladsOutput(output, method = "density", id_par = "lambda0")



## Rate through time and diversity through time
plot_CladsOutput(output, method = "density", id_par = "rate_12")
plot_CladsOutput(output, method = "density", id_par = "div_33")


