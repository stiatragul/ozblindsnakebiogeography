# install_biogeobears.R

install.packages("rexpokit")
install.packages("cladoRcpp")


# Install BioGeoBEARS from Github
library(devtools)
devtools::install_github(repo="nmatzke/BioGeoBEARS")

# Check for previous install
find.package(package="BioGeoBEARS")

# Check version
packageVersion("rexpokit")
# ‘0.26.6.0.9001’ or higher
packageVersion("cladoRcpp")
# ‘0.15’ or higher
packageVersion("BioGeoBEARS")
# ‘1.1’ or higher
citation(package="BioGeoBEARS")


# SIMPLE BIOGEOBEARS SETUP
# ======================================================
# Install optimx
# install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")

# Also get snow (for parallel processing)
# install.packages("snow")
# library(snow)

# Install phylobase
# install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com")


