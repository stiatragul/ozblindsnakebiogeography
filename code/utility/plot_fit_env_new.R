# plot_fit_env_new.R
# modified function of plot_fit_env from RPANDA.

library(RPANDA); library(dplyr)
library(pspline)

plot_fit_env_new <- function (fit.env, env_data, tot_time, name_env) 
{
  if (!inherits(fit.env, "fit.env")) 
    stop("object is not of class \"fit.env\"")
  t <- seq(0, tot_time, length.out = 100)
  # dev.new()
  plot(-t, sapply(t, fit.env$f.lamb), type = "l", xlab = "time", 
       ylab = "speciation rate", main = "Fitted speciation rate")
  df <- smooth.spline(env_data[, 1], env_data[, 2])$df
  spline_result <- sm.spline(env_data[, 1], env_data[, 2], 
                             df = df)
  env_func <- function(t) {
    predict(spline_result, t)
  }
  # dev.new()
  plot(env_func(t), sapply(t, fit.env$f.lamb), xlab = paste(name_env, "Environmental data", sep=" "), 
       ylab = "speciation rate", main = "Fitted speciation rate")
  if ("f.mu" %in% attributes(fit.env)$names) {
    # dev.new()
    plot(-t, sapply(t, fit.env$f.mu), type = "l", xlab = "time", 
         ylab = "extinction rate", main = "Fitted extinction rate")
    # dev.new()
    plot(env_func(t), sapply(t, fit.env$f.mu), xlab = paste(name_env, "Environmental data", sep=" "), 
         ylab = "extinction rate", main = "Fitted extinction rate")
    r <- function(t) {
      fit.env$f.lamb(t) - fit.env$f.mu(t)
    }
    # dev.new()
    plot(-t, sapply(t, r), type = "l", xlab = "time", ylab = "net diversification rate", 
         main = "Fitted net diversification rate")
    # dev.new()
    plot(env_func(t), sapply(t, r), xlab = paste(name_env, "Environmental data", sep=" "), 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
  }
  else {
    # dev.new()
    plot(-t, sapply(t, fit.env$f.lamb), type = "l", xlab = "time", 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
    # dev.new()
    plot(env_func(t), apply(t, fit.env$f.lamb), xlab = paste(name_env, "Environmental data", sep=" "), 
         ylab = "net diversification rate", main = "Fitted net diversification rate")
  }
}

# CUSTOM DF FUNCTION ------------------------------------------------------

fit_env_df_make <- function(env_data, tot_time, fit.env) {
  
  # Get time
  t <- seq(0, tot_time, length.out = 100)
  
  # Fitted speciation rate per time
  speciation_rate <- sapply(t, fit.env$f.lamb)
  
  # Fit spline to environmental data
  df <- smooth.spline(env_data[, 1], env_data[, 2])$df
  spline_result <- pspline::sm.spline(env_data[, 1], env_data[, 2], df = df)
  
  # Fit environment data
  env_func <- function(t) {
    predict(spline_result, t)
  }
  Environmental_data <- env_func(t)
  
  # If model has extinction calculate extinction
  if ("f.mu" %in% attributes(fit.env)$names) {
    
    # Fitted Extinction
    extinction_rate <- sapply(t, fit.env$f.mu)
    
    net_diversification_rate_calc <- function(t) {
      fit.env$f.lamb(t) - fit.env$f.mu(t)
    }
    
    # Fitted speciation - extinction
    net_div_rate <- sapply(t, net_diversification_rate_calc)
    
    # Returns a dataframe with values for easy plotting
    result_output1 <- data.frame(age = t, rev_age = -t, 
                                 speciation = speciation_rate,
                                 extinction = extinction_rate, 
                                 environment_data = Environmental_data)
    return(result_output1)
  }
  else {
    result_output2 <- data.frame(age = t, rev_age = -t,
                                 speciation = speciation_rate,
                                 environment_data = Environmental_data)
    return(output2)
  }
}



# Prep work env -----------------------------------------------------------

# PREP function to work
load(file = "data/intermediate_data/diversification_analyses/env_data_list.RData")
trees <- ape::read.tree(file = "data/intermediate_data/diversification_analyses/blindsnake.trees", tree.names = c("st", "b"))

fos_tree <- phytools::force.ultrametric(trees[[2]],"extend")
fos_tree$edge.length <- fos_tree$edge.length * 100
plot(fos_tree)

# Drop Anilios splendidus
fos_tree <- ape::drop.tip(phy = fos_tree, tip = "Anilios_splendidus")

# Crown age
tot_time <- max(node.age(fos_tree)$ages)

# Limit year to max total time
subset_env_age <- function(df) {
  df <- df %>% filter(time_mya < 26)
}

env_data_list <- lapply(env_data_list, subset_env_age)

## Fraction of diversity represented. In our tree we have 
## 48 Described species of Anilios + 2 more grypus undescribed + 1 more ligatus - A. splendidus
## So our "true" species at the moment is 50 and we have 41 tips in our tree

length(fos_tree$tip.label)
frac = 38/50