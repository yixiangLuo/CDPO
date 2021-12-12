#!/usr/bin/env Rscript

library(here)

source(here("R", "data.R"))
source(here("R", "ATE_methods.R"))
source(here("R", "mod_knockoff.R"))

source(here("expr", "experiment.R"))
source(here("expr", "utils.R"))
source(here("expr", "plot.R"))


# source(here("expr", "settings", "test.R"))


args <- commandArgs(trailingOnly=TRUE)
if (length(args) > 1) {
  stop("One and only one experiment at a time.")
} else if (length(args) == 1) {
  source(here("expr", "settings", args[1]))
}

file_name <- here("data", "temp", paste0("progress-", experiment, ".txt"))
if (file.exists(file_name)) {
  file.remove(file_name)
  invisible()
}


results <- lapply(1:length(X_types), function(iter_i){
  X_type <- X_types[iter_i]
  X_title <- X_titles[iter_i]
  
  n_level <- n_levels[iter_i]
  levels <- rep(n_level, p)
  levels[1 : (n_level * floor(p/n_level))] <- rep(1:n_level, each = floor(p/n_level))
  
  result <- get_fdp_power(n, p, levels, pi1, sig_sign, noise, alpha,
                          fig_x_var, sample_size, n_cores,
                          experiment, X_title, X_seed = iter_i+1000)
  
  save(result, alphas, fig_x_var,
       file = here("data", "temp", paste0(experiment, "-", X_type, ".Rdata")))
  
  return(result)
})

names(results) <- X_types

save(results, alphas, fig_x_var, file = here("data", paste0(experiment, ".Rdata")))
for(X_type in X_types){
  file.remove(here("data", "temp", paste0(experiment, "-", X_type, ".Rdata")))
}

method_names <- c("Hajek", "DR", "Match", "knockoff")
method_colors <- c("#984ea3", "dodgerblue3", "red", "#333333")
method_shapes <- rep(19, 4)

draw_fdp_power_curve(experiment, X_types, sample_size,
                     method_names, method_colors, method_shapes,
                     error_bar = F)

