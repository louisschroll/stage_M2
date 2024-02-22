# HEADER ------------------------------------------------------------------------
#
# Script name:  plot_functions.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-16
#
# Script description:
#
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space


plot_PPC <- function(ppc.out, datasets_names=NULL){
  nb_datasets <- length(ppc.out$fit.y)
  if (is.null(datasets_names)){
    datasets_names = 1:nb_datasets
  }
  obs_values <- unlist(ppc.out$fit.y)
  sim_values <- unlist(ppc.out$fit.y.rep)
  
  nb_points <- length(obs_values)/nb_datasets
  
  ppc_df <- tibble(obs_val = obs_values, 
                   sim_val = sim_values, 
                   data_source = rep(datasets_names, each = nb_points))
  
  plot <- ggplot(data = ppc_df, aes(x=obs_val, y=sim_val, col = as.factor(data_source))) +
    facet_wrap(~data_source, scales = "free") +
    geom_point(shape = 20) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.75) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Posterior predictive check",
         subtitle = "Calculated for each dataset",
         x = "Observed values",
         y = "Simulated values")
  
  return(plot)
}