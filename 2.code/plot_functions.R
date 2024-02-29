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
  
  x_max <- ppc_df %>% group_by(data_source) %>% summarise(maxi = max(obs_val)) %>% pull(maxi)
  x_min <- ppc_df %>% group_by(data_source) %>% summarise(mini = min(obs_val)) %>% pull(mini)
  y_max <- ppc_df %>% group_by(data_source) %>% summarise(maxi = max(sim_val)) %>% pull(maxi)
  y_min <- ppc_df %>% group_by(data_source) %>% summarise(mini = min(sim_val)) %>% pull(mini)
  
  pvalues_df <- tibble(label = compute_bayesian_pvalue(ppc.out),
                       data_source = datasets_names,
                       x = 0.1*(x_max - x_min) + x_min, y = 0.9*(y_max-y_min)+y_min)
  
  plot <- ggplot(data = ppc_df, aes(x=obs_val, y=sim_val, col = obs_val<sim_val)) +
    geom_point(shape = 20) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.75) +
    facet_wrap(~data_source, scales = "free") +
    geom_text(data = pvalues_df, mapping = aes(x = x, y = y, label = label), color = "black") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Posterior predictive check",
         subtitle = "Calculated for each dataset",
         x = "Observed values",
         y = "Simulated values")
  
  return(plot)
}
