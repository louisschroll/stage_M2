# HEADER -----------------------------------------------------------------------
#
# Script name:  plot_functions.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-16
#
# Script description:
#
#
# ------------------------------------------------------------------------------

# ---- Functions to plot Posterior Predicitve Checks ----
plot_PPC <- function(ppc.out, datasets_names=NULL, test = NULL, group = NULL){
  nb_datasets <- length(ppc.out$fit.y)
  if (is.null(datasets_names)){
    datasets_names = 1:nb_datasets
  }
  
  obs_values <- unlist(ppc.out$fit.y)
  sim_values <- unlist(ppc.out$fit.y.rep)
  
  nb_points <- length(obs_values)/nb_datasets
  
  ppc_df <- tibble(obs_val = obs_values, 
                   sim_val = sim_values, 
                   data_source = rep(as.factor(datasets_names), each = nb_points))
  
  x_max <- ppc_df %>% group_by(data_source) %>% summarise(maxi = max(obs_val)) %>% pull(maxi)
  x_min <- ppc_df %>% group_by(data_source) %>% summarise(mini = min(obs_val)) %>% pull(mini)
  y_max <- ppc_df %>% group_by(data_source) %>% summarise(maxi = max(sim_val)) %>% pull(maxi)
  y_min <- ppc_df %>% group_by(data_source) %>% summarise(mini = min(sim_val)) %>% pull(mini)
  data_source <- ppc_df %>% group_by(data_source) %>% summarise(maxi = max(obs_val)) %>% pull(data_source)
  
  pvalues_df <- tibble(label = compute_bayesian_pvalue(ppc.out),
                       data_source = as.factor(datasets_names)) %>% 
    full_join(tibble(data_source = data_source,
                     x = 0.1*(x_max - x_min) + x_min, y = 0.9*(y_max-y_min)+y_min),
              by = join_by(data_source))
                       
  
  plot <- ggplot(data = ppc_df, aes(x=obs_val, y=sim_val, col = obs_val<sim_val)) +
    geom_point(shape = 20) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.75) +
    facet_wrap(~data_source, scales = "free") +
    geom_text(data = pvalues_df, mapping = aes(x = x, y = y, label = label), color = "black") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste0("Posterior predictive check (", test," ", group, ")"),
         subtitle = "Calculated for each dataset",
         x = "Observed values",
         y = "Simulated values")
  
  return(plot)
}


# ---- Functions to plot the trace and density of parameters ----
make_coeff_plot <- function(model_result, param_to_plot = NULL, ncol = 1){
  nchains <- model_result$n.chains
  niter <- (model_result$n.samples - model_result$n.burn)/model_result$n.thin
  
  beta_df <- prepare_coefficients_df(model_result$beta.samples, nchains, niter)
  alpha_df <- prepare_coefficients_df(model_result$alpha.samples, nchains, niter)
  
  if (is.null(param_to_plot))
    param_to_plot <- c(unique(beta_df$param), unique(alpha_df$param))
  
  coeff_plots <- list(
    alpha_traceplot = make_traceplot(alpha_df, param_to_plot, ncol),
    alpha_density = make_density_plot(alpha_df, param_to_plot, ncol),
    beta_traceplot = make_traceplot(beta_df, param_to_plot, ncol),
    beta_density = make_density_plot(beta_df, param_to_plot, ncol)
  )
  return(coeff_plots)
}


prepare_coefficients_df <- function(coeff_df, nchains, niter){
  coeff_df %>%  
    janitor::clean_names() %>% 
    as_tibble() %>% 
    add_column(iteration = rep(1:niter, nchains)) %>%
    add_column(chain = rep(1:nchains, each = nrow(.)/nchains)) %>%
    pivot_longer(c(-iteration, -chain), values_to = "value", names_to = "param")
}


make_traceplot <- function(clean_coeff_df, param_to_plot, ncol = 1){
  clean_coeff_df %>% 
    filter(param %in% param_to_plot) %>% 
    ggplot(aes(x = iteration, y = value, color = as.factor(chain))) +
    facet_wrap(~param, scales = "free", ncol = ncol) +
    geom_line(alpha = 0.75) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Traceplot")
}


make_density_plot <- function(clean_coeff_df, param_to_plot, ncol = 1){
  clean_coeff_df %>% 
    filter(param %in% param_to_plot) %>% 
    ggplot(aes(x = value, color = as.factor(chain))) +
    facet_wrap(~param, scales = "free", ncol = ncol) +
    geom_density() +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = "Density plot")
}
