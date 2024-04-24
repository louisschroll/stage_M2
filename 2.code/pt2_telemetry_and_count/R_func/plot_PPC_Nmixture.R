#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/R_func/plot_PPC_Nmixture.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-23
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------


plot_PPC_Nmixture <- function(samplesNmixture, datasets_names=NULL){
  
  tibble_Nmix <- retrieve_fit_values(samplesNmixture)
  # nb_dataset <- max(as.numeric(tibble_Nmix$data_source))
  # if (is.null(datasets_names)){
  #   datasets_names <- as.character(1:nb_dataset)
  # }
  # tibble_Nmix <- tibble_Nmix %>%
  #   mutate(data_source = rep(datasets_names, nrow(.) / nb_dataset))
  # 
  pvalues_df <- tibble_Nmix %>%
    group_by(data_source) %>%
    summarise(x_max = max(fit),
              x_min = min(fit),
              y_max = max(fit_rep),
              y_min = min(fit_rep)) %>%
    mutate(x = 0.1 * (x_max - x_min) + x_min,
           y = 0.9 * (y_max - y_min) + y_min) %>%
    full_join(compute_pvalue_nmix(tibble_Nmix = tibble_Nmix),
              by = join_by(data_source))
  
  plot <- ggplot(data = tibble_Nmix, aes(x = fit, y = fit_rep, color = fit<fit_rep)) +
    geom_point(shape = 20) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.75) +
    facet_wrap(~data_source, scales = "free") +
    geom_text(data = pvalues_df, mapping = aes(x = x, y = y, label = pvalue), color = "black") +
    theme_bw() +
    theme(legend.position = "none") +
    labs(title = paste0("Posterior predictive check"),
         subtitle = "Calculated for each dataset",
         x = "Observed values",
         y = "Simulated values")
  
  return(plot)
}


compute_pvalue_nmix <- function(tibble_Nmix = NULL, samples_Nmix = NULL){
  if (is.null(tibble_Nmix)){
    tibble_Nmix <- retrieve_fit_values(samplesNmixture)
  }
  pval_df <- tibble_Nmix %>% 
    group_by(data_source) %>% 
    summarise(pvalue = round(sum(fit_rep > fit) / length(fit), 2))
  
  return(pval_df)
}


retrieve_fit_values <- function(samplesNmixture){
  samplesNmixture %>% map( ~ {
    .x %>% as_tibble() %>% janitor::clean_names()
  }) %>%
    bind_rows() %>%
    select(starts_with("fit")) %>%
    mutate(id = 1:nrow(.)) %>%
    pivot_longer(-id) %>%
    mutate(type = ifelse(str_detect(string = name, pattern = "fit_rep"), "fit_rep", "fit"),
           data_source = str_sub(name,-1)) %>%
    select(-name) %>%
    pivot_wider(names_from = type, values_from = value) 
}
