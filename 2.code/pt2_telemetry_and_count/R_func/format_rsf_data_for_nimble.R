#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  2.code/pt2_telemetry_and_count/R_func/format_rsf_data_for_nimble.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-10
#'
#' Script description:
#' This function creates 3 lists: inits, constants and data in the good format for 
#' the RSF model written in nimble.
#' @df_RSF: a data frame containing the RSF data already prepared, this is the 
#' output of the script prepare_data_for_RSF.R
#' @selected_cov: a character vector with the covariates to use in a model
#' -------------------------------------------------------------------------------


format_rsf_data_for_nimble <- function(df_RSF, selected_cov){
  # Add weight (this is done to have more pseudo-absence while saving computational ressources)
  weight <- df_RSF$case
  weight[weight==0] <- 1000
  
  # Define constants
  n.occ.cov <- length(selected_cov) + 1
  nindividual = n_distinct(df_RSF$individual_id)
  XN <- df_RSF %>% 
    select(all_of(selected_cov)) %>% 
    mutate(intercept = 1) %>% 
    relocate(intercept) %>% 
    as.matrix()
  
  constants <-  list(npoints = nrow(df_RSF),
                     idind = get_new_id(df_RSF, id_column = "individual_id"),
                     XN = XN,
                     nindividual = nindividual,
                     n.occ.cov = n.occ.cov,
                     w = weight)
  
  # Define data
  data <- list(kase = df_RSF$case)
  
  # Inits
  inits <-  list(beta_pop = rep(0, n.occ.cov),
                 sd_pop = runif(n.occ.cov, 0, 5),
                 beta_ind = matrix(0, ncol = n.occ.cov, nrow = nindividual)
  )
  return(list(data = data, constants = constants, inits = inits))
}
