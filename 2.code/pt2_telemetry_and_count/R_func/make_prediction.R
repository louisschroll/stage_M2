#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  make_prediction.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-03
#'
#' Script description:
#' This function predicts space use intensity by birds using the output from the
#' N-mixture or the RSF. 
#' @mcmc.ouput: a mcmc.list, output from a nimble model (RSF, N-mixture in our case),
#' make sure that occurence coefficient names begin with 'beta'
#' @grid: a sf object, the grid of the study area containing the covariate values
#' @selected_cov: a character vector containing the names of the covariates used 
#' in the model
#' @output: a sf object with the same format than grid, containg the values of 
#' the selected covariates and the mean and sd values of the predicted space use 
#' intensity. 
#' -------------------------------------------------------------------------------


make_prediction <- function(mcmc.output, grid, selected_cov, include_intercept=TRUE, rsf_intercept=NULL){
  new_grid <- grid %>% select(all_of(selected_cov))
  
  # Get coefficients values from mcmc output
  if (class(mcmc.output)=="mcmc"){
    coeff_values <- mcmc.output %>% as_tibble() %>%  select(starts_with("beta")) %>% as.matrix()
  } else {
    coeff_values <- map(mcmc.output, ~{.x %>% as_tibble() %>% select(starts_with("beta"))}) %>% 
      bind_rows() %>% 
      select(-all_of(c(rsf_intercept))) %>% 
      as.matrix()
  }
  
  # Get covariates values within each cells of grid
  if (include_intercept){
    occurence_covs <- new_grid %>% 
      st_drop_geometry() %>% 
      mutate(intercept = 1) %>% 
      relocate(intercept) %>% 
      as.matrix()
  } else {
    occurence_covs <- new_grid %>% 
      st_drop_geometry() %>% 
      as.matrix()
  }
  
  # Multiply each rows of occurence_covs matrix by each rows of coeff_values
  psi_pred <- apply(occurence_covs, 1, function(occurence_row){
    coeff_values %*% occurence_row
  })
  
  # Store results
  new_grid$mean_psi <- apply(psi_pred, 2, mean)
  new_grid$sd_psi <- apply(psi_pred, 2, sd)
  return(new_grid)
}

