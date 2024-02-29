# HEADER ------------------------------------------------------------------------
#
# Script name:  
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-27
#
# Script description:
#
#
# -------------------------------------------------------------------------------


source("2.code/format_data_for_spOccupancy.R")
source("2.code/model_selection_functions.R")

predict_distribution <- function(data.int = NULL, data_list = NULL, grid, species, selected_cov, add_spatial=FALSE){
  # Wrapper for intPGOcc() function of spOccupancy
  # data_list:
  # selected_cov: a character vector with the covariates to include in the model
  if (is.null(data.int))
    data.int <- get_data_for_spOccupancy(data_list, grid, species)

  occ.formula <- writeFormula(selected_cov)
  
  nb_datasets <- length(data.int$y)
  total_sites_nb <- nrow(data.int$occ.covs)
  det.formula <- as.list(rep("~ scale(transect_length) + session", nb_datasets)) %>%  
    map(as.formula)
    
    n.samples <- 10000
    n.burn <- 1500
    n.thin <- 5
    n.chains <- 3
    if (nb_datasets>1){
      inits.list <- list(alpha = as.list(rep(0, nb_datasets)),
                         beta = 0, 
                         z = rep(1, total_sites_nb))
      
      prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                         alpha.normal = list(mean = as.list(rep(0, nb_datasets)), 
                                             var = as.list(rep(2.72, nb_datasets))))
      
      model_result <- intPGOcc(occ.formula = occ.formula,
                               det.formula = det.formula, 
                               data = data.int,
                               inits = inits.list,
                               n.samples = n.samples, 
                               priors = prior.list, 
                               verbose = FALSE, 
                               n.burn = n.burn, 
                               n.thin = n.thin, 
                               n.chains = n.chains)
    }
    
  if (nb_datasets==1){
    data <- list(y = data.int$y[[1]], 
                 occ.covs = data.int$occ.covs,
                 det.covs = data.int$det.covs[[1]])
    
    inits.list <- list(alpha = 0, 
                       beta = 0, 
                       z = apply(data$y, 1, max, na.rm = TRUE))
    
    prior.list <- list(alpha.normal = list(mean = 0, var = 2.72), 
                        beta.normal = list(mean = 0, var = 2.72))
    
    model_result <- PGOcc(occ.formula = occ.formula, 
                 det.formula = det.formula[[1]], 
                 data = data, 
                 inits = inits.list, 
                 n.samples = n.samples, 
                 priors = prior.list, 
                 verbose = F, 
                 n.burn = n.burn, 
                 n.thin = n.thin, 
                 n.chains = n.chains)
  }
  grid_pred <- grid %>% 
    as_tibble() %>% 
    select(all_of(selected_cov))
  
  X.0 <- cbind(1, grid_pred)
  out.int.pred <- predict(model_result, as.matrix(X.0))
  
  mean.psi = apply(out.int.pred$psi.0.samples, 2, mean)
  sd.psi = apply(out.int.pred$psi.0.samples, 2, sd)
  
  psi <- ggplot() + 
    geom_sf(data = grid, aes(fill = mean.psi), lwd = 0.1) +
    scale_fill_viridis_c() + 
    labs(title = 'Occupancy') +
    theme_bw()
  
  sdpsi <- ggplot() + 
    geom_sf(data = grid, aes(fill = sd.psi), lwd = 0.1) +
    scale_fill_viridis_c(option = "B") + 
    labs(title = 'Occupancy SD') +
    theme_bw()
  
  return(list(res = model_result, psi = psi, sdpsi = sdpsi))
}


