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

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# load packages
library(tidyverse)
library(sf)
library(spOccupancy)
library(openxlsx)
library(patchwork)

source("2.code/format_data_for_spOccupancy.R")
source("2.code/model_selection_functions.R")

# load data
load("1.data/all_seabirds_counts.rdara")
load("1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

species <- c("oceanite tempete", "oceanite ind")

data.int <- get_data_for_spOccupancy(data_list, grid, species)

predict_distribution <- function(data.int, selected_cov, add_spatial=FALSE){
  # Wrapper for intPGOcc() function of spOccupancy
  # data.int: a list containing the data with the correct format for intPGOcc()
  # selected_cov: a character vector with the covariates to include in the model
  
  data_copie <- data.int
  
  occ.formula <- writeFormula(selected_cov)
  
  nb_datasets <- length(data_copie$y)
  total_sites_nb <- nrow(data.int$occ.covs)
  det.formula <- as.list(rep("~ scale(transect_length) + session", nb_datasets)) %>%  
    map(as.formula)
  
    inits.list <- list(alpha = as.list(rep(0, nb_datasets)),
                       beta = 0, 
                       z = rep(1, total_sites_nb))
    
    prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                       alpha.normal = list(mean = as.list(rep(0, nb_datasets)), 
                                           var = as.list(rep(2.72, nb_datasets))))
    
    n.samples <- 10000
    n.burn <- 1500
    n.thin <- 5
    
    model_result <- intPGOcc(occ.formula = occ.formula,
                             det.formula = det.formula, 
                             data = data_copie,
                             inits = inits.list,
                             n.samples = n.samples, 
                             priors = prior.list, 
                             n.omp.threads = 1, 
                             verbose = FALSE, 
                             n.report = 2000, 
                             n.burn = n.burn, 
                             n.thin = n.thin, 
                             n.chains = 3)
    
  grid_pred <- grid %>% 
    as_tibble() %>% 
    select(all_of(selected_cov))
  
  X.0 <- cbind(1, grid_pred)
  out.int.pred <- predict(model_result, X.0)
  
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
  
  return(list(psi = psi, sdpsi = sdpsi))
}

selected_cov <- c("mean_winter_SST", "log_dist_to_shore", "sd_VEL")
aa <- predict_distribution(data.int, selected_cov)
aa$sdpsi
