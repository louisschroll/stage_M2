#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-09
#'
#' Script description:
#' Run N-mixture for the selected species and save the mcmc outputs (in 3.results
#' /mcmc_outputs/) as well as the prediction (in 3.results.prediction_grid/)
#'
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# cluster
adress <- "/lustre/schrolll/"
# local 
# adress <- ""

# Load the data
load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/grid_cells.rdata"))

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
                  samm = list(obs_data = samm_obs, effort_data = samm_eff)
)

# Load packages
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(nimble, warn.conflicts = FALSE)

# Load functions
path_to_Rfunc <- paste0(adress, "2.code/pt2_telemetry_and_count/R_func")
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

# Select the species
species_vector <- c("sterne_pierregarin_R",
                    "mouette_rieuse_HR")

# Indicate the covar to use
best_cov <- list(
  sterne_caugek_winter = c("mean_CHL", "mean_SSH", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_summer = c("log_bathymetry", "mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  sterne_pierregarin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  mouette_rieuse_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", "mean_autumn_SST")
)

# Run the model and save the results
n.iter = 100000
n.burnin = 0.1 * n.iter
n.chains = 3

for (i in seq_along(species_vector)){
  species <- species_vector[i]
  print(species)
  covar <- best_cov[[species]]
  data_nmix <- prepare_data_Nmix(data_list = data_list, 
                                 grid = grid, 
                                 species = species, 
                                 selected_cov = covar)
  samplesNmixture <- run_Nmixture(data_nmix = data_nmix, 
                                  n.iter = n.iter, 
                                  n.burnin = n.burnin, 
                                  n.chains = n.chains, 
                                  compute_pvalues = T)
  ## save mcmc output in a .rdata
  save(samplesNmixture, file = paste0(adress, "3.results/mcmc_outputs/Nmix_output_", species, ".rdata"))
  
  ## predict and save prediction as a sf object in a .rdata
  grid_nmix <- make_prediction(samplesNmixture, grid, selected_cov=covar)
  save(grid_nmix, file = paste0(adress, "3.results/prediction_grid/grid_Nmix_", species, ".rdata"))
}



