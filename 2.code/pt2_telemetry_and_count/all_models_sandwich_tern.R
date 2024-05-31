#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-09
#'
#' Script description:
#' Run N-mixture, RSF, and the integrated RSF-Nmixture for the sandwich tern, 
#' and save the mcmc outputs and the predicted grid
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# cluster
adress <- "/lustre/schrolll/"
# local 
# adress <- ""

load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/grid_cells.rdata"))

# Load package
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(nimble, warn.conflicts = FALSE)

path_to_Rfunc <- paste0(adress, "2.code/pt2_telemetry_and_count/R_func")
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

# Prepare data
species_vector <- c("sterne_caugek_R", "sterne_caugek_HR")

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
                  samm = list(obs_data = samm_obs, effort_data = samm_eff)
)

best_cov <- list(
  goeland_leucophee_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", 
                           "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_R = c("log_dist_to_shore", "mean_winter_SST", "mean_summer_SST", "mean_autumn_SST"),

  petit_puffin_HR = c("log_dist_to_shore", "log_bathymetry", 'mean_CHL', "sd_SAL", "mean_SSH"),
  petit_puffin_R = c("mean_CHL", "mean_SSH"),
  
  puffin_de_scopoli_R = c("log_mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  sterne_caugek_HR = c("mean_CHL", "mean_SSH", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_R = c("log_bathymetry", "mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL")
)

month_to_keep <- list(
  goeland_leucophee_HR = c(1:3, 9:12),
  goeland_leucophee_R = 4:8,
  
  petit_puffin_HR = c(1:3, 9:12),
  petit_puffin_R = 4:8,
  
  puffin_de_scopoli_R = 1:12,
  
  sterne_caugek_HR = c(1:3, 9:12),
  sterne_caugek_R = 4:8
)

file_rsf_data <- list(
  goeland_leucophee_HR = "1.data/RSF_data_yellow_legged_gull.rdata",
  goeland_leucophee_R = "1.data/RSF_data_yellow_legged_gull.rdata",
  
  petit_puffin_HR = "1.data/RSF_puffin_yelkouan.rdata",
  petit_puffin_R = "1.data/RSF_puffin_yelkouan.rdata",
  
  puffin_de_scopoli_R = "1.data/RSF_data_scopoli_shearwaters.rdata",
  
  sterne_caugek_HR = "1.data/RSF_data_sandwich_tern.rdata",
  sterne_caugek_R = "1.data/RSF_data_sandwich_tern.rdata"
)


n.iter = 50000
n.burnin = 0.1 * n.iter
n.chains = 3

for (i in seq_along(species_vector)){
  species <- species_vector[i]
  covar <- best_cov[[species]]
  # 1/3 - N-mixture 
  data_nmix <- prepare_data_Nmix(data_list = data_list,
                                 grid = grid,
                                 species = species,
                                 selected_cov = covar)
  # samplesNmixture <- run_Nmixture(data_nmix = data_nmix,
  #                                 n.iter = n.iter,
  #                                 n.burnin = n.burnin,
  #                                 n.chains = n.chains)
  # grid_nmix <- make_prediction(mcmc.output = samplesNmixture, 
  #                              grid = grid, 
  #                              selected_cov = covar)

  # 2/3 - RSF 
  load(paste0(adress, file_rsf_data[[species]]))
  data_rsf <- format_rsf_data_for_nimble(df_RSF %>% filter(month %in% month_to_keep[[species]]), 
                                         covar)
  # samplesRSF <- run_RSF(data_rsf, 
  #                       n.iter = n.iter, 
  #                       n.burnin=n.burnin, 
  #                       n.chains=n.chains)
  # grid_RSF <- make_prediction(mcmc.output = samplesRSF, 
  #                             grid, 
  #                             selected_cov = covar, 
  #                             include_intercept = F, 
  #                             rsf_intercept = "beta_pop[1]")

  # 3/3 - Integrated RSF and N-mixture
  samplesint <- run_integrated_Nmix_RSF(data_nmix = data_nmix, 
                                        data_rsf = data_rsf, 
                                        nmix_model = "NB",
                                        n.iter = n.iter, 
                                        n.burnin = n.burnin,
                                        n.chains = n.chains)
  save(samplesint, 
       file = paste0(adress, "3.results/mcmc_outputs/int_output_", species, ".rdata"))
  
  grid_int <- make_prediction(mcmc.output = samplesint, 
                              grid = grid, 
                              selected_cov = covar, 
                              include_intercept = F, 
                              rsf_intercept = c("beta_pop[1]", "beta0_nmix"))
  save(grid_int,
       file = paste0(adress, "3.results/prediction_grid/grid_int_", species, ".rdata"))
  
  # Save everything
  # save(samplesNmixture, samplesRSF, samplesint, 
  #      file = paste0(adress, "3.results/mcmc_outputs/mcmc_output_", species, ".rdata"))
  # save(grid_nmix, grid_RSF, grid_int,
  #      file = paste0(adress, "3.results/prediction_grid/grid_", species, ".rdata"))
}

