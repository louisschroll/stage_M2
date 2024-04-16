#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-09
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# cluster
adress <- "/lustre/schrolll/"
# local 
adress <- ""

load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/grid_cells.rdata"))

# Load package
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(nimble, warn.conflicts = FALSE)

path_to_Rfunc <- paste0(adress, "2.code/pt2_telemetry_and_count/R_func")
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

# Prepare data
species_vector <- c("sterne_caugek_R", "sterne_caugek_HR",
                    "goeland_leucophee_HR", "goeland_leucophee_R")

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  #pnm = list(obs_data = pnm_obs, effort_data = pnm_eff),
                  #samm = list(obs_data = samm_obs, effort_data = samm_eff)
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


n.iter = 1500
n.burnin = 0.1 * n.iter
n.chains = 2

for (i in seq_along(species_vector)){
  species <- species_vector[i]
  covar <- best_cov[[species]]
  data_nmix <- prepare_data_Nmix(data_list, grid, species, selected_cov=covar)
  samplesNmixture <- run_Nmixture(data_nmix = data_nmix, 
                                  n.iter = n.iter, 
                                  n.burnin = n.burnin, 
                                  n.chains = n.chains, 
                                  compute_pvalues = TRUE)
}

sss <- map(samplesNmixture, ~{as_tibble(.x)}) %>% 
  bind_rows() %>% 
  janitor::clean_names()

MCMCvis::MCMCsummary(samplesNmixture)

plot(sss$fit_2, sss$fit_rep_2)
abline(a = 0, b = 1, col = "red")
