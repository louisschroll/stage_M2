#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt1_spOccupancy/save_occupancy_results.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-12
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
# adress <- ""

# load packages
library(tidyverse)
library(sf)
library(spOccupancy)

# load data
load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/covariates_data.rdata"))

source(paste0(adress, "2.code/pt1_spOccupancy/prediction_functions.R"))
source(paste0(adress, "2.code/pt1_spOccupancy/format_data_for_spOccupancy.R"))
source(paste0(adress, "2.code/pt1_spOccupancy/model_selection_functions.R"))

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

spOccupancy_res_grid <- grid %>% select(grid_c)

data_list = list(
  pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
  samm = list(obs = samm_obs, eff = samm_eff),
  pnm = list(obs = pnm_obs, eff = pnm_eff),
  migralion = list(obs = migralion_obs, eff = migralion_eff)
)

best_covs <- list(
  fou_de_bassan_winter = c("dist_to_shore", "log_bathymetry"),
  
  goeland_leucophee_winter = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", 
                           "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_summer = c("log_dist_to_shore", "mean_winter_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_winter = c("mean_CHL", "sd_SAL", "mean_winter_SST", "mean_autumn_SST"),
  mouette_melanocephale_summer = c("mean_autumn_SST"),
  
  mouette_pygmee_winter = c("log_dist_to_shore", "log_bathymetry", "mean_CHL", "sd_SAL", 
                        "mean_SSH", "log_sd_VEL", "mean_autumn_SST", 
                        "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  mouette_rieuse_winter = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", "mean_autumn_SST"),
  
  oceanite_tempete_summer = c("dist_to_shore", "log_bathymetry", "sd_SAL", "sd_SSH", "log_sd_VEL"),
  
  petit_puffin_winter = c("log_dist_to_shore", "log_bathymetry", 'mean_CHL', "sd_SAL", "mean_SSH"),
  petit_puffin_summer = c("mean_CHL", "mean_SSH"),
  
  puffin_de_scopoli_summer = c("log_mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  sterne_caugek_winter = c("mean_CHL", "mean_SSH", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_summer = c("log_bathymetry", "mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  sterne_pierregarin_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  macareux_moine_winter = c("log_dist_to_shore", "mean_SSH"),
  
  labbe_summer = c("mean_SSH", "sd_SSH", "sd_SAL", "mean_autumn_SST", "mean_winter_SST")
)

species_list <- migralion_obs %>%
  filter(!is.na(species_name)) %>%
  pull(species_name) %>%
  unique()

n_iter <- length(species_list)

for (i in 1:n_iter){
  
  species <- species_list[i]
  print(species)
  selected_cov <- best_covs[[species]]
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  model_result <- run_model_without_kfold(data.int = data.int, 
                                          grid=grid, 
                                          species=species, 
                                          selected_cov=selected_cov,
                                          n.samples = 50000)
  save(model_result, file = paste0(adress, "3.results/spOccupancy_outputs/spOccupancy_",species,".RData"))
  
  grid2 <- put_results_in_grid(grid, model_result = model_result, selected_cov) %>% 
    select(mean.psi, sd.psi)
  
  spOccupancy_res_grid[[paste0("mean_psi_", species)]] <- grid2$mean.psi
  spOccupancy_res_grid[[paste0("sd_psi_", species)]] <- grid2$sd.psi
  }

#plot(spOccupancy_res_grid)

save(spOccupancy_res_grid, file = paste0(adress, "3.results/grid_spOccupancy_results.RData"))
#st_write(spOccupancy_res_grid, dsn = "3.results/grid_spOccupancy_results.shp")
