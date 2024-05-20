# load package
library(tidyverse)
library(sf)

local_path <- "~/stage_M2/"

# load data
load(paste0(local_path, "1.data/grid_cells.rdata"))

grid_nmix_all_sp <- grid %>% select(grid_c)

best_covs <- list(
  fou_de_bassan_HR = c("dist_to_shore", "log_bathymetry"),
  
  goeland_leucophee_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", 
                               "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_R = c("log_dist_to_shore", "mean_winter_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_HR = c("mean_CHL", "sd_SAL", "mean_winter_SST", "mean_autumn_SST"),
  mouette_melanocephale_R = c("mean_autumn_SST"),
  
  mouette_pygmee_HR = c("log_dist_to_shore", "log_bathymetry", "mean_CHL", "sd_SAL", 
                            "mean_SSH", "log_sd_VEL", "mean_autumn_SST", 
                            "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  mouette_rieuse_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", "mean_autumn_SST"),
  
  oceanite_tempete = c("dist_to_shore", "log_bathymetry", "sd_SAL", "sd_SSH", "log_sd_VEL"),
  
  petit_puffin_HR = c("log_dist_to_shore", "log_bathymetry", 'mean_CHL', "sd_SAL", "mean_SSH"),
  petit_puffin_R = c("mean_CHL", "mean_SSH"),
  
  puffin_de_scopoli_R = c("log_mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  sterne_caugek_HR = c("mean_CHL", "mean_SSH", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_R = c("log_bathymetry", "mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  sterne_pierregarin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  macareux_moine_HR = c("log_dist_to_shore", "mean_SSH"),
  
  labbe = c("mean_SSH", "sd_SSH", "sd_SAL", "mean_autumn_SST", "mean_winter_SST")
)

species_nmix <- c(
  "fou_de_bassan_HR",
  "mouette_melanocephale_HR",
  "mouette_melanocephale_R",
  "mouette_pygmee_HR",
  "mouette_rieuse_HR",
  "oceanite_tempete",
  "sterne_pierregarin_R",
  #"macareux_moine_HR",
  "labbe"
)

path_nmix <- paste0(local_path, 
                    "3.results/mcmc_outputs/Nmix_output_", 
                    species_nmix, 
                    ".rdata")

for (i in seq_along(species_nmix)){
  species <- species_nmix[i]
  print(species)
  load(path_nmix[i])
  grid_nmix <- make_prediction(mcmc.output = samplesNmixture, 
                           grid = grid, 
                           selected_cov = best_covs[[species]]) %>% 
    select(mean_psi, sd_psi) %>% 
    mutate(mean_psi = mean_psi / max(mean_psi))
  
  grid_nmix_all_sp[[paste0("mean_psi_", species)]] <- grid_nmix$mean_psi
  grid_nmix_all_sp[[paste0("sd_psi_", species)]] <- grid_nmix$sd_psi
}


species_integre <- c(
  "goeland_leucophee_HR",
  "goeland_leucophee_R",
  "petit_puffin_HR",
  "petit_puffin_R",
  "puffin_de_scopoli_R",
  "sterne_caugek_HR",
  "sterne_caugek_R"
)

path_modelint <- paste0(local_path, 
                        "3.results/mcmc_outputs/mcmc_output_", 
                        species_integre, 
                        ".rdata")


for (i in seq_along(species_integre)){
  species <- species_integre[i]
  load(path_modelint[i])
  grid_int <- make_prediction(samplesint, 
                              grid, 
                              selected_cov = best_covs[[species]], 
                              rsf_intercept = "beta_pop[1]") %>% 
    select(mean_psi, sd_psi) %>% 
    mutate(mean_psi = mean_psi / max(mean_psi))
  
  grid_nmix_all_sp[[paste0("mean_psi_", species)]] <- grid_int$mean_psi
  grid_nmix_all_sp[[paste0("sd_psi_", species)]] <- grid_int$sd_psi
}

save(grid_nmix_all_sp, file = paste0(adress, "3.results/prediction_grid/grid_nmix_all_sp.RData"))
