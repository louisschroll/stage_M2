

# Load packages
library(tidyverse)
library(sf)
library(spOccupancy)

# cluster adress
adress <- "/lustre/schrolll/"
# Local adress
# adress <- ""

# Load functions
source(paste0(adress, "2.code/pt1_spOccupancy/format_data_for_spOccupancy.R"))
source(paste0(adress, "2.code/pt1_spOccupancy/model_selection_functions.R"))

# Load data
load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/grid.rdata"))

species_list <- c("labbe", "macareux_moine_HR")
# Write data in a list
data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

cov_list <- c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", 
                      "mean_autumn_SST", "mean_SST", "sd_SST", "concavity", 
                      "dist_to_shore", "bathymetry",
                      "mean_CHL", "mean_SSH", "sd_SSH", "sd_VEL")

cov_combination_list <- generateAllCombinations(cov_list)

# Test quadratic and log effect for all species
for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  result_df <- test_all_models(cov_combination_list, data.int)
  save(result_df, file = paste0(adress, "3.results/test_all_models_", species, ".rdata"))
}