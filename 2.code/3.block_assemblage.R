# HEADER ------------------------------------------------------------------------
#
# Script name:  model_selection_step2.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-19
#
# Script description:
# This script aim at assembling the blocks obtained at step 1 and assessing 
# whether the model is performing better when adding the blocks. 
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space
species <- "sterne caugek"

# load packages
library(tidyverse)
library(sf)
library(spOccupancy)
library(openxlsx)

source("~/stage_M2/2.code/format_data_for_spOccupancy.R")
source("~/stage_M2/2.code/model_selection_functions.R")

# load data
load("~/stage_M2/1.data/all_seabirds_counts.rdata")
load("~/stage_M2/1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

test_blocks_assemblage <- function(data.int, best_static_covs, best_dyn_covs, best_SST_covs, species){
  models_to_test <- list(best_static_covs, 
                         best_dyn_covs,
                         best_SST_covs,
                         c(best_static_covs, best_dyn_covs), 
                         c(best_static_covs, best_SST_covs),
                         c(best_dyn_covs, best_SST_covs),
                         c(best_static_covs, best_dyn_covs, best_SST_covs))
  
  df_grouping <- test_all_models(models_to_test, data.int)
  # Create a new workbook
  wb <- createWorkbook()
  addSelectionSheet(wb, sheet_name = "blocks", df = df_grouping, datasets_nb=length(data.int$y))
  # Save the workbook
  file_path <- paste0("3.results/model_selection/", 
                      str_replace(species, " ","_"), 
                      "_3_blocks_assemblage.xlsx")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  print(paste("Results in", file_path))
}

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

species_list <- migralion_obs %>%
  filter(!is.na(species_name)) %>%
  pull(species_name) %>%
  unique() %>% 
  str_subset("sterne", negate = T) %>% 
  str_subset("goeland", negate = T)

best_static_covs <- list(
  sterne_caugek_HR = c("log_dist_to_shore", "log_bathymetry"),
  sterne_caugek_R = c("log_bathymetry"),
  
  goeland_leucophee_HR = c("log_dist_to_shore", "log_bathymetry"),
  goeland_leucophee_R = c("log_dist_to_shore", "log_bathymetry"),
  
  petit_puffin_HR = c("log_dist_to_shore", "log_bathymetry"),
  petit_puffin_R = c("dist_to_shore", "concavity"),
  
  mouette_melanocephale_HR = c("dist_to_shore"),
  mouette_melanocephale_R = c("log_bathymetry"),

  mouette_pygmee_HR = c("log_dist_to_shore", "log_bathymetry")
)

best_SST_covs <- list(
  sterne_caugek_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  goeland_leucophee_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  goeland_leucophee_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  petit_puffin_HR = c("mean_winter_SST", "mean_autumn_SST"),
  petit_puffin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_HR = c("mean_winter_SST", "mean_autumn_SST"),
  mouette_melanocephale_R = c("mean_SST"),
  
  mouette_pygmee_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST")
)

best_dyn_covs <- list(
  sterne_caugek_HR = c("mean_CHL", "mean_SSH"),
  sterne_caugek_R = c("mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  goeland_leucophee_HR = c("mean_CHL"),
  goeland_leucophee_R = c("mean_CHL"),

  petit_puffin_HR = c("mean_CHL", "sd_SAL", "mean_SSH"),
  petit_puffin_R = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  mouette_melanocephale_HR = c("mean_CHL", "sd_SAL"),
  mouette_melanocephale_R = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  mouette_pygmee_HR = c("sd_SAL", "log_sd_VEL")
  )

for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_blocks_assemblage(data.int = data.int, 
                         best_static_covs = best_static_covs[[species]], 
                         best_dyn_covs = best_dyn_covs[[species]], 
                         best_SST_covs = best_SST_covs[[species]], 
                         species = species)
}

