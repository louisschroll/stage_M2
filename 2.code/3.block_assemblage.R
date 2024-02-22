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
load("~/stage_M2/1.data/all_seabirds_counts.rdara")
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
  
  df_grouping <- test_all_models(models_to_test, data.int, det.formula)
  # Create a new workbook
  wb <- createWorkbook()
  addSelectionSheet(wb, sheet_name = "blocks", df = df_grouping)
  # Save the workbook
  file_path <- paste0(str_replace(species, " ","_"), "_3_blocks_assemblage.xlsx")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  print(paste("Results in", file_path))
}

# ----- Hors repro -----
data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs2, eff = migralion_eff2))

species_list <- c("sterne caugek", "mouette pygmee", "goeland leucophee", "puffin yelkouan")
# "mouette melanocephale", "puffin de scopoli", "oceanite tempete",
# "mouette pygmee")

best_static_cov_list <- list(
  sterne_caugek = c("log_dist_to_shore", "concavity", "slope", "bathymetry"),
  goeland_leucophee = c("dist_to_shore", "concavity", "slope + I(slope)^2", "bathymetry"),
  puffin_yelkouan = c("dist_to_shore", "concavity", "slope", "bathymetry")
)

best_SST_cov_list <- list(
  sterne_caugek = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
                    "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST"),
  goeland_leucophee = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
                        "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST"),
  puffin_yelkouan = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
                      "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST")
)

best_dyn_cov_list <- list(
  sterne_caugek = c("mean_CHL", "mean_VEL", "log_sd_VEL"),
  goeland_leucophee = c("mean_CHL", "mean_VEL", "sd_VEL"),
  puffin_yelkouan = c("mean_CHL", "mean_VEL", "sd_VEL")
)

for (i in seq_along(species_list)){
  species <- species_list[i]
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_blocks_assemblage(data.int, best_static_covs, best_dyn_covs, best_SST_covs, species)
}

# ----- repro ------
