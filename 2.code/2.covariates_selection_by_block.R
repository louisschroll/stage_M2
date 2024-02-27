# HEADER ------------------------------------------------------------------------
#
# Script name:  model_selection_step1.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-19
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

source("2.code/format_data_for_spOccupancy.R")
source("2.code/model_selection_functions.R")

# load data
load("1.data/all_seabirds_counts.rdara")
load("1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

test_covariates_blocks <- function(data.int, static_covs, sst_covs, dyn_covs, species){
  # Static covariates selection
  static_cov_combination <- generateAllCombinations(static_covs)
  df_static_cov <- test_all_models(static_cov_combination, data.int)

  # temperature cov selection
  SST_combination <- generateAllCombinations(sst_covs) %>%
    c(list("mean_SST", c("mean_SST", "sd_SST")))
  df_SST_cov <- test_all_models(SST_combination, data.int)

  # Dynamic cov selection
  dyn_cov_combination <- generateAllCombinations(dyn_covs) 
  df_dyn_cov <- test_all_models(dyn_cov_combination, data.int)
  
  # Create a new workbook
  wb <- createWorkbook()
  addSelectionSheet(wb, sheet_name = "static_covs", df = df_static_cov, datasets_nb=length(data.int$y))
  addSelectionSheet(wb, sheet_name = "SST_covs", df = df_SST_cov, datasets_nb=length(data.int$y))
  addSelectionSheet(wb, sheet_name = "dynamic_covs", df = df_dyn_cov, datasets_nb=length(data.int$y))
  
  # Save the workbook
  file_path <- paste0("3.results/model_selection/",
                      str_replace(species, " ","_"), 
                      "_2_covariates_blocks.xlsx")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  print(paste("results in", file_path))
}

# 1. ------ Hors repro --------
data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

species_list <- c("sterne caugek", "goeland leucophee", "petit puffin", "puffin de scopoli")
                  # "mouette melanocephale", "oceanite tempete")

static_cov_list <- list(
  sterne_caugek = c("log_dist_to_shore", "concavity", "slope", "log_bathymetry"),
  goeland_leucophee = c("log_dist_to_shore", "concavity", "log_slope", "log_bathymetry"),
  petit_puffin = c("log_dist_to_shore", "log_slope", "log_bathymetry"),
  petit_puffin = c("log_dist_to_shore", "slope", "bathymetry")
)

SST_cov_list <- list(
  sterne_caugek = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  petit_puffin = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  grand_puffin = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST")
)

dyn_cov_list <- list(
  sterne_caugek = c("mean_CHL", "mean_VEL", "log_sd_VEL"),
  goeland_leucophee = c("mean_CHL", "mean_VEL", "log_sd_VEL"),
  petit_puffin = c("mean_CHL", "mean_VEL", "log_sd_VEL"),
  grand_puffin = c("mean_CHL")
)


for (i in seq_along(species_list)){
  species <- species_list[i]
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_covariates_blocks(data.int = data.int, 
                         static_covs = static_cov_list[[i]], 
                         sst_covs = SST_cov_list[[i]], 
                         dyn_covs = dyn_cov_list[[i]], 
                         species = species)
}

# 1. ------ during repro (summer) --------
# data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
#                  samm = list(obs = samm_obs, eff = samm_eff),
#                  pnm = list(obs = pnm_obs, eff = pnm_eff),
#                  migralion = list(obs = migralion_obs, eff = migralion_eff))
# 
# species_list <- c("sterne caugek", "goeland leucophee", "puffin yelkouan")
# # "mouette melanocephale", "puffin de scopoli", "oceanite tempete",
# # "mouette pygmee")
# 
# static_cov_list <- list(
#   sterne_caugek = c("log_dist_to_shore", "concavity", "slope", "bathymetry"),
#   goeland_leucophee = c("dist_to_shore", "concavity", "slope + I(slope)^2", "bathymetry"),
#   puffin_yelkouan = c("dist_to_shore", "concavity", "slope", "bathymetry")
# )
# 
# SST_cov_list <- list(
#   sterne_caugek = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
#                     "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST"),
#   goeland_leucophee = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
#                         "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST"),
#   puffin_yelkouan = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST",
#                       "sd_winter_SST", "sd_spring_SST", "sd_summer_SST", "sd_autumn_SST")
# )
# 
# dyn_cov_list <- list(
#   sterne_caugek = c("mean_CHL", "mean_VEL", "log_sd_VEL"),
#   goeland_leucophee = c("mean_CHL", "mean_VEL", "sd_VEL"),
#   puffin_yelkouan = c("mean_CHL", "mean_VEL", "sd_VEL")
# )
# 
# 
# for (i in seq_along(species_list)){
#   species <- species_list[i]
#   data.int <- get_data_for_spOccupancy(data_list, grid, species)
#   test_covariates_blocks(data.int = data.int, 
#                          static_covs = static_cov_list[[i]], 
#                          sst_covs = SST_cov_list[[i]], 
#                          dyn_covs = dyn_cov_list[[i]], 
#                          species = species)
# }
