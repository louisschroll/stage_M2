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

source("2.code/pt1_spOccupancy/format_data_for_spOccupancy.R")
source("2.code/pt1_spOccupancy/model_selection_functions.R")

# load data
load("1.data/all_seabirds_counts.rdata")
load("1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

test_covariates_blocks <- function(data.int, static_covs, sst_covs, dyn_covs, species){
  df_null_model <- test_all_models("1", data.int)$comparison_df
  # Static covariates selection
  static_cov_combination <- generateAllCombinations(static_covs)
  static_cov_combination <- static_cov_combination[2:length(static_cov_combination)]
  df_static_cov <- df_null_model %>% 
    bind_rows(test_all_models(static_cov_combination, data.int)$comparison_df)

  # temperature cov selection
  SST_combination <- generateAllCombinations(sst_covs) %>%
    c(list("mean_SST", "sd_SST", c("mean_SST", "sd_SST")))
  SST_combination <- SST_combination[2:length(SST_combination)]
  df_SST_cov <- df_null_model %>% 
    bind_rows(test_all_models(SST_combination, data.int)$comparison_df)

  # Dynamic cov selection
  dyn_cov_combination <- generateAllCombinations(dyn_covs) 
  dyn_cov_combination <- dyn_cov_combination[2:length(dyn_cov_combination)]
  df_dyn_cov <- df_null_model %>% 
    bind_rows(test_all_models(dyn_cov_combination, data.int)$comparison_df)
  
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

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

species_list <- migralion_obs %>%
  filter(!is.na(species_name)) %>%
  pull(species_name) %>%
  unique() 

species_list <- c(species_list[8:length(species_list)], "goeland_leucophee_summer")

SST_cov_list <- list(
  sterne_caugek_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  sterne_caugek_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  goeland_leucophee_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  petit_puffin_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  petit_puffin_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_pygmee_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  mouette_melanocephale_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  fou_de_bassan_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  sterne_pierregarin_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  oceanite_tempete_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_rieuse_winter = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  puffin_de_scopoli_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  labbe_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  macareux_moine_summer = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST")
)

static_cov_list <- list(
  sterne_caugek_winter = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  sterne_caugek_summer = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  goeland_leucophee_winter = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  goeland_leucophee_summer = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  petit_puffin_winter = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "log_sd_SSH"),
  petit_puffin_summer = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  mouette_pygmee_winter = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  mouette_melanocephale_winter = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  mouette_melanocephale_summer = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  fou_de_bassan_winter = c("dist_to_shore", "log_bathymetry", "mean_SSH", "log_sd_SSH"),
  
  sterne_pierregarin_summer = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  oceanite_tempete_summer = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  mouette_rieuse_winter = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  puffin_de_scopoli_summer = c("dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH"),
  
  labbe_summer = c("dist_to_shore", "bathymetry", "mean_SSH", "sd_SSH"),
  
  macareux_moine_summer = c("log_dist_to_shore", "log_bathymetry", "mean_SSH", "sd_SSH")
)


dyn_cov_list <- list(
  sterne_caugek_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  sterne_caugek_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  goeland_leucophee_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  goeland_leucophee_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  petit_puffin_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  petit_puffin_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL"),

  mouette_pygmee_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  mouette_melanocephale_winter = c("mean_CHL", "sd_SAL", "sd_VEL"),
  mouette_melanocephale_summer = c("log_mean_CHL", "sd_SAL", "sd_VEL"),
  
  fou_de_bassan_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  sterne_pierregarin_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  oceanite_tempete_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  mouette_rieuse_winter = c("mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  puffin_de_scopoli_summer = c("log_mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  labbe_summer = c("mean_CHL", "sd_SAL", "sd_VEL"),
  
  macareux_moine_summer = c("mean_CHL", "sd_SAL", "log_sd_VEL")
)


for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_covariates_blocks(data.int = data.int, 
                         static_covs = static_cov_list[[species]], 
                         sst_covs = SST_cov_list[[species]], 
                         dyn_covs = dyn_cov_list[[species]], 
                         species = species)
}

