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
load("1.data/all_seabirds_counts.rdata")
load("1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

test_covariates_blocks <- function(data.int, static_covs, sst_covs, dyn_covs, species){
  df_null_model <- test_all_models("1", data.int)
  # Static covariates selection
  static_cov_combination <- generateAllCombinations(static_covs)
  static_cov_combination <- static_cov_combination[2:length(static_cov_combination)]
  df_static_cov <- df_null_model %>% 
    bind_rows(test_all_models(static_cov_combination, data.int))

  # temperature cov selection
  SST_combination <- generateAllCombinations(sst_covs) %>%
    c(list("mean_SST", "sd_SST", c("mean_SST", "sd_SST")))
  SST_combination <- SST_combination[2:length(SST_combination)]
  df_SST_cov <- df_null_model %>% 
    bind_rows(test_all_models(SST_combination, data.int))

  # Dynamic cov selection
  dyn_cov_combination <- generateAllCombinations(dyn_covs) 
  dyn_cov_combination <- dyn_cov_combination[2:length(dyn_cov_combination)]
  df_dyn_cov <- df_null_model %>% 
    bind_rows(test_all_models(dyn_cov_combination, data.int))
  
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

species_list <- c("goeland_leucophee_R", "goeland_leucophee_HR")

# species_list <- migralion_obs %>%
#   filter(!is.na(species_name)) %>%
#   pull(species_name) %>%
#   unique() #%>% 
#  #str_subset("mouette", negate = F)

static_cov_list <- list(
  sterne_caugek_HR = c("log_dist_to_shore", "concavity", "log_bathymetry"),
  sterne_caugek_R = c("log_dist_to_shore", "concavity", "log_bathymetry"),
  
  goeland_leucophee_HR = c("log_dist_to_shore", "concavity", "log_bathymetry"),
  goeland_leucophee_R = c("log_dist_to_shore", "concavity", "log_bathymetry"),
  
  petit_puffin_HR = c("log_dist_to_shore", "log_bathymetry"),
  petit_puffin_R = c("dist_to_shore", "concavity", "log_bathymetry"),
  
  mouette_pygmee_HR = c("log_dist_to_shore", "log_bathymetry"),
  
  mouette_melanocephale_HR = c("dist_to_shore", "log_bathymetry"),
  mouette_melanocephale_R = c("dist_to_shore", "log_bathymetry")
)

SST_cov_list <- list(
  sterne_caugek_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  sterne_caugek_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  goeland_leucophee_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  petit_puffin_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  petit_puffin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),

  mouette_pygmee_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_HR = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  mouette_melanocephale_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST")
)

dyn_cov_list <- list(
  sterne_caugek_HR = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  sterne_caugek_R = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  goeland_leucophee_HR = c("mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  goeland_leucophee_R = c("mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  petit_puffin_HR = c("mean_CHL", "sd_SAL", "mean_SSH", "log_sd_SSH", "log_sd_VEL"),
  petit_puffin_R = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "log_sd_VEL"),

  mouette_pygmee_HR = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  mouette_melanocephale_HR = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "sd_VEL"),
  mouette_melanocephale_R = c("mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH", "sd_VEL")
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

