# HEADER ------------------------------------------------------------------------
#
# Script name:  model_selection_step3.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-19
#
# Script description:
# Assess and compare performance of the model choosed after step 1 and 2 with
# and without adding spatial auto-correlation.  
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

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

test_adding_sp_effect <- function(data.int, species, best_covs){
  selected_cov = best_covs[[species]]
  data_copie <- data.int
  data_copie$occ.covs <- data_copie$occ.covs %>% select(selected_cov)

  without_spatial <- run_int_model(data_copie, selected_cov=selected_cov, add_spatial=FALSE)
  with_spatial_exp <- run_int_model(data_copie, selected_cov=selected_cov, add_spatial=TRUE, spatial_model="exponential")
  with_spatial_sphere <- run_int_model(data_copie, selected_cov=selected_cov, add_spatial=TRUE, spatial_model="spherical")
  with_spatial_gaussian <- run_int_model(data_copie, selected_cov=selected_cov, add_spatial=TRUE, spatial_model="gaussian")

  stat_df <- bind_rows(get_model_stat(without_spatial),
                      get_model_stat(with_spatial_exp),
                     get_model_stat(with_spatial_sphere),
                     get_model_stat(with_spatial_gaussian))

  beta_df <- bind_rows(get_beta_values(without_spatial, selected_cov, model_nb=1),
                     get_beta_values(with_spatial_exp, selected_cov, model_nb=2),
                     get_beta_values(with_spatial_sphere, selected_cov, model_nb=3),
                     get_beta_values(with_spatial_gaussian, selected_cov, model_nb=4))%>% 
    pivot_wider(names_from = covar, values_from = c(beta, sd_beta), names_sep = "_") 

  comparison_df <- bind_cols(c("without spatial", "with spatial exp", 
                             "with sp shpere", "with sp gaussian"), 
                           stat_df, beta_df) 

  # Create a new workbook
  wb <- createWorkbook()
  addSelectionSheet(wb, sheet_name = "spatial_effect", df = comparison_df, datasets_nb=length(data.int$y))
  # Save the workbook
  file_path <- paste0("3.results/model_selection/", 
                    str_replace(species, " ","_"), 
                    "_4_add_spatial_effects.xlsx")
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

best_covs <- list(
  fou_de_bassan_HR = c("dist_to_shore", "log_bathymetry"),
  
  goeland_leucophee_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", 
                           "mean_spring_SST", "mean_summer_SST", "mean_autumn_SST"),
  goeland_leucophee_R = c("log_dist_to_shore", "mean_winter_SST", "mean_summer_SST", "mean_autumn_SST"),
  
  mouette_melanocephale_HR = c("mean_CHL", "sd_SAL", "mean_winter_SST", "mean_autumn_SST"),
  mouette_melanocephale_R = c("mean_autumn_SST"),
  
  mouette_pygmee_HR = c("log_dist_to_shore", "log_bathymetry", "mean_CHL", "sd_SAL", "mean_SSH", "log_sd_VEL", "mean_autumn_SST", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  mouette_rieuse_HR = c("log_dist_to_shore", "log_bathymetry", "mean_winter_SST", "mean_autumn_SST"),
  
  oceanite_tempete = c("dist_to_shore", "log_bathymetry", "sd_SAL", "sd_SSH", "log_sd_VEL"),
  
  petit_puffin_HR = c("log_dist_to_shore", "log_bathymetry", 'mean_CHL', "sd_SAL", "mean_SSH"),
  petit_puffin_R = c("mean_CHL", "mean_SSH"),
  
  puffin_de_scopoli_R = c("log_mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  sterne_caugek_HR = c("mean_CHL", "mean_SSH", "mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  sterne_caugek_R = c("log_bathymetry", "mean_CHL", "mean_SSH", "sd_SSH", "log_sd_VEL"),
  
  sterne_pierregarin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST")
)


for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species, add_coords = TRUE)
  test_adding_sp_effect(data.int, species, best_covs)
}
