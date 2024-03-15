# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/5.test_data_integration.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-03-11
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
library(cowplot)

source("2.code/format_data_for_spOccupancy.R")
source("2.code/model_selection_functions.R")
source("2.code/prediction_functions.R")

# load data
load("1.data/all_seabirds_counts.rdata")
load("1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))


best_covs <- list(
  sterne_caugek_HR = c("log_dist_to_shore", "log_bathymetry",
                       "mean_winter_SST", "mean_spring_SST", "mean_summer_SST",
                       "mean_CHL", "mean_SSH"),
  sterne_caugek_R = c("log_bathymetry",
                      "mean_winter_SST", "mean_spring_SST", "mean_summer_SST",
                      "mean_CHL", "mean_SSH", "sd_SSH"),
  
  goeland_leucophee_HR = c(),
  goeland_leucophee_R = c(),
  
  petit_puffin_HR = c("mean_CHL", "sd_SAL", "mean_SSH", "mean_winter_SST", "mean_autumn_SST"),
  petit_puffin_R = c("dist_to_shore", "mean_CHL", "sd_SAL", "log_sd_VEL"),
  
  mouette_melanocephale_HR = c("mean_CHL", "sd_SAL", "mean_winter_SST", "mean_autumn_SST"),
  mouette_melanocephale_R = c("log_bathymetry", "mean_SST"),
  
  mouette_pygmee_HR = c("log_dist_to_shore", "log_bathymetry", "sd_SAL", "log_sd_VEL")
)


test_data_integration <- function(data_list, grid, species, selected_cov){
  dataset_names <- get_dataset_names(data_list, species)
  datasets_combi <- generateAllCombinations(dataset_names)[-1]
  
  model_nb <- length(datasets_combi)
  model_names <- map_dfr(datasets_combi, ~{ 
    names <- .x
    rep(paste(names, collapse = "_"), length(names)) %>% as_tibble()
  }) %>% 
    pull(value) 
  
  run_the_good_model <- function(data_list, grid, species, selected_cov, datasets_names){
    data.int <- get_data_for_spOccupancy(data_list = data_list[datasets_names], 
                                         grid, 
                                         species)
    if (length(data.int$y) == 1)
      model_result <- run_non_int_model(data.int, selected_cov)
    else
      model_result <- run_int_model(data.int, selected_cov)
    return(model_result)
  }
  model_res_list <- map(datasets_combi, 
                        function(x) run_the_good_model(data_list, grid, species, selected_cov, x))
  
  # comparison_df <- unique(model_names) %>%
  #   bind_cols(map_dfr(model_res_list, ~{
  #                model_result <- .x
  #                get_model_stat(model_result)
  #  }))
  
  waic_df <- map_dfr(model_res_list, ~{
    waicOcc(.x) %>% round(1)
  }) %>% 
    select(WAIC) %>% 
    mutate(model = model_names,
           dataset = unique(model_names) %>% str_split(pattern = '_') %>% unlist()) %>% 
    arrange(dataset)
  
  beta_df <- map_dfr(model_res_list, ~{
    get_beta_values(model_result = .x, selected_cov, model_nb=1)
  }) %>% 
    mutate(model = rep(unique(model_names), each = length(selected_cov))) %>% 
    pivot_wider(names_from = covar, values_from = c(beta, sd_beta), names_sep = "_") 
  
  predictive_maps_list <- map(model_res_list, 
                  function(x) make_predictive_map(x, grid = grid, selected_cov = selected_cov))
  names(predictive_maps_list) <- unique(model_names)
  
  return(list(waic_df = waic_df, beta_df = beta_df, maps = predictive_maps_list))
}


save_stat <- function(integration_res, species){
  # Create a new workbook
  workbook <- createWorkbook()
  addWorksheet(workbook, "beta_values")
  # Add sheets
  writeData(workbook, sheet = "beta_values", x = integration_res$beta_df)
  addWorksheet(workbook, "WAIC")
  writeData(workbook, sheet = "WAIC", x = integration_res$waic_df)
  # Save the workbook
  file_path <- paste0("3.results/model_selection/", 
                      species, 
                      "_5_test_integration.xlsx")
  saveWorkbook(workbook, file_path, overwrite = TRUE)
  print(paste("Results in", file_path))
}


species_list <- migralion_obs %>%
  filter(!is.na(species_name)) %>%
  pull(species_name) %>%
  unique() %>% 
  str_subset(pattern = "goeland", negate = T)

for (species in species_list){
  print(species)
  integration_res <- test_data_integration(data_list, grid, species, selected_cov = best_covs[[species]])
  save_stat(integration_res, species)
  save_maps_as_pdf(integration_res$map, 
                  paste0("3.results/model_selection/maps_integration_test/", 
                         species, "_5_test_int_plot.pdf"))
}

