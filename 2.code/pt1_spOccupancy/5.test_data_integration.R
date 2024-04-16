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

source("2.code/pt1_spOccupancy/format_data_for_spOccupancy.R")
source("2.code/pt1_spOccupancy/model_selection_functions.R")
source("2.code/pt1_spOccupancy/prediction_functions.R")

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
  
  sterne_pierregarin_R = c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST"),
  
  labbe = c("mean_SSH", "sd_SSH", "sd_SAL", "mean_autumn_SST", "mean_winter_SST"),
  
  macareux_moine_HR = c("log_dist_to_shore", "mean_SSH", "mean_autumn_SST", "I(mean_autumn_SST)^2", "mean_CHL")
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
  stat_df <- map_dfr(model_res_list, ~{
    model_result <- .x
    tibble(
      pavlue_CS1 = ppcOcc(model_result, fit.stat = "chi-squared", group = 1) %>% compute_bayesian_pvalue(),
      pavlue_CS2 = ppcOcc(model_result, fit.stat = "chi-squared", group = 2) %>% compute_bayesian_pvalue(),
      pavlue_FT1 = ppcOcc(model_result, fit.stat = "freeman-tukey", group = 1) %>% compute_bayesian_pvalue(),
      pavlue_FT2 = ppcOcc(model_result, fit.stat = "freeman-tukey", group = 2) %>% compute_bayesian_pvalue(),
      WAIC = waicOcc(model_result) %>% bind_rows() %>% pull(WAIC) %>% round(0),
      CV_deviance = model_result$k.fold.deviance %>% round(0)
    ) 
  }) %>% 
    mutate(model = model_names,
           dataset = unique(model_names) %>% str_split(pattern = '_') %>% unlist()) %>% 
    arrange(dataset) %>% 
    relocate(dataset, model)
  
  
  beta_df <- map_dfr(model_res_list, ~{
    get_beta_values(model_result = .x, selected_cov, model_nb=1)
  }) %>% 
    mutate(model = rep(unique(model_names), each = length(selected_cov))) %>% 
    pivot_wider(names_from = covar, values_from = c(beta, sd_beta), names_sep = "_") 
  
  predictive_maps_list <- map(model_res_list, 
                  function(x) make_predictive_map(x, grid = grid, selected_cov = selected_cov))
  names(predictive_maps_list) <- unique(model_names)
  
  return(list(stat_df = stat_df, beta_df = beta_df, maps = predictive_maps_list))
}


save_stat <- function(integration_res, species){
  # Create a new workbook
  workbook <- createWorkbook()
  
  # Add sheets
  addWorksheet(workbook, "perf")
  writeData(workbook, sheet = "perf", x = integration_res$stat_df)
  for (cols in 3:6) {
    conditionalFormatting(workbook, sheet = "perf", cols = cols, rows = 2:(nrow(integration_res$stat_df)+1),
                          type = "between", rule = c(0.1, 0.9), 
                          style = createStyle(bgFill = "lightgreen", fontColour = "darkgreen"))
  }
  
  addWorksheet(workbook, "beta_values")
  writeData(workbook, sheet = "beta_values", x = integration_res$beta_df)
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
  unique()

species_list <- c("labbe", "macareux_moine_HR")

for (species in species_list){
  print(species)
  integration_res <- test_data_integration(data_list, grid, species, selected_cov = best_covs[[species]])
  save_stat(integration_res, species)
  save_maps_as_pdf(integration_res$map, 
                  paste0("3.results/model_selection/maps_integration_test/", 
                         species, "_5_test_int_plot.pdf"))
}

