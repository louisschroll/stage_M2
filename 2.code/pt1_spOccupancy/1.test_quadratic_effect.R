#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  1.test_quadratic_effect.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-02
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(sf)
library(spOccupancy)
library(openxlsx)

# cluster adress
# adress <- "/lustre/schrolll/"
# Local adress
adress <- ""

# Load functions
source(paste0(adress, "2.code/pt1_spOccupancy/format_data_for_spOccupancy.R"))
source(paste0(adress, "2.code/pt1_spOccupancy/model_selection_functions.R"))

# Load data
load(paste0(adress, "1.data/all_seabirds_counts.rdata"))
load(paste0(adress, "1.data/grid.rdata"))

# Function to test quadratic and log effects
test_quad_and_log <- function(data.int, species){
  
  cov_list <- as.list(c("mean_winter_SST", "mean_spring_SST", "mean_summer_SST", 
                        "mean_autumn_SST", "mean_SST", "sd_SST", "concavity", 
                        "dist_to_shore", "bathymetry",
                        "mean_CHL", "mean_SSH", "sd_SSH", "sd_VEL"))

  models_to_test <- map(cov_list, 
                        function(x) list(x,  c(x, paste0("I(",x,")^2")), paste0("log_",x)))
  
  test_and_write_models <- function(cov_combination, data.int, df_null_model){
    df <- test_all_models(cov_combination, data.int)$comparison_df
    addSelectionSheet(wb, sheet_name = cov_combination[[1]], 
                      df = bind_rows(df_null_model, df), 
                      datasets_nb=length(data.int$y))
  }
  # Create a new workbook
  wb <- createWorkbook()
  # Fill the workbook
  df_null_model <- test_all_models("1", data.int)$comparison_df
  map(models_to_test, ~ test_and_write_models(cov_combination = .x, data.int = data.int, df_null_model))
  
  # Save the workbook
  file_path <- paste0(adress, "3.results/model_selection/", 
                      str_replace(species, " ","_"), 
                      "_1_test_log_quad.xlsx")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  print(paste("results are in", file_path))
}

# Get species names
# species_list <- migralion_obs %>%
#   filter(!is.na(species_name)) %>%
#   pull(species_name) %>%
#   unique() %>% 
#   str_subset("macareux", negate = F)
species_list <- c("labbe", "macareux_moine_HR")
# Write data in a list
data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

# Test quadratic and log effect for all species
for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_quad_and_log(data.int, species)
}

 