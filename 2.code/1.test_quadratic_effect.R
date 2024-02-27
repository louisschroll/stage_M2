# HEADER ------------------------------------------------------------------------
#
# Script name:  1.test_quadratic_effect.R
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

test_quad_and_log <- function(data.int, species){
  
  cov_list <- as.list(c("mean_winter_SST", "mean_spring_SST",
                        "mean_summer_SST", "mean_autumn_SST",
                        "sd_winter_SST", "sd_spring_SST",
                        "sd_summer_SST", "sd_autumn_SST",
                        "mean_SST", "sd_SST",
                        "concavity", "slope", "dist_to_shore",
                        "bathymetry",
                        "mean_CHL", "mean_VEL", "sd_VEL"))

  models_to_test <- map(cov_list, 
                        function(x) list("1", x,  c(x, paste0("I(",x,")^2")), paste0("log_",x)))
  
  test_and_write_models <- function(cov_combination, data.int){
    df <- test_all_models(cov_combination, data.int)
    addSelectionSheet(wb, sheet_name = cov_combination[[2]], df = df, datasets_nb=length(data.int$y))
  }
  # Create a new workbook
  wb <- createWorkbook()
  
  # Fill the workbook
  map(models_to_test, ~ test_and_write_models(cov_combination = .x, data.int = data.int))
  
  # Save the workbook
  file_path <- paste0("3.results/model_selection/", 
                      str_replace(species, " ","_"), 
                      "_1_test_log_quad.xlsx")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  print(paste("results are in", file_path))
}

# ----- Hors repro -----
species_list <- c("sterne caugek", "mouette pygmee", "goeland leucophee", "petit puffin",
                   "mouette melanocephale", "puffin de scopoli", "oceanite tempete")

migralion_obs2 <- migralion_obs %>% filter(session != "prenup_2022")
migralion_eff2 <- migralion_eff %>% filter(session != "prenup_2022")

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs2, eff = migralion_eff2))


for (species in species_list){
  print(species)
  data.int <- get_data_for_spOccupancy(data_list, grid, species)
  test_quad_and_log(data.int, species)
}

 