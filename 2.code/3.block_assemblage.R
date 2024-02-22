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
source("~/stage_M2/2.code/generateAllCombinations.R")
source("~/stage_M2/2.code/model_selection_functions.R")

# load data
load("~/stage_M2/1.data/all_seabirds_counts.rdara")
load("~/stage_M2/1.data/covariates_data.rdata")

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

migralion_obs2 <- migralion_obs #%>% filter(session != "prenup_2022")
migralion_eff2 <- migralion_eff #%>% filter(session != "prenup_2022")

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs2, eff = migralion_eff2))

data.int <- get_data_for_spOccupancy(data_list, grid, species)


best_static_cov <- c("bathymetry")
best_dyn_cov <- c("winter_SST")
best_SST_cov <- c("mean_CHL")

models_to_test <- list(best_static_cov, 
                      best_dyn_cov,
                      best_SST_cov,
                      c(best_static_cov, best_dyn_cov), 
                      c(best_static_cov, best_SST_cov),
                      c(best_dyn_cov, best_SST_cov),
                      c(best_static_cov, best_dyn_cov, best_SST_cov))

df_grouping <- test_all_models(models_to_test, data.int, det.formula)

