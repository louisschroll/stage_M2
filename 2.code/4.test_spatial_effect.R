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

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs, eff = migralion_eff))

data.int <- get_data_for_spOccupancy(data_list, grid, species, add_coords = TRUE)

str(data.int)

det.formula <- list(pelmed = ~ scale(transect_length) + session, 
                    pnm = ~ scale(transect_length) + session,
                    samm = ~ scale(transect_length) + session,
                    migralion = ~ scale(transect_length) + session)

best_covs <- c()

without_spatial <- run_int_model(data.int, selected_cov=best_covs, det.formula, add_spatial=FALSE)
with_spatial <- run_int_model(data.int, selected_cov=best_covs, det.formula, add_spatial=TRUE)

comparison_df <- bind_rows(get_model_stat(without_spatial),
                           get_model_stat(with_spatial))

