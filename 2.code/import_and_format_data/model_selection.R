# HEADER ------------------------------------------------------------------------
#
# Script name:  
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-21
#
# Script description:
#
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(sf)
library(spOccupancy)
library(openxlsx)

# Load functions
source("~/stage_M2/2.code/format_data_for_spOccupancy.R")
source("~/stage_M2/2.code/generateAllCombinations.R")
source("~/stage_M2/2.code/model_selection_functions.R")

# load data
load("~/stage_M2/1.data/all_seabirds_counts.rdara")
load("~/stage_M2/1.data/covariates_data.rdata")
grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

## ------------- Species: Sterne caugek ---------------
species <- "sterne caugek"

data_list = list(pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
                 samm = list(obs = samm_obs, eff = samm_eff),
                 pnm = list(obs = pnm_obs, eff = pnm_eff),
                 migralion = list(obs = migralion_obs2, eff = migralion_eff2))

data.int <- get_data_for_spOccupancy(data_list, grid, species)

det.formula <- list(pelmed = ~ scale(transect_length) + session, 
                    pnm = ~ scale(transect_length) + session,
                    samm = ~ scale(transect_length) + session,
                    migralion = ~ scale(transect_length) + session)
# Step 1: test log and quadratic effect

# Step 2: test covariates by block

# Step 3: test adding the block together

# Step 4: test adding spatial autocorrelation

# Step 5: assess integration relevance



