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

str(data.int)

det.formula <- list(pelmed = ~ scale(transect_length) + session, 
                    pnm = ~ scale(transect_length) + session,
                    samm = ~ scale(transect_length) + session,
                    migralion = ~ scale(transect_length) + session)

# Static covariates selection
static_cov <- c("dist_to_shore", "concavity", "slope", "bathymetry")
cov_combination <- generateAllCombinations(static_cov) %>% 
  addQuadraticEffect("bathymetry")
df_static_cov <- test_all_models(cov_combination, data.int, det.formula)

# temperature cov selection
SST_cov <- c("winter_SST", "spring_SST", "summer_SST", "autumn_SST")
SST_combination <- generateAllCombinations(SST_cov) %>% c(list("mean_SST"))
df_SST_cov <- test_all_models(SST_combination, data.int, det.formula)

# Dynamic cov selection
dyn_cov <- c("mean_CHL")
dyn_cov_combination <- generateAllCombinations(dyn_cov) 
df_dyn_cov <- test_all_models(dyn_cov_combination, data.int, det.formula)

# Create a new workbook
wb <- createWorkbook()
addSelectionSheet(wb, sheet_name = "static_cov", df = df_static_cov)
addSelectionSheet(wb, sheet_name = "SST_cov", df = df_SST_cov)
addSelectionSheet(wb, sheet_name = "dynamic_cov", df = df_dyn_cov)

# Save the workbook
saveWorkbook(wb, paste0("covariates_selection_", str_replace(species, " ","_"), ".xlsx"), overwrite = TRUE)


