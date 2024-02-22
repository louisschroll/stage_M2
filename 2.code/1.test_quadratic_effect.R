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


species <- "sterne caugek"
species_list <- c("sterne caugek", "goeland leucophee", "puffin yelkouan",
                  "mouette melanocephale", "puffin de scopoli", "oceanite tempete")

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

cov_list <- as.list(c("mean_winter_SST", "mean_spring_SST", 
                      "mean_summer_SST", "mean_autumn_SST", 
                      #"sd_winter_SST", "sd_spring_SST", 
                      #"sd_summer_SST", "sd_autumn_SST", 
                      "mean_SST", "sd_SST", 
                      "concavity", "slope", "dist_to_shore", 
                      #"bathymetry", 
                      "mean_CHL", "mean_VEL", "sd_VEL"))

get_model_list <- function(x){
  list("1", x,  c(x, paste0("I(",x,")^2")), paste0("log_",x))
}
model_to_test <- map(cov_list, get_model_list)

# Create a new workbook
wb <- createWorkbook()
for (i in seq_along(model_to_test)){
  print(paste("Covariates", i, "/", length(model_to_test)))
  df <- test_all_models(model_to_test[[i]], data.int, det.formula)
  addSelectionSheet(wb, sheet_name = cov_list[[i]], df = df)
}

# Save the workbook
saveWorkbook(wb, paste0("quadratic_effect_", str_replace(species, " ","_"), ".xlsx"), overwrite = TRUE)


