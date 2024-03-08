# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/predictions_tests.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-28
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
library(patchwork)

# load data
load("1.data/all_seabirds_counts.rdata")
load("1.data/covariates_data.rdata")
load("1.data/countLL.rdata")

source("2.code/prediction_functions.R")
source("2.code/plot_functions.R")
grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))


data_list = list(
   pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
   samm = list(obs = samm_obs, eff = samm_eff),
   #pnm = list(obs = pnm_obs, eff = pnm_eff),
  migralion = list(obs = migralion_obs, eff = migralion_eff)
  )

species <- "petit_puffin_R"

data.int <- get_data_for_spOccupancy(data_list, grid, species)
# "mean_winter_SST", "log_dist_to_shore", "log_sd_VEL", "mean_CHL"
selected_cov <- c("sd_SSH",
                  #"sd_SAL", 
                  "mean_SSH", 
                  "sd_VEL"
                   )
aa <- predict_distribution(data.int, grid=grid, species=species, selected_cov=selected_cov)

aa$psi +
  # geom_sf(data = migralion_obs %>% filter(species_name == species)) +
  # geom_sf(data = pelmed_obs %>% filter(species_name == species)) +
  # geom_sf(data = pnm_obs %>% filter(species_name == species)) +
  # geom_sf(data = migralion_eff) + 
  # geom_sf(data = pelmed_eff)
  #geom_sf(data = countLL %>% filter(Species == "Sterne caugek"))
  #geom_sf(data = pnm_eff) +
  #geom_sf(data = samm_eff) +
 

waicOcc(aa$res)

pavlue_CS1 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 1)
pavlue_CS2 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 2)
pavlue_FT1 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 1)
pavlue_FT2 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 2)

# ppc.out <- list(fit.y = list(pavlue_CS1$fit.y, pavlue_CS2$fit.y, pavlue_FT1$fit.y, pavlue_FT2$fit.y),
#             fit.y.rep = list(pavlue_CS1$fit.y.rep, pavlue_CS2$fit.y.rep, pavlue_FT1$fit.y.rep, pavlue_FT2$fit.y.rep))
ii = 1
ppc.out <- list(fit.y = list(pavlue_CS1$fit.y[[ii]], pavlue_CS2$fit.y[[ii]],
                             pavlue_FT1$fit.y[[ii]], pavlue_FT2$fit.y[[ii]]),
                fit.y.rep = list(pavlue_CS1$fit.y.rep[[ii]], pavlue_CS2$fit.y.rep[[ii]],
                                 pavlue_FT1$fit.y.rep[[ii]], pavlue_FT2$fit.y.rep[[ii]]))
compute_bayesian_pvalue(ppc.out)
plot_PPC(ppc.out)

