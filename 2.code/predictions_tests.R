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
load("1.data/all_seabirds_counts.rdara")
load("1.data/covariates_data.rdata")

source("2.code/prediction_functions.R")
source("2.code/plot_functions.R")
grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))

pnm_obs2 <- pnm_obs %>% filter(session %in% c("2020_printemps", "2021_printemps"))
pnm_eff2 <- pnm_eff %>% filter(session %in% c("2020_printemps", "2021_printemps"))

data_list = list(
  # pelmed = list(obs = pelmed_obs2, eff = pelmed_eff2)
  # samm = list(obs = samm_obs, eff = samm_eff),
   pnm = list(obs = pnm_obs2, eff = pnm_eff2)
 # migralion = list(obs = migralion_obs2, eff = migralion_eff2)
  )

species <- "goeland_leucophee_HR"

data.int <- get_data_for_spOccupancy(data_list, grid, species)
# "mean_winter_SST", "log_dist_to_shore", "log_sd_VEL", "mean_CHL"
selected_cov <- c("log_dist_to_shore", "log_bathymetry", "mean_CHL",
                  "mean_winter_SST", "mean_spring_SST", "mean_summer_SST", "sd_SST")
aa <- predict_distribution(data.int, grid=grid, species=species, selected_cov=selected_cov)
aa$psi
waicOcc(aa$res)

pavlue_CS1 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 1)
pavlue_CS2 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 2)
pavlue_FT1 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 1)
pavlue_FT2 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 2)

# ppc.out <- list(fit.y = list(pavlue_CS1$fit.y, pavlue_CS2$fit.y, pavlue_FT1$fit.y, pavlue_FT2$fit.y),
#            fit.y.rep = list(pavlue_CS1$fit.y.rep, pavlue_CS2$fit.y.rep, pavlue_FT1$fit.y.rep, pavlue_FT2$fit.y.rep))
ii = 4
ppc.out <- list(fit.y = list(pavlue_CS1$fit.y[[ii]], pavlue_CS2$fit.y[[ii]],
                             pavlue_FT1$fit.y[[ii]], pavlue_FT2$fit.y[[ii]]),
                fit.y.rep = list(pavlue_CS1$fit.y.rep[[ii]], pavlue_CS2$fit.y.rep[[ii]],
                                 pavlue_FT1$fit.y.rep[[ii]], pavlue_FT2$fit.y.rep[[ii]]))
compute_bayesian_pvalue(ppc.out)
plot_PPC(ppc.out)

effect_strength2 <- aa$res$beta.samples %>% as_tibble() %>% 
  select(all_of(selected_cov)) %>% 
  summarise_all(mean) %>% 
  pivot_longer(all_of(selected_cov), names_to = "covar", values_to = "beta") %>% 
  mutate(model = "M2")

sd_beta <- aa$res$beta.samples %>% as_tibble() %>% 
  select(all_of(selected_cov)) %>% 
  summarise_all(sd) %>% 
  pivot_longer(all_of(selected_cov), names_to = "covar", values_to = "sd_beta") %>% 
  mutate(model = "M2")

effect_strength2 %>% 
  full_join(sd_beta, by = join_by(covar, model))

effect_strength %>% bind_rows(effect_strength2) %>% 
  complete(cov, model, fill = list(beta = NA)) %>% 
  select(cov, covar, beta) %>% 
  pivot_wider(names_from = cov, values_from = c(covar, beta))
             