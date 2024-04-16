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

source("2.code/pt1_spOccupancy/prediction_functions.R")
source("2.code/pt1_spOccupancy/plot_functions.R")
grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs))


data_list = list(
  pelmed = list(obs = pelmed_obs, eff = pelmed_eff),
  samm = list(obs = samm_obs, eff = samm_eff),
  pnm = list(obs = pnm_obs, eff = pnm_eff),
  migralion = list(obs = migralion_obs, eff = migralion_eff)
  )

species <- "oceanite_tempete"

data.int <- get_data_for_spOccupancy(data_list, grid, species)

selected_cov <- c("dist_to_shore", "log_bathymetry", "sd_SSH", "log_sd_VEL")
aa <- run_model_without_kfold(data.int = data.int, grid=grid, species=species, selected_cov=selected_cov)
mappp <- make_predictive_map(aa, grid, selected_cov)
mappp$psi 

fff <- make_coeff_plot(aa)
fff$beta_traceplot
grid_caugek <- put_results_in_grid(grid, model_result = aa$res, selected_cov) %>% 
  select(mean.intensity, sd.intensity)
plot(grid_caugek)
save(grid_caugek, file = "results_sterne_caugek_hors_repro.rdata")
  # geom_sf(data = migralion_obs %>% filter(species_name == species)) +
  # geom_sf(data = pelmed_obs %>% filter(species_name == species)) +
  # geom_sf(data = pnm_obs %>% filter(species_name == species)) +
  # geom_sf(data = migralion_eff) + 
  # geom_sf(data = pelmed_eff)
  #geom_sf(data = countLL %>% filter(Species == "Sterne caugek"))
  #geom_sf(data = pnm_eff) +
  #geom_sf(data = samm_eff) 
 
# ggplot(migralion_obs %>% filter(nom_fr == "sterne pierregarin")) +
#   geom_sf() +
#   facet_wrap(~session) +
#   geom_sf(data = contour_golfe) +
#   geom_sf(data = migralion_eff)
# 
# ggplot(pelmed_obs %>% filter(nom_fr == "sterne pierregarin" & session != "2020")) +
#   geom_sf() +
#   facet_wrap(~session) +
#   geom_sf(data = contour_golfe) +
#   geom_sf(data = pelmed_eff %>% filter(session != "2020"))
  
waicOcc(aa$res)

pavlue_CS1 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 1)
pavlue_CS2 <- ppcOcc(aa$res, fit.stat = "chi-squared", group = 2)
pavlue_FT1 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 1)
pavlue_FT2 <- ppcOcc(aa$res, fit.stat = "freeman-tukey", group = 2)

ppc.out <- list(fit.y = list(pavlue_CS1$fit.y, pavlue_CS2$fit.y, pavlue_FT1$fit.y, pavlue_FT2$fit.y),
            fit.y.rep = list(pavlue_CS1$fit.y.rep, pavlue_CS2$fit.y.rep, pavlue_FT1$fit.y.rep, pavlue_FT2$fit.y.rep))
# ii = 3
# ppc.out <- list(fit.y = list(pavlue_CS1$fit.y[[ii]], pavlue_CS2$fit.y[[ii]],
#                              pavlue_FT1$fit.y[[ii]], pavlue_FT2$fit.y[[ii]]),
#                 fit.y.rep = list(pavlue_CS1$fit.y.rep[[ii]], pavlue_CS2$fit.y.rep[[ii]],
#                                  pavlue_FT1$fit.y.rep[[ii]], pavlue_FT2$fit.y.rep[[ii]]))
compute_bayesian_pvalue(ppc.out)
plot_PPC(ppc.out)


# ------ graphe livrable
library(cowplot)

species <- "sterne_caugek_HR"

data.int <- get_data_for_spOccupancy(data_list, grid, species)

selected_cov_sterne <- c("log_dist_to_shore", "log_bathymetry",
                  "mean_winter_SST", "mean_spring_SST", "mean_summer_SST",
                  "mean_CHL", "mean_SSH")
aa <- run_and_predict(data.int = data.int, grid=grid, species=species, selected_cov=selected_cov_sterne)

aa$psi 



species <- "petit_puffin_HR"

data.int <- get_data_for_spOccupancy(data_list, grid, species)

selected_cov_puffin <- c("mean_CHL", "sd_SAL", "mean_SSH", "mean_winter_SST", "mean_autumn_SST")
bb <- run_and_predict(data.int = data.int, grid=grid, species=species, selected_cov=selected_cov_puffin)

bb$psi 

species <- "mouette_pygmee_HR"

data.int <- get_data_for_spOccupancy(data_list, grid, species)

selected_cov_mouette <- c("log_dist_to_shore", "log_bathymetry", "sd_SAL", "log_sd_VEL")
cc <- run_and_predict(data.int = data.int, grid=grid, species=species, selected_cov=selected_cov_mouette)

cc$psi 


size = 13

sterne <- plot_grid(plotlist = make_predictive_map(aa$res, grid, selected_cov_sterne), 
                    labels = "Sterne caugek", 
                    label_y = 1.02, label_x = 0, scale = 0.9, label_size = size)

puffin <- plot_grid(plotlist = make_predictive_map(bb$res, grid, selected_cov_puffin), 
                    labels = "Petits puffins", 
                    label_y = 1.02,label_x = 0 , scale = 0.9, label_size = size)

mouette <- plot_grid(plotlist = make_predictive_map(cc$res, grid, selected_cov_mouette), 
                     labels = "Mouette pygmée", 
                     label_y = 1.02, label_x = 0, scale = 0.9, label_size = size)


plot_final <- plot_grid(plotlist = list(sterne, mouette, puffin), nrow = 3)
plot_final
