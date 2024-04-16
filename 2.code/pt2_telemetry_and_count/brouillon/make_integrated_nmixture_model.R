# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/make_integrated_nmixture_model.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-03-28
#
# Script description:
#
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

load("~/stage_M2/1.data/all_seabirds_counts.rdata")
load("~/stage_M2/1.data/covariates_data.rdata")

# Load package
library(tidyverse)
library(sf)
library(nimble)
source("~/stage_M2/2.code/pt2_telemetry_and_count/prepare_data_Nmix.R")

# Prepare data
species <- "sterne_caugek_R"

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff), 
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff))

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs)) 

selected_cov <- c("mean_CHL", "mean_SST", "dist_to_shore")

data_nmix <- prepare_data_Nmix(data_list, grid, species, selected_cov)



# check convergence
mcmcplots::traplot(mcmc.output)
mcmcplots::denplot(mcmc.output)
coda::effectiveSize(mcmc.output)

MCMCvis::MCMCsummary(object = mcmc.output, round = 2,  params = c("beta"))

new_grid <- make_prediction(mcmc.output, grid, selected_cov)

plots <- plot_prediction(new_grid, add_colonies = T, species_colony = "Sterne caugek")
library(patchwork)
plots$mean_psi_plot + plots$sd_psi_plot


samples <- do.call(rbind, mcmc.output)    # single matrix of samples
waic <- calculateWAIC(samples, Rmodelo)
print(waic)



