
# Load package
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(patchwork)

species_name <- "labbe"
load(paste0("~/stage_M2/3.results/mcmc_outputs/Nmix_output_", species_name, ".rdata"))

# check convergence
mcmcplots::traplot(samplesNmixture, parms = c("alpha", "beta"))
mcmcplots::denplot(samplesNmixture, parms = c("alpha", "beta"))

MCMCvis::MCMCsummary(object = samplesNmixture, round = 2, params = c("alpha", "beta"))

# Look at parameter values
MCMCvis::MCMCplot(object = samplesNmixture, params = c("beta"))


# Compare coefficients for occupancy and N-mix
load("~/stage_M2/3.results/spOccupancy_outputs/spOccupancy_fou_de_bassan_HR.RData")
load("1.data/grid_cells.rdata")

samples_spOcc <- model_result$beta.samples

selected_cov <- c("dist_to_shore", "log_bathymetry")
coeff_df <- gather_coeff_values(sampleNmixture = sampleNmixture, samples_spOcc = samples_spOcc, selected_cov)

coeff_df <- coeff_df %>% 
  mutate(covariates = case_when(
    covariates == "dist_to_shore" ~ "Distance to coast",
    covariates == "log_bathymetry" ~ "Bathymetry"
  ))

plot_coeff_by_model2(coeff_df)


# Look at the maps
load("~/stage_M2/3.results/grid_spOccupancy_results.RData")
grid_occupancy <- spOccupancy_res_grid %>% select(all_of(c("grid_c", "mean_psi_fou_de_bassan_HR", "sd_psi_fou_de_bassan_HR")))
names(grid_occupancy) <- c("mean_psi", "sd_psi", "grid_c")
occupancy_plot <- plot_occupancy(grid_occupancy = grid_occupancy)

grid_nmix <- make_prediction(samplesNmixture, grid = grid, selected_cov = selected_cov)
plots_nmix <- plot_prediction(grid_nmix, plot_title = "N-mixture")

(occupancy_plot$mean_psi_plot + occupancy_plot$sd_psi_plot) / 
  (plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot)

