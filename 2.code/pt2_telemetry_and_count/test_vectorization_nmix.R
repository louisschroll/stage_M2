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
source("~/stage_M2/2.code/pt2_telemetry_and_count/prepare_data_Nmix_V2.R")

# Prepare data
species <- "sterne_caugek_R"

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff), 
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff))

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs)) 

selected_cov <- c("mean_CHL", "mean_SST", "dist_to_shore")

data_nmix <- prepare_data_Nmix_V2(data_list, grid, species, selected_cov)

# Nimble model
int.Nmixture.model <- nimbleCode({
  # Priors
  for(i in 1:n.occ.cov){
    beta[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    for (nd in 1:ndatasets){
      alpha[i, nd] ~ dnorm(0,1)
    }
  }
  
  # Likelihood
  # State process
  for(i in 1:nsites_total){
    log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
    N[i] ~ dpois(lambda[i])
  }
  
  # Observation process
  for (i in 1:nsampled_points){
    logit(p[i]) <- alpha[1, dataset_nb[i]] + alpha[2, dataset_nb[i]] * transect_length[i]
    nobs[i] ~ dbin(p[i], N[site_id[i]])
  }
})


parameters.to.save <- c("beta", "alpha1","alpha2", "lambda", "p1", "N")

# Nombre d'iterations, burn-in et nombre de chaine
n.iter <- 100000
n.burnin <- 10000
n.chains <- 2

# In one step
mcmc.output <- nimbleMCMC(code = int.Nmixture.model,
                          data = data_nmix$data,
                          constants = data_nmix$constants,
                          inits = data_nmix$inits,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains,
                          samplesAsCodaMCMC = TRUE)

# In several step
Rmodelo <- nimbleModel(code = int.Nmixture.model, 
                       constants = data_nmix$constants,
                       data = data_nmix$data, 
                       inits = data_nmix$inits)

Rmodelo$initializeInfo()
Rmodelo$calculate()
# configure model
confo <- configureMCMC(Rmodelo)
#confo$printSamplers()
#confo$removeSampler(c("a[1]", "a[2]", "a[3]", "b[1]", "b[2]"))
#confo$addSampler("a[1]", "a[2]", type = "RW_block")
#confo$addSampler("b[1]", "b[2]", type = "RW_block")
## Build and compile MCMC
Rmcmco <- buildMCMC(confo, monitors = parameters.to.save)
Cmodelo <- compileNimble(Rmodelo)
Cmcmco <- compileNimble(Rmcmco, project = Cmodelo)

# Run
mcmc.output <- runMCMC(Cmcmco, 
                       niter = n.iter, 
                       nburnin = n.burnin, 
                       nchains = n.chains, 
                       samplesAsCodaMCMC = TRUE)

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





