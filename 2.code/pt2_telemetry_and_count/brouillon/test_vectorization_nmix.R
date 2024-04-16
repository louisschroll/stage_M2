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
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(nimble, warn.conflicts = FALSE)
source("~/stage_M2/2.code/pt2_telemetry_and_count/prepare_data_Nmix.R")

# Prepare data
species <- "sterne_caugek_R"

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff),
                  #migralion = list(obs_data = migralion_obs, effort_data = migralion_eff),
                  pnm = list(obs_data = pnm_obs, effort_data = pnm_eff)
                  )

grid <- covariates_data %>% 
  mutate(cells_id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs)) 

selected_cov <- c("mean_CHL", "mean_SST", "dist_to_shore")

data_nmix <- prepare_data_Nmix(data_list, grid, species, selected_cov)

# Nimble code
int.Nmixture.code <- nimbleCode({
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


parameters.to.save <- c("beta", "alpha", "lambda", "p", "N")

# Nombre d'iterations, burn-in et nombre de chaine
n.iter <- 100000
n.burnin <- 10000
n.chains <- 2

# In one step
# mcmc.output <- nimbleMCMC(code = int.Nmixture.model,
#                           data = data_nmix$data,
#                           constants = data_nmix$constants,
#                           inits = data_nmix$inits,
#                           monitors = parameters.to.save,
#                           niter = n.iter,
#                           nburnin = n.burnin,
#                           nchains = n.chains,
#                           samplesAsCodaMCMC = TRUE)

# In several step
int.Nmixture.model <- nimbleModel(code = int.Nmixture.code, 
                       constants = data_nmix$constants,
                       data = data_nmix$data, 
                       inits = data_nmix$inits)

int.Nmixture.model$initializeInfo()
int.Nmixture.model$calculate()
# configure model
confo <- configureMCMC(int.Nmixture.model, monitors = parameters.to.save)
#confo$printSamplers()
#confo$removeSampler(c("a[1]", "a[2]", "a[3]", "b[1]", "b[2]"))
#confo$addSampler("a[1]", "a[2]", type = "RW_block")
#confo$addSampler("b[1]", "b[2]", type = "RW_block")
## Build and compile MCMC
Rmcmco <- buildMCMC(confo)
Cmodelo <- compileNimble(int.Nmixture.model)
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


# PPC

## Ensure we have the nodes needed to simulate new datasets
dataNodes <- int.Nmixture.model$getNodeNames(dataOnly = TRUE)
parentNodes <- int.Nmixture.model$getParents(dataNodes, stochOnly = TRUE)  # `getParents` is new in nimble 0.11.0
## Ensure we have both data nodes and deterministic intermediates (e.g., lifted nodes)
simNodes <- int.Nmixture.model$getDependencies(parentNodes, self = FALSE)


ppSamplerNF <- nimbleFunction(
  setup = function(model, mcmc) {
    dataNodes <- model$getNodeNames(dataOnly = TRUE)
    parentNodes <- model$getParents(dataNodes, stochOnly = TRUE)
    cat("Stochastic parents of data are:", paste(parentNodes, collapse = ','), ".\n")
    simNodes <- model$getDependencies(parentNodes, self = FALSE)
    vars <- mcmc$mvSamples$getVarNames()  # need ordering of variables in mvSamples / samples matrix
    cat("Using posterior samples of:", paste(vars, collapse = ','), ".\n")
    n <- length(model$expandNodeNames(dataNodes, returnScalarComponents = TRUE))
  },
  run = function(samples = double(2)) {
    nSamp <- dim(samples)[1]
    ppSamples <- matrix(nrow = nSamp, ncol = n)   
    for(i in 1:nSamp) {
      values(model, vars) <<- samples[i, ]
      model$simulate(simNodes, includeData = TRUE)
      ppSamples[i, ] <- values(model, dataNodes)
    }
    returnType(double(2))       
    return(ppSamples)
  })

## Create the sampler for this model and this MCMC.
ppSampler <- ppSamplerNF(int.Nmixture.model, Cmcmco)

cppSampler <- compileNimble(ppSampler, project = int.Nmixture.model)

colnames(samples)  
identical(colnames(samples), int.Nmixture.model$expandNodeNames(Cmcmco$mvSamples$getVarNames()))

ppSamples <- cppSampler$run(samples)

obsMin <- min(data_nmix$data$nobs)
ppMin <- apply(ppSamples, 1, min)
dim(ppSamples)
dim(samples)
# ## Check with plot in Gelman et al. (3rd edition), Figure 6.3
hist(ppMin, #xlim = c(-50, 20),
     main = "Discrepancy = min(y)", 
     xlab = "min(y_rep)")
abline(v = obsMin, col = 'red')
