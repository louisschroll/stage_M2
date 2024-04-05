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


load("1.data/RSF_data_sandwich_tern.rdata")
df_RSF <- df_RSF %>% filter(individual_id %in% 1:3)

# Nimble model
int.Nmix.RSF.model <- nimbleCode({
  # Priors
  beta0_nmix ~ dnorm(0,1) # different intercepts for RSF and Nmix. RSF intercept is beta_pop[1]
  for(i in 1:n.occ.cov){
    beta_pop[i] ~ dnorm(0,1)
    sd_pop[i] ~ dunif(0,1e2)
    
    for(j in 1:nindividual){
      # add individual heterogeneity
      beta_ind[j, i] ~ dnorm(beta_pop[i], sd = sd_pop[i])
    }
  }
  
  for(i in 1:2){
    alpha1[i] ~ dnorm(0,1)
    alpha2[i] ~ dnorm(0,1)
    alpha3[i] ~ dnorm(0,1)
  }
  
  # likelihood
  # RSF
  for(t in 1:npoints){
    logit(omega[t]) <- sum(beta_ind[idind[t], 1:n.occ.cov] * XN_rsf[t, 1:n.occ.cov])
    kase[t] ~ dbinom(omega[t], w[t])
  }
  # N mixture
  # state process
  for(i in 1:nsites_total){
    log(lambda[i]) <- beta0_nmix + sum(beta_pop[2:n.occ.cov] * XN_nmix[i,2:n.occ.cov])
    N[i] ~ dpois(lambda[i])
  }
  ## observation process
  ### Data source 1
  logit(p1[1:nsites1,1:nreplicates1]) <- alpha1[1] + alpha1[2] * transect_length1[1:nsites1, 1:nreplicates1]
  for(i in 1:nsites1){
    for(j in 1:nreplicates1){
      nobs1[i,j] ~ dbin(p1[i,j], N[site_id1[i]])
    }
  }
  ### Data source 2
  logit(p2[1:nsites2,1:nreplicates2]) <- alpha2[1] + alpha2[2] * transect_length2[1:nsites2, 1:nreplicates2]
  for(i in 1:nsites2){
    for(j in 1:nreplicates2){
      nobs2[i,j] ~ dbin(p2[i,j], N[site_id2[i]])
    }
  }
  ### Data source 3
  logit(p3[1:nsites3,1:nreplicates3]) <- alpha3[1] + alpha3[2] * transect_length3[1:nsites3, 1:nreplicates3]
  for(i in 1:nsites3){
    for(j in 1:nreplicates3){
      nobs3[i,j] ~ dbin(p3[i,j], N[site_id3[i]])
    }
  }
})



rsf_list <- format_rsf_data_for_nimble(df_RSF, selected_cov)

int.inits <- c(rsf_list$inits, data_nmix$inits[2:length(data_nmix$inits)],
               list(beta0_nmix = rnorm(1,0,1)))

names(rsf_list$constants)[names(rsf_list$constants) == "XN"] <- "XN_rsf"
rsf_list$constants %>% str()
data_nmix$constants[which(names(data_nmix$constants)%in%c("n.occ.cov"))] <- NULL
names(data_nmix$constants)[names(data_nmix$constants)=="XN"] <- "XN_nmix"
data_nmix$constants %>% str()

int.constants <- c(rsf_list$constants, data_nmix$constants)

int.data <- c(rsf_list$data, data_nmix$data)
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

# Nimble pre run
int.model <- nimbleModel(code = int.Nmix.RSF.model, 
                         constants = int.constants, 
                         data = int.data, 
                         inits = int.inits)
int.model$initializeInfo()
int.model$calculate()

# configure model
model_configuration <- configureMCMC(int.model)
model_configuration$printMonitors()

# Build and compile MCMC
Rmcmc <- buildMCMC(model_configuration)
Cmodel <- compileNimble(int.model)
Cmcmc <- compileNimble(Rmcmc, project = Cmodel)

# Run
samplesint <- runMCMC(Cmcmc, 
                      niter = n.iter, 
                      nburnin = n.burnin, 
                      thin = 1, 
                      nchains = n.chains, 
                      samplesAsCodaMCMC = TRUE)

# check convergence
mcmcplots::traplot(samplesint)
mcmcplots::denplot(samplesint)
coda::effectiveSize(samplesint)

# make prediction maps
new_grid <- make_prediction(mcmc.output = samplesint, grid, selected_cov)

plot <- plot_prediction(new_grid)
plot$mean_psi_plot + geom_point(data = df_RSF %>% filter(case==1), aes(x=X, y=Y), alpha = 0.1)
plot$sd_psi_plot
