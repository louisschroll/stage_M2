# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/make_RSF_model.R
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

library(nimble)

load("1.data/data_RSF_caugek.rdata")
#### RSF #####

##### code #### 
rsf.caugek <- nimbleCode({
  # priors
  betaprof_pop ~ dnorm(0,1)
  tauprof_pop ~ dunif(0,1e2)
  betaqprof_pop ~ dnorm(0,1)
  tauqprof_pop ~ dunif(0,1e2)
  betadcol_pop ~ dnorm(0,1)
  taudcol_pop ~ dunif(0,1e2)
  
  beta0_pop ~ dnorm(0,1)
  
  for( i in 1:nindividual){
    # habitat cov
    beta_prof[i] ~ dnorm(betaprof_pop, sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop, sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop, sd = taudcol_pop)
    # movement cov
    beta_0[i] ~ dnorm(beta0_pop,1e2)
  }
  
  # likelihood
  for(t in 1:npoints){
    
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * prof[t] + 
      beta_qprof[idind[t]] *prof [t] * prof[t] +
      beta_dcol[idind[t]] * dcol[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
})


# --- bundle and run ----
weight <- df_RSF$case
weight[weight==0] <- 1000
# constants
nindividual = n_distinct(df_RSF$individual_id)
constants.ni <-  list(npoints = nrow(df_RSF),
                      idind = df_RSF$individual_id,
                      prof = df_RSF$log_bathymetry,
                      dcol = df_RSF$dist_to_shore,
                      nindividual = nindividual,
                      w = weight)

# data
data.ni <- list(kase = df_RSF$case)

# Inits
inits.ni <-  list(beta_0 = rep(0, nindividual),
                  beta_prof =  rep(0, nindividual),
                  beta_qprof=  rep(0, nindividual),
                  beta_dcol =  rep(0, nindividual),
                  betaprof_pop = 0 ,
                  beta0_pop = 0 ,
                  tauprof_pop = runif(1,0,5),
                  betaqprof_pop = 0 ,
                  tauqprof_pop = runif(1,0,5),
                  betadcol_pop = 0 ,
                  taudcol_pop = runif(1,0,5))


# Nimble pre run
Rmodel2 <- nimbleModel(code= rsf.caugek, constants = constants.ni, data = data.ni, inits = inits.ni)
Rmodel2$initializeInfo()
Rmodel2$calculate() # - 58854147

# configure model
conf2 <- configureMCMC(Rmodel2)
conf2$printMonitors()
#
## Build and compile MCMC
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)


# Run
samplesRSF <- runMCMC(Cmcmc2, niter = 11000, nburnin = 1000, thin = 1, nchains = 1, samplesAsCodaMCMC = TRUE)  ## DT: use runMCMC

# check convergence
mcmcplots::traplot(samplesRSF)
mcmcplots::denplot(samplesRSF)
coda::effectiveSize(samplesRSF)

# store results
# when only 1 chain
res_prof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaprof_pop"))

res_tauprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("tauprof_pop"))

res_qprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaqprof_pop"))

res_dcol <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betadcol_pop"))

# make prediction maps
if(F){

gridpred <- grid 
dist_to_coast <- unlist(map2(gridpred$x, st_union(gdlmap), st_distance))


pred_bathy <- tibble(prof =  grid$log_bathymetry, dcoast = grid$dist_to_shore) %>%
  mutate(pred = prof* mean(res_prof$betaprof_pop) +
            dcoast* mean(res_dcol$betadcol_pop))

ggplot(data = pred_bathy) +
  geom_point(aes( x = prof, y = pred))

ggplot() + geom_sf(data = grid, aes(fill = pred_bathy$pred))
}
