load("~/stage_M2/1.data/all_seabirds_counts.rdata")
load("~/stage_M2/1.data/covariates_data.rdata")

# Load package
library(tidyverse)
library(sf)
library(nimble)
source("~/stage_M2/2.code/pt2_telemetry_and_count/functions_for_nmixture.R")

# Prepare data
species <- "sterne_caugek_HR"

data <- list(obs_data = pelmed_obs, effort_data = pelmed_eff)

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs)) 

selected_cov <- c("mean_CHL", "mean_SST", "dist_to_shore")

data_nmix <- prepare_data_for_Nmixture(data, grid, species, selected_cov) 
data2 <- list(obs_data = migralion_obs %>% rename(effectif = Effectif), effort_data = migralion_eff)
data_nmix2 <- prepare_data_for_Nmixture(data2, grid, species, selected_cov)

# Nimble model
Nmixture.model <- nimbleCode({
  
  # priors
  for(i in 1:n.occ.cov){
    a[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    b[i] ~ dnorm(0,1)
  }
  
  # likelihood
  for(ds in 1:ndatasets){
    logit(p[1:nsites[[ds]],1:nreplicates[[ds]]]) <- b[1] + b[2] * transect_length[[ds]][1:nsites[[ds]], 1:nreplicates[[ds]]]
    for(i in 1:nsites[[ds]]){
      # state process
      log(lambda[site_id[i]]) <- sum(a[1:n.occ.cov] * XN[site_id[i],1:n.occ.cov])
      N[site_id[i]] ~ dpois(lambda[site_id[i]])
      # observation process
      for(j in 1:nreplicates){
        #logit(p[[ds]][i,j]) <- b[1] + b[2] * transect_length[[ds]][i, j]
        nobs[[ds]][i,j] ~ dbin(p[[ds]][i,j],N[site_id[i]])
      }
    }
  }
    
})

# Prepare data
XN = data_nmix$occurence.cov
n.occ.cov <- ncol(XN)
my.constants <- list(XN = cbind(XN$intersect, XN$mean_CHL, XN$mean_SST, XN$dist_to_shore),
                     transect_length = scale(data_nmix$detection.cov$transect_length),
                     nsites = nrow(data_nmix$effectif),
                     nreplicates = ncol(data_nmix$effectif),
                     n.occ.cov = n.occ.cov)

my.data <- list(nobs = data_nmix$effectif)

initial.values <- list(a = rnorm(n.occ.cov,0,1), 
                       b = rnorm(2,0,1), 
                       N = apply(data_nmix$effectif,1,sum)+1)

parameters.to.save <- c("a", "b", "lambda", "p", "N")

# Nombre d'iterations, burn-in et nombre de chaine
n.iter <- 100000
n.burnin <- 10000
n.chains <- 3

# In one step
mcmc.output <- nimbleMCMC(code = Nmixture.model,
                          data = my.data,
                          constants = my.constants,
                          inits = initial.values,
                          monitors = parameters.to.save,
                          niter = n.iter,
                          nburnin = n.burnin,
                          nchains = n.chains, 
                          samplesAsCodaMCMC = TRUE)

# In several step
Rmodelo <- nimbleModel(code = Nmixture.model, 
                       constants = my.constants,
                       data = my.data, 
                       inits = initial.values)

Rmodelo$initializeInfo()
Rmodelo$calculate()

# configure model
confo <- configureMCMC(Rmodelo)
#confo$printSamplers()
#confo$removeSampler(c("a[1]", "a[2]", "a[3]", "b[1]", "b[2]"))
#confo$addSampler("a[1]", "a[2]", type = "RW_block")
#confo$addSampler("b[1]", "b[2]", type = "RW_block")
## Build and compile MCMC
Rmcmco <- buildMCMC(confo)
Cmodelo <- compileNimble(Rmodelo)
Cmcmco <- compileNimble(Rmcmco, project = Cmodelo)

# Run
mcmc.output <- runMCMC(Cmcmco, 
                       niter = n.iter, 
                       nburnin = n.burnin, 
                       nchains = n.chains, 
                       samplesAsCodaMCMC = TRUE)

# check convergence
mcmcplots::denplot(mcmc.output)
coda::effectiveSize(mcmc.output)

MCMCvis::MCMCsummary(object = mcmc.output, round = 2,  params = c("a"))

MCMCvis::MCMCtrace(object = mcmc.output,
                   pdf = FALSE, 
                   ind = TRUE, 
                   Rhat = TRUE, 
                   n.eff = TRUE, 
                   params = "a")

# store results
res_b0ipp <- rbind(mcmc.output$chain1, mcmc.output$chain2) %>%
  as_tibble()%>%
  dplyr::select(starts_with("a[1]"))

res_profipp <- rbind(mcmc.output$chain1, mcmc.output$chain2) %>%
  as_tibble()%>%
  select(starts_with("a[2]"))

res_qprofipp <- rbind(mcmc.output$chain1, mcmc.output$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[3]"))

res_dcolipp <- rbind(mcmc.output$chain1, mcmc.output$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[4]"))

coeff_values <- cbind(res_profipp,res_qprofipp, res_dcolipp)

make_prediction <- function(grid, coeff_values){
  omega <- apply(coeff_values, #%>% slice_sample(size = 50000, replace = F),
                 1, function(x){x * grid2}) %>% 
    map(~{exp(rowSums(.x))})
  
  gg <- matrix(unlist(omega), ncol = length(omega))
  sd_pred <- apply(gg, 1, sd)
  mean_pred <- rowMeans(gg)
  CV_pred <- sd_pred / mean_pred
}

grid2 <- grid %>% as_tibble() %>% select(all_of(selected_cov)) 
make_prediction(grid, coeff_values)

ggplot() + geom_sf(data = grid, aes(fill = log(mean_pred))) +
  scale_fill_viridis_c()

ggplot() + geom_sf(data = grid, aes(fill = log(sd_pred))) +
  scale_fill_viridis_c(option = "B")
