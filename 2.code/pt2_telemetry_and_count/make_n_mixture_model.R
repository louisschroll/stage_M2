# HEADER ------------------------------------------------------------------------
#
# Script name:  
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-03-21
#
# Script description:
#
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# load package
library(nimble)

#### N-mixture model ####

# fit in Bayesian to be consistent

# format data
constants.o <- list(prof = scale(ypel$depth.x)[,1],
                    dcol = scale(ypel$dist_coast)[,1],
                    seff = scale(cbind(ypel$eff2017,
                                       ypel$eff2018,
                                       ypel$eff2019,
                                       ypel$eff2020,
                                       ypel$eff2021)),
                    nsites = nrow(ypel),
                    nocc = 5)

yy <- ypel %>%
  dplyr::select(starts_with("y20")) %>%
  st_drop_geometry()

data.o <- list(nobs = yy)

Ninit <- apply(yy,1,sum)+1
inits.o <- list(a = rnorm(4,0,1), b = rnorm(2,0,1), N = Ninit)

Rmodelo <- nimbleModel(code = ipp.ni, constants = constants.o,
                      data = data.o, inits = inits.o)

# start the nimble process
#code
nmixtureModel <- nimbleCode({
  
  # priors
  for(i in 1:4){
    a[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    b[i] ~ dnorm(0,1)
  }
  
  # latent occu
  log(lambda[1:nsites]) <- a[1] + a[2] * prof[1:nsites] +
    a[3] * prof[1:nsites] * prof[1:nsites] +
    a[4] * dcol[1:nsites]
  
  # observation process
  logit(p[1:nsites,1:nocc]) <- b[1] + b[2] * seff[1:nsites, 1:nocc]
  
  for(i in 1:nsites){
    # likelihood
    N[i] ~ dpois(lambda[i])
    
    for(j in 1:nocc){
      nobs[i,j] ~ dbin(p[i,j],N[i])
    }
  }
})

# nimble process
Rmodelo$initializeInfo()
Rmodelo$calculate() # - 6515

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
resIPP <- runMCMC(Cmcmco, niter = 510000, nburnin = 50000, nchains = 2, samplesAsCodaMCMC = TRUE)

# check convergence
mcmcplots::denplot(resNmix)
coda::effectiveSize(resNmix)


# store results
res_b0ipp <- rbind(resIPP$chain1, resIPP$chain2) %>%
  as_tibble()%>%
  dplyr::select(starts_with("a[1]"))

res_profipp <- rbind(resIPP$chain1, resIPP$chain2) %>%
  as_tibble()%>%
  dplyr::select(starts_with("a[2]"))

res_qprofipp <- rbind(resIPP$chain1, resIPP$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[3]"))

res_dcolipp <- rbind(resIPP$chain1, resIPP$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[4]"))
