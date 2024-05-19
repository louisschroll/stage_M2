

run_integrated_Nmix_RSF <- function(data_nmix, 
                                    data_rsf, 
                                    nmix_model = "Poisson",
                                    n.iter = 100000, 
                                    n.burnin = 10000,
                                    n.chains = 3) {
  
  # Combine data in list (data, inits and constants)
  # INITS
  int.inits <- c(data_rsf$inits, # keep beta and sd inits from RSF
                 data_nmix$inits[2:length(data_nmix$inits)], # add inits for alpha and N
                 list(beta0_nmix = rnorm(1,0,1))) # add inits for N-mix intercept
  # CONSTANTS
  ## Rename occurence covariates to have different names
  names(data_rsf$constants)[names(data_rsf$constants) == "XN"] <- "XN_rsf"
  names(data_nmix$constants)[names(data_nmix$constants)=="XN"] <- "XN_nmix"
  
  data_nmix$constants[which(names(data_nmix$constants) %in% c("n.occ.cov"))] <- NULL
  int.constants <- c(data_rsf$constants, data_nmix$constants)
  # DATA
  int.data <- c(data_rsf$data, data_nmix$data)
  
  # Nimble code
  if (nmix_model == "Poisson") {
    int.Nmix.RSF.code <- nimbleCode({
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
        for (nd in 1:ndatasets){
          alpha[i, nd] ~ dnorm(0,1)
        }
      }
      
      # likelihood
      ## RSF
      for(t in 1:npoints){
        logit(omega[t]) <- sum(beta_ind[idind[t], 1:n.occ.cov] * XN_rsf[t, 1:n.occ.cov])
        kase[t] ~ dbinom(omega[t], w[t])
      }
      ## N mixture
      ### State process
      for(i in 1:nsites_total){
        log(lambda[i]) <- beta0_nmix + sum(beta_pop[2:n.occ.cov] * XN_nmix[i,2:n.occ.cov])
        N[i] ~ dpois(lambda[i])
      }
      ### Observation process
      for (i in 1:nsampled_points){
        logit(p[i]) <- alpha[1, dataset_nb[i]] + alpha[2, dataset_nb[i]] * transect_length[i]
        nobs[i] ~ dbin(p[i], N[site_id[i]])
      }
    })
    
  } else if (nmix_model == "NB"){
    int.inits$kappa <- 0.5
    int.Nmix.RSF.code <- nimbleCode({
      # Priors
      kappa ~ dunif(min = 0.01, max = 100)
      
      for(i in 1:2){
        for (nd in 1:ndatasets){
          alpha[i, nd] ~ dnorm(0,1)
        }
      }
      
      ## intercept
      beta0_nmix ~ dnorm(0,1) # N mixture intercept
      beta_pop[1] ~ dnorm(0, 1) # RSF intercept
      for(j in 1:nindividual){
        # add individual heterogeneity
        beta_ind[j, 1] ~ dnorm(beta_pop[1], 100)
      }
      
      for(i in 2:n.occ.cov){
        beta_pop[i] ~ dnorm(0, 1)
        sd_pop[i-1] ~ dunif(0, 1e2)
        
        for(j in 1:nindividual){
          # add individual heterogeneity
          beta_ind[j, i] ~ dnorm(beta_pop[i], sd = sd_pop[i-1])
        }
      }
      
      # Likelihood
      ## RSF
      for(t in 1:npoints){
        logit(omega[t]) <- sum(beta_ind[idind[t], 1:n.occ.cov] * XN_rsf[t, 1:n.occ.cov])
        kase[t] ~ dbinom(omega[t], w[t])
      }
      ## N mixture
      ### State process
      for(i in 1:nsites_total){
        log(lambda[i]) <- beta0_nmix + sum(beta_pop[2:n.occ.cov] * XN_nmix[i, 2:n.occ.cov])
        succprob <- kappa / (kappa + lambda[i])
        N[i] ~ dnegbin(prob = succprob, size = kappa)
      }
      
      ### Observation process
      for (i in 1:nsampled_points){
        logit(p[i]) <- alpha[1, dataset_nb[i]] + alpha[2, dataset_nb[i]] * transect_length[i]
        nobs[i] ~ dbin(p[i], N[site_id[i]])
      }
    })
  }
  
  # Nimble pre run
  int.model <- nimbleModel(code = int.Nmix.RSF.code, 
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
  return(samplesint)
}


