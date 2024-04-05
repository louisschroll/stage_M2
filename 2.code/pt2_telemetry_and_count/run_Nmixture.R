#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/run_Nmixture.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-04
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------


run_Nmixture <- function(data_nmix, n.iter = 100000, n.burnin = 10000, n.chains = 2){
  ndataset <- length(data_nmix$data)
  
  if(ndataset == 1){
    # Nimble model for 1 data source
  } else if (ndataset == 2){
    # Nimble model for 2 data sources
    Nmixture.model <- nimbleCode({
      # Priors
      for(i in 1:n.occ.cov){
        beta[i] ~ dnorm(0,1)
      }
      
      for(i in 1:2){
        alpha1[i] ~ dnorm(0,1)
        alpha2[i] ~ dnorm(0,1)
      }
      
      # Likelihood
      ## state process
      for(i in 1:nsites_total){
        log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
        N[i] ~ dpois(lambda[i])
      }
      ## observation process
      ### data source 1
      logit(p1[1:nsites1,1:nreplicates1]) <- alpha1[1] + alpha1[2] * transect_length1[1:nsites1, 1:nreplicates1]
      for(i in 1:nsites1){
        for(j in 1:nreplicates1){
          nobs1[i,j] ~ dbin(p1[i,j], N[site_id1[i]])
        }
      }
      ### data source 2
      logit(p2[1:nsites2,1:nreplicates2]) <- alpha2[1] + alpha2[2] * transect_length2[1:nsites2, 1:nreplicates2]
      for(i in 1:nsites2){
        for(j in 1:nreplicates2){
          nobs2[i,j] ~ dbin(p2[i,j], N[site_id2[i]])
        }
      }
    })
  } else if (ndataset == 3){
    # Nimble model for 3 data sources
    Nmixture.model <- nimbleCode({
      # Priors
      for(i in 1:n.occ.cov){
        beta[i] ~ dnorm(0,1)
      }
      
      for(i in 1:2){
        alpha1[i] ~ dnorm(0,1)
        alpha2[i] ~ dnorm(0,1)
        alpha3[i] ~ dnorm(0,1)
      }
      
      # Likelihood
      ## state process
      for(i in 1:nsites_total){
        log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
        N[i] ~ dpois(lambda[i])
      }
      ## Observation process
      ### data source 1
      logit(p1[1:nsites1,1:nreplicates1]) <- alpha1[1] + alpha1[2] * transect_length1[1:nsites1, 1:nreplicates1]
      for(i in 1:nsites1){
        for(j in 1:nreplicates1){
          nobs1[i,j] ~ dbin(p1[i,j], N[site_id1[i]])
        }
      }
      ### data source 2
      logit(p2[1:nsites2,1:nreplicates2]) <- alpha2[1] + alpha2[2] * transect_length2[1:nsites2, 1:nreplicates2]
      for(i in 1:nsites2){
        for(j in 1:nreplicates2){
          nobs2[i,j] ~ dbin(p2[i,j], N[site_id2[i]])
        }
      }
      ### data source 3
      logit(p3[1:nsites3,1:nreplicates3]) <- alpha3[1] + alpha3[2] * transect_length3[1:nsites3, 1:nreplicates3]
      for(i in 1:nsites3){
        for(j in 1:nreplicates3){
          nobs3[i,j] ~ dbin(p3[i,j], N[site_id3[i]])
        }
      }
    })
  } else if (ndataset == 4){
    # Nimble model for 4 data sources
    Nmixture.model <- nimbleCode({
      # Priors
      for(i in 1:n.occ.cov){
        beta[i] ~ dnorm(0,1)
      }
      
      for(i in 1:2){
        alpha1[i] ~ dnorm(0,1)
        alpha2[i] ~ dnorm(0,1)
        alpha3[i] ~ dnorm(0,1)
        alpha4[i] ~ dnorm(0,1)
      }
      
      # Likelihood
      ## state process
      for(i in 1:nsites_total){
        log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
        N[i] ~ dpois(lambda[i])
      }
      ## Observation process
      ### data source 1
      logit(p1[1:nsites1,1:nreplicates1]) <- alpha1[1] + alpha1[2] * transect_length1[1:nsites1, 1:nreplicates1]
      for(i in 1:nsites1){
        for(j in 1:nreplicates1){
          nobs1[i,j] ~ dbin(p1[i,j], N[site_id1[i]])
        }
      }
      ### data source 2
      logit(p2[1:nsites2,1:nreplicates2]) <- alpha2[1] + alpha2[2] * transect_length2[1:nsites2, 1:nreplicates2]
      for(i in 1:nsites2){
        for(j in 1:nreplicates2){
          nobs2[i,j] ~ dbin(p2[i,j], N[site_id2[i]])
        }
      }
      ### data source 3
      logit(p3[1:nsites3,1:nreplicates3]) <- alpha3[1] + alpha3[2] * transect_length3[1:nsites3, 1:nreplicates3]
      for(i in 1:nsites3){
        for(j in 1:nreplicates3){
          nobs3[i,j] ~ dbin(p3[i,j], N[site_id3[i]])
        }
      }
      ### data source 4
      logit(p4[1:nsites4,1:nreplicates4]) <- alpha4[1] + alpha4[2] * transect_length4[1:nsites4, 1:nreplicates4]
      for(i in 1:nsites4){
        for(j in 1:nreplicates4){
          nobs4[i,j] ~ dbin(p4[i,j], N[site_id4[i]])
        }
      }
    })
  }
  
  parameters.to.save <- c("beta", paste0("alpha", 1:ndataset), "lambda", "N", paste0("p", 1:ndataset))
  
  # In one step
  # mcmc.output <- nimbleMCMC(code = int.Nmixture.model,
  #                           data = my.data,
  #                           constants = my.constants,
  #                           inits = initial.values,
  #                           monitors = parameters.to.save,
  #                           niter = n.iter,
  #                           nburnin = n.burnin,
  #                           nchains = n.chains, 
  #                           samplesAsCodaMCMC = TRUE)
  
  # In several step
  Rmodelo <- nimbleModel(code = int.Nmixture.model, 
                         constants = data_nmix$constants,
                         data = data_nmix$data, 
                         inits = data_nmix$inits)
  
  Rmodelo$initializeInfo()
  Rmodelo$calculate()
  # configure model
  confo <- configureMCMC(Rmodelo)
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
  return(mcmc.output)
}