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


run_Nmixture <- function(data_nmix, n.iter = 100000, n.burnin = 10000, n.chains = 3, compute_pvalues = F, zero_inflated = F, overdispersion = F){
  # Nimble code
  if (compute_pvalues){
    data_nmix$constants$npoints_dataset <- c(0, table(data_nmix$constants$dataset_nb) %>% unname())
    # Version with bayesian p-values computation
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
        
        # Compute discrepancy for real and simulated data
        ## Expected count at this site and this survey
        exp_count[i] <- N[site_id[i]] * p[i] 
        
        ## Discrepancy for the real data
        ## (small value added to denominator to avoid potential divide by zero)
        E[i] <- pow((nobs[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
        
        ## Simulate new count from model
        nobs.rep[i] ~ dbin(p[i], N[site_id[i]])
        
        ## Discrepancy for the simulated data
        E.rep[i] <- pow((nobs.rep[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
      }
      
      # chi-squared test statistics
      for (nd in 1:ndatasets){
        fit[nd] <- sum(E[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
        fit.rep[nd] <- sum(E.rep[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
      }      
    })
    parameters.to.save <- c("beta", "alpha", "fit", "fit.rep")
  } else if (zero_inflated) {
    data_nmix$constants$npoints_dataset <- c(0, table(data_nmix$constants$dataset_nb) %>% unname())
    data_nmix$inits$omega <- runif(1, 0.5, 0.7)
    data_nmix$inits$z <- ifelse(data_nmix$inits$N == 1, 0, 1)
    # Version with bayesian p-values computation
    int.Nmixture.code <- nimbleCode({
      # Priors
      omega ~ dbeta(1, 1)
      
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
        z[i] ~ dbern(omega)
        
        log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
        N[i] ~ dpois(lambda[i] * z[i])
      }
      
      # Observation process
      for (i in 1:nsampled_points){
        logit(p[i]) <- alpha[1, dataset_nb[i]] + alpha[2, dataset_nb[i]] * transect_length[i]
        nobs[i] ~ dbin(p[i], N[site_id[i]])
        
        # Compute discrepancy for real and simulated data
        ## Expected count at this site and this survey
        exp_count[i] <- N[site_id[i]] * p[i] 
        
        ## Discrepancy for the real data
        ## (small value added to denominator to avoid potential divide by zero)
        E[i] <- pow((nobs[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
        
        ## Simulate new count from model
        nobs.rep[i] ~ dbin(p[i], N[site_id[i]])
        
        ## Discrepancy for the simulated data
        E.rep[i] <- pow((nobs.rep[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
      }
      
      # chi-squared test statistics
      for (nd in 1:ndatasets){
        fit[nd] <- sum(E[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
        fit.rep[nd] <- sum(E.rep[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
      }      
    })
    parameters.to.save <- c("beta", "alpha", "fit", "fit.rep")
  } else if (overdispersion) {
    data_nmix$constants$npoints_dataset <- c(0, table(data_nmix$constants$dataset_nb) %>% unname())
    data_nmix$inits$omega <- runif(1, 0.5, 0.7)
    data_nmix$inits$sd.p = runif(1, 0, 0.25)
    data_nmix$inits$z <- ifelse(data_nmix$inits$N == 1, 0, 1)
    # Version with bayesian p-values computation
    int.Nmixture.code <- nimbleCode({
      # Priors
      omega ~ dbeta(1, 1)
      # tau.p <- pow(sd.p, -2)
      sd.p ~ dunif(0, 3)
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
        z[i] ~ dbern(omega)
        
        log(lambda[i]) <- sum(beta[1:n.occ.cov] * XN[i,1:n.occ.cov])
        N[i] ~ dpois(lambda[i] * z[i])
      }
      
      # Observation process
      for (i in 1:nsampled_points){
        logit(p[i]) <- lp[i]
        mu.lp[i] <- alpha[1, dataset_nb[i]] + alpha[2, dataset_nb[i]] * transect_length[i]
        lp[i] ~ dnorm(mu.lp[i], sd = sd.p)
        
        nobs[i] ~ dbin(p[i], N[site_id[i]])
        
        # Compute discrepancy for real and simulated data
        ## Expected count at this site and this survey
        exp_count[i] <- N[site_id[i]] * p[i] 
        
        ## Discrepancy for the real data
        ## (small value added to denominator to avoid potential divide by zero)
        E[i] <- pow((nobs[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
        
        ## Simulate new count from model
        nobs.rep[i] ~ dbin(p[i], N[site_id[i]])
        
        ## Discrepancy for the simulated data
        E.rep[i] <- pow((nobs.rep[i] - exp_count[i]), 2) / (exp_count[i] + 0.005)
      }
      
      # chi-squared test statistics
      for (nd in 1:ndatasets){
        fit[nd] <- sum(E[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
        fit.rep[nd] <- sum(E.rep[(1 + npoints_dataset[nd]):(npoints_dataset[nd+1] + npoints_dataset[nd])])
      }      
    })
    parameters.to.save <- c("beta", "alpha", "fit", "fit.rep")
  } else {
    # Version without p-values computation
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
    
    parameters.to.save <- c("beta", "alpha")
  }
  
  # Define model
  int.Nmixture.model <- nimbleModel(code = int.Nmixture.code, 
                                    constants = data_nmix$constants,
                                    data = data_nmix$data, 
                                    inits = data_nmix$inits)
  
  int.Nmixture.model$initializeInfo()
  int.Nmixture.model$calculate()
  # Configure model
  confo <- configureMCMC(int.Nmixture.model, monitors = parameters.to.save)
  
  #confo$removeSampler(c("a[1]", "a[2]", "a[3]", "b[1]", "b[2]"))
  #confo$addSampler("a[1]", "a[2]", type = "RW_block")
  confo$removeSampler(paste0("beta[", 1:data_nmix$constants$n.occ.cov, "]"))
  confo$addSampler(paste0("beta[", 1:data_nmix$constants$n.occ.cov, "]"), 
                   type = "RW_block")
  confo$printSamplers()
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
  return(mcmc.output)
}
