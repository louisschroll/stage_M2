#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  run_RSF.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-03
#'
#' Script description:
#' This function run a Resource Selection Function using nimble and return the 
#' resulting mcmc object.
#' @df_RSF: a data frame, which must contain at least indivual_id (an integer),
#' the case (0 or 1) and the covariates data. 
#' @selected_cov: a character vector indicating the covariate the use for the 
#' modelling
#' @output: a mcmc or mcmc.list object with the values taken by the parameters
#' -------------------------------------------------------------------------------


run_RSF <- function(df_RSF, selected_cov, niter = 11000, nburnin = 1000, thin = 1, nchains = 2){
  # Nimble model
  rsf.nimble <- nimbleCode({
    # Priors
    for(i in 1:n.occ.cov){
      beta_pop[i] ~ dnorm(0,1)
      sd_pop[i] ~ dunif(0,1e2)
      
      for(j in 1:nindividual){
        # add individual heterogeneity
        beta_ind[j, i] ~ dnorm(beta_pop[i], sd = sd_pop[i])
      }
    }
    # Likelihood
    for(t in 1:npoints){
      logit(omega[t]) <- sum(beta_ind[idind[t], 1:n.occ.cov] * XN[t, 1:n.occ.cov])
      kase[t] ~ dbinom(omega[t], w[t])
    }
  })
  
  # Add weight (this is done to have more pseudo-absence while saving computational ressources)
  weight <- df_RSF$case
  weight[weight==0] <- 1000
  
  # Define constants
  n.occ.cov <- length(selected_cov)
  nindividual = n_distinct(df_RSF$individual_id)
  
  constants <-  list(npoints = nrow(df_RSF),
                        idind = df_RSF$individual_id,
                        XN = df_RSF %>% select(all_of(selected_cov)) %>% as.matrix(),
                        nindividual = nindividual,
                        n.occ.cov = n.occ.cov,
                        w = weight)
  
  # Define data
  data <- list(kase = df_RSF$case)
  
  # Inits
  inits <-  list(beta_pop = rep(0, n.occ.cov),
                 sd_pop = runif(n.occ.cov, 0, 5),
                 beta_ind = matrix(0, ncol = n.occ.cov, nrow = nindividual)
  )
  
  # Nimble pre run
  rsf.model <- nimbleModel(code= rsf.nimble, constants = constants, data = data, inits = inits)
  rsf.model$initializeInfo()
  rsf.model$calculate()
  
  # configure model
  model_configuration <- configureMCMC(rsf.model)
  model_configuration$printMonitors()
  
  # Build and compile MCMC
  Rmcmc <- buildMCMC(model_configuration)
  Cmodel <- compileNimble(rsf.model)
  Cmcmc <- compileNimble(Rmcmc, project = Cmodel)
  
  # Run
  samplesRSF <- runMCMC(Cmcmc, 
                        niter = niter, 
                        nburnin = nburnin, 
                        thin = thin, 
                        nchains = nchains, 
                        samplesAsCodaMCMC = TRUE)
  return(samplesRSF)
}

