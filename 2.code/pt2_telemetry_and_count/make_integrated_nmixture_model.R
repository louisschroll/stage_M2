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
source("~/stage_M2/2.code/pt2_telemetry_and_count/functions_for_nmixture.R")

# Prepare data
species <- "sterne_caugek_R"

data_list <- list(pelmed = list(obs_data = pelmed_obs, effort_data = pelmed_eff), 
                  migralion = list(obs_data = migralion_obs, effort_data = migralion_eff))

grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed_obs)) 

selected_cov <- c("mean_CHL", "mean_SST", "dist_to_shore")

data_nmix <- prepare_data_for_int_nmix(data_list, grid, species, selected_cov)

# Nimble model
int.Nmixture.model <- nimbleCode({
  # priors
  for(i in 1:n.occ.cov){
    a[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    b1[i] ~ dnorm(0,1)
    b2[i] ~ dnorm(0,1)
  }
  
  # likelihood
  # state process
  for(i in 1:nsites_total){
    log(lambda[i]) <- sum(a[1:n.occ.cov] * XN[i,1:n.occ.cov])
    N[i] ~ dpois(lambda[i])
  }
  
    logit(p1[1:nsites1,1:nreplicates1]) <- b1[1] + b1[2] * transect_length1[1:nsites1, 1:nreplicates1]
    for(i in 1:nsites1){
      # observation process
      for(j in 1:nreplicates1){
        nobs1[i,j] ~ dbin(p1[i,j], N[site_id1[i]])
      }
    }
  
    logit(p2[1:nsites2,1:nreplicates2]) <- b2[1] + b2[2] * transect_length2[1:nsites2, 1:nreplicates2]
    for(i in 1:nsites2){
      # observation process
      for(j in 1:nreplicates2){
        nobs2[i,j] ~ dbin(p2[i,j], N[site_id2[i]])
      }
    }
})

# Prepare data
n.occ.cov <- length(selected_cov) + 1
nsites_total <- data_nmix$constants$site_id %>% unlist() %>% n_distinct()

constants <- list(XN = data_nmix$constants$XN,
                  site_id1 = data_nmix$constants$site_id$pelmed,
                  site_id2 = data_nmix$constants$site_id$migralion,
                  transect_length1 = data_nmix$constants$transect_length$pelmed,
                  transect_length2 = data_nmix$constants$transect_length$migralion,
                  nsites1 = data_nmix$constants$nsites[1],
                  nsites2 = data_nmix$constants$nsites[2],
                  nreplicates1 = data_nmix$constants$nreplicates[1],
                  nreplicates2 = data_nmix$constants$nreplicates[2],
                  n.occ.cov = data_nmix$constants$n.occ.cov,
                  nsites_total = nsites_total)



ff <- data_nmix$data$pelmed %>% mutate(cells_id = data_nmix$constants$site_id$pelmed) %>% 
  full_join(data_nmix$data$migralion %>% mutate(cells_id = data_nmix$constants$site_id$migralion), by = join_by(cells_id)) %>% 
  arrange(by = cells_id) %>% 
  select(-cells_id)
  
initial.values <- list(a = rnorm(n.occ.cov,0,1), 
                       b1 = rnorm(2,0,1),
                       b2 = rnorm(2,0,1),
                       N = apply(ff, 1, function(x){sum(x, na.rm = T)}) + 1)

parameters.to.save <- c("a", "b1", "lambda", "p1", "N")

# Nombre d'iterations, burn-in et nombre de chaine
n.iter <- 100000
n.burnin <- 10000
n.chains <- 2

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
                       constants = constants,
                       data = list(nobs1 = data_nmix$data$pelmed, nobs2 = data_nmix$data$migralion), 
                       inits = initial.values)

Rmodelo$initializeInfo()
Rmodelo$calculate()
Rmodelo$logProb_nobs2 %>% as_tibble() %>% print(n = 15)
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
mcmcplots::traplot(mcmc.output)
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

coeff_values <- cbind(res_b0ipp, res_profipp, res_qprofipp, res_dcolipp)


new_grid <- make_prediction(mcmc.output, grid, selected_cov)

load("~/stage_M2/1.data/contour_golfe_du_lion.rdata")
load("~/stage_M2/1.data/countLL.rdata")

plot_prediction <- function(new_grid, add_colonies = F, species_colony=NULL){
  mean_psi_plot <- ggplot() + 
    geom_sf(data = new_grid, aes(fill = mean_psi), color = NA) +
    scale_fill_distiller(palette = "Spectral",
                         guide = guide_colorbar(ticks = FALSE,
                                                barwidth = 9,
                                                barheight = 0.7)) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = "Intensity of space use") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  if (add_colonies){
    mean_psi_plot <- mean_psi_plot +
      geom_sf(data = countLL %>% filter(Species == species_colony) %>% st_crop(st_bbox(contour_golfe)),
              pch = 16, size = 2)
  }
    
  sd_psi_plot <- ggplot() + geom_sf(data = new_grid, aes(fill = sd_psi), color = NA) +
    scale_fill_viridis_c(option = "inferno", 
                         guide = guide_colorbar(ticks = FALSE,
                                                barwidth = 9,
                                                barheight = 0.7)) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = "SD intensity") +
    theme_bw() +
    theme(legend.position = "top", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(list(mean_psi_plot = mean_psi_plot, sd_psi_plot = sd_psi_plot))
}


plots <- plot_prediction(new_grid, add_colonies = T, species_colony = "Sterne caugek")
library(patchwork)
plots$mean_psi_plot + plots$sd_psi_plot










# Nimble model
int.Nmixture.model <- nimbleCode({
  
  # priors
  for(i in 1:n.occ.cov){
    a[i] ~ dnorm(0,1)
  }
  
  for(i in 1:2){
    b[i] ~ dnorm(0,1)
  }
  
  # likelihood
  for(ds in 1:ndatasets){
    logit(p[[ds]][1:nsites[ds],1:nreplicates[ds]]) <- b[1] + b[2] * transect_length[[ds]][1:nsites[ds], 1:nreplicates[ds]]
    for(i in 1:nsites[ds]){
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

