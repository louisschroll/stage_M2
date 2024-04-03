
setwd("~/stage_M2/1.data")

# load packages
library(tidyverse)
library(sf)

# # load the grid and the boundaries
# load("gdlhex.rdata") # the grid is loaded in the object sea.gdl2
# load("gdlmap.rdata") # the boundary is loaded in gdlmap
# 
# # Check that coordinate reference system is the same
# st_crs(sea.gdl2$x) == st_crs(gdlmap$geometry)

# # Compute distance to coast for each grid cell
# dist_to_coast <- unlist(map2(sea.gdl2$x, st_union(gdlmap), st_distance))

if(FALSE){
  ggplot() +
    geom_sf(data = sea.gdl2, aes(fill = dist_to_coast), lwd = 0.1) + 
    scale_fill_viridis_c(direction = -1 ) +
    geom_sf(data = gdlmap) 
}

# save(covariates_data, file = "covariates_data.rdata")

# ---- load data and filter CAUGEK ----
# Pelmed 2017 -> 2020
load("pelmed.rdata")
pelmed <- pelmed_obs %>% 
  filter(nom_fr == "Sterne caugek")

remove(pelmed_obs)

# Fermes pilotes 2011, 2013, 2016 - 2018
load("pilotebiotope.rdata")
ferme <- efglx %>% 
  filter(nomCite == c("Thalasseus sandvicensis (Latham, 1787)","Sterne caugek"))

remove(efglx)

# Migralion lot 4
load("prenup2022.rdata")
load("postnup2022.rdata")

migralion <- prenup22_obs %>%
  select(Espece, Effectif, geometry, Date_UTC) %>% 
  bind_rows(postnup2022_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC)) %>% 
  filter(Espece == "Sterne caugek")

remove(prenup22_obs, postnup2022_obs)

# PNM Golfe du Lion
load("megaobs.rdata")
pnm <- obs_oiseaux %>% 
  filter(espece == "Sterne Caugek")

# SAMM 2011/2012 - 2018-2019
load("birdSAMM.rdata")
samm <- birdsamm %>% 
  filter(nom_fr == "Sterne caugek")

#### Integrated model
library(spOccupancy)
library(coda)
library(stars)
library(tidyverse)
set.seed(102)
source("~/stage_M2/2.code/1.data_formating_functions.R")

#---- Load grid cells  ----

load("gdlmap.rdata")
load("~/stage_M2/1.data/covariates_data.rdata")
plot(covariates_data)

# name grid cells
grid <- covariates_data %>% 
  mutate(id = 1:nrow(covariates_data)) %>% 
  st_transform(st_crs(pelmed))


# --- prepare data & check outliers ---- 

# store which cells of `grid` are sampled, for each dataset
sitesocc <- extract_sampled_cells(do.pel = T, do.sam = T, do.pnm = T, do.mig = T, grid = grid)

sites.int <- list(pelmed = sitesocc$pelmed, 
                  samm = sitesocc$samm, 
                  pnm = sitesocc$pnm, 
                  migralion = sitesocc$migralion)

# Give the correct format to the data (as defined in spOccupancy)
pelmed_list = format_pelmed_data(pelmed_obs=pelmed, pelmed_eff=pelmed_eff, sitesocc=sitesocc)
samm_list = format_samm_data(samm_obs = samm, samm_eff = effsamm, sitesocc=sitesocc)
pnm_list = format_pnm_data(pnm_obs = pnm, pnm_eff = transect, sitesocc=sitesocc)
migralion_list = format_migralion_data(migralion_obs = migralion, 
                                       migralion_eff = postnup2022_eff, sitesocc=sitesocc)


y.int <- list(pelmed = pelmed_list$y, 
              samm = samm_list$y,
              pnm = pnm_list$y,
              migralion = migralion_list$y)
str(y.int)

det.covs.int <- list(pelmed = pelmed_list$det.covs, 
                     samm = samm_list$det.covs,
                     pnm = pnm_list$det.covs,
                     migralion = migralion_list$det.covs)

occ.covs.int <- sitesocc$grid %>% 
  as_tibble() %>% 
  select("bathymetry", "dist_to_shore", "slope", "winter_SST",
         "spring_SST", "summer_SST", "autumn_SST", "concavity",
         "mean_CHL")

str(occ.covs.int)

data.int <- list(y = y.int, 
                 occ.covs = occ.covs.int, 
                 det.covs = det.covs.int, 
                 sites = sites.int)
str(data.int)

occ.formula.int <- as.formula(paste("~" , 
                                    paste(#"bathymetry", 
                                          "dist_to_shore", 
                                          "slope", 
                                          "winter_SST",
                                          "spring_SST", 
                                          "summer_SST", 
                                          "autumn_SST", 
                                          "concavity",
                                          "mean_CHL", sep = " + ")))

det.formula.int = list(pelmed = ~ scale(transect_length), 
                       samm = ~ scale(transect_length),
                       pnm = ~ scale(transect_length),
                       migralion = ~ scale(transect_length))


# Total number of sites
J <- nrow(data.int$occ.covs)
inits.list <- list(alpha = list(0, 0, 0, 0),
                   beta = 0, 
                   z = rep(1, J))

prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72)))

n.samples <- 8000
n.burn <- 3000
n.thin <- 25

# Approx. run time: < 15 sec
out.int <- intPGOcc(occ.formula = occ.formula.int,
                    det.formula = det.formula.int, 
                    data = data.int,
                    inits = inits.list,
                    n.samples = n.samples, 
                    priors = prior.list, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    n.report = 2000, 
                    n.burn = n.burn, 
                    n.thin = n.thin, 
                    n.chains = 3) 

summary(out.int)
plot(out.int$beta.samples, density = T)

ppc.int.out <- ppcOcc(out.int, 'freeman-tukey', group = 2)
summary(ppc.int.out)

waicOcc(out.int)

out.int.k.fold <- intPGOcc(occ.formula = occ.formula.int,
                           det.formula = det.formula.int, 
                           data = data.int,
                           inits = inits.list,
                           n.samples = n.samples, 
                           priors = prior.list, 
                           n.omp.threads = 1, 
                           verbose = FALSE, 
                           n.report = 2000, 
                           n.burn = n.burn, 
                           n.thin = n.thin, 
                           n.chains = 1,
                           k.fold = 4) 
out.int.k.fold.small <- intPGOcc(occ.formula = ~ 1, 
                                 det.formula = list(hbef = ~ 1, neon = ~ 1), 
                                 data = data.int,
                                 inits = inits.list,
                                 n.samples = n.samples, 
                                 priors = prior.list, 
                                 n.omp.threads = 1, 
                                 verbose = FALSE, 
                                 n.report = 2000, 
                                 n.burn = n.burn, 
                                 n.thin = n.thin, 
                                 n.chains = 1,
                                 k.fold = 4) 
# Summarize the CV results
out.int.k.fold$k.fold.deviance
out.int.k.fold.small$k.fold.deviance


out.int.k.fold.hbef <- intPGOcc(occ.formula = occ.formula.int,
                                det.formula = det.formula.int, 
                                data = data.int,
                                inits = inits.list,
                                n.samples = n.samples, 
                                priors = prior.list, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                n.report = 2000, 
                                n.burn = n.burn, 
                                n.thin = n.thin, 
                                k.fold = 4, 
                                k.fold.data = 1) 

# Look at CV results again 
# Single data source model
out.k.fold$k.fold.deviance
# Integrated model
out.int.k.fold.hbef$k.fold.deviance


# Make sure to standardize using mean and sd from fitted model
# depth.pred <- (grid$depth.sc - mean(data.int$occ.covs[, 1])) / sd(data.int$occ.covs[, 1])

# Make sure to standardize using mean and sd from fitted model
depth.pred <- (grid$bathymetry)
dist.pred <- (grid$dist_to_shore)
slope.pred <- (grid$slope)

grid_pred <- grid %>% 
  as_tibble() %>% 
  #mutate(bathymetry_quadra = bathymetry^2) %>% 
  select("dist_to_shore", "slope", "winter_SST",
         "spring_SST", "summer_SST", "autumn_SST", "concavity",
         "mean_CHL")

X.0 <- cbind(1, grid_pred)
out.int.pred <- predict(out.int, X.0)
mean.psi = apply(out.int.pred$psi.0.samples, 2, mean)
sd.psi = apply(out.int.pred$psi.0.samples, 2, sd)

psi <- ggplot() + 
  geom_sf(data = grid, aes(fill = mean.psi), lwd = 0.1) +
  scale_fill_viridis_c() + 
  labs(title = 'Occupancy') +
  theme_bw()

sdpsi <- ggplot() + 
  geom_sf(data = grid, aes(fill = sd.psi), lwd = 0.1) +
  scale_fill_viridis_c(option = "B") + 
  labs(title = 'Occupancy') +
  theme_bw()

library(patchwork)
psi + sdpsi


# ----- With spatial autocorrelation spIntPGOcc() -----
coords = st_coordinates(st_centroid(sitesocc$grid))

sites.int <- list(pelmed = sitesocc$pelmed, 
                  samm = sitesocc$samm, 
                  pnm = sitesocc$pnm, 
                  migralion = sitesocc$migralion)


data.int <- list(y = y.int, 
                 occ.covs = occ.covs.int, 
                 det.covs = det.covs.int, 
                 sites = sites.int,
                 coords = coords)

dist.int <- dist(data.int$coords)
min.dist <- min(dist.int)
max.dist <- max(dist.int)
J <- nrow(data.int$det.covs)

inits.list <- list(alpha = list(0, 0, 0, 0),
                   beta = 0, 
                   z = rep(1, J),
                   sigma.sq = 2,
                   phi = 3 / mean(dist.int), 
                   w = rep(0, J))

prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72)),
                   sigma.sq.ig = c(2, 2), 
                   phi.unif = c(3 / max.dist, 1))

n.batch = 4
batch.length = 1000
n.burn = 1000
tuning <- list(phi = .2)

out.sp.int <- spIntPGOcc(occ.formula = occ.formula.int, 
                    det.formula = det.formula.int, 
                    data = data.int,  
                    inits = inits.list, 
                    n.batch = n.batch, 
                    batch.length = batch.length, 
                    accept.rate = 0.43, 
                    priors = prior.list, 
                    cov.model = "exponential", 
                    tuning = tuning, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    NNGP = T, 
                    n.report = 1000, 
                    n.burn = 50, 
                    n.thin = 2,
                    n.chains = 2)


summary(out.sp.int)
waicOcc(out.sp.int)

ppc.sp.int.out <- ppcOcc(out.sp.int, 'freeman-tukey', group = 1)
summary(ppc.sp.int.out)

# Make sure to standardize using mean and sd from fitted model

depth.pred <- grid$depth.sc[,1]
dist.pred = scale(grid$dist_to_coast)
X.0 <- cbind(1, depth.pred, depth.pred^2, dist.pred)
coords.0 <- st_coordinates(st_centroid(grid))
out.int.pred <- predict(out.sp.int, X.0, coords.0)

mean.psi = apply(out.int.pred$psi.0.samples, 2, mean)
sd.psi = apply(out.int.pred$psi.0.samples, 2, sd)

psi2 <- ggplot() + 
  geom_sf(data = grid, aes(fill = mean.psi), lwd = 0.1) +
  scale_fill_viridis_c() + 
  labs(title = 'Occupancy') +
  theme_bw()

sdpsi2 <- ggplot() + 
  geom_sf(data = grid, aes(fill = sd.psi), lwd = 0.1) +
  scale_fill_viridis_c(option = "B") + 
  labs(title = 'Occupancy') +
  theme_bw()

library(patchwork)
psi2 + sdpsi2

psi + psi2

ggplot(data = world) +
  geom_sf()

