# Integrated occupancy model for caugek terns

library(spOccupancy)
library(tidyverse)
library(sf)
library(lubridate)
library(here)


# load functions
source("~/stage_M2/2.code/code_valentin/fun_data.R")

# ---- load data and filter CAUGEK ----
# Pelmed 2017 -> 2020
load("~/stage_M2/1.data/pelmed.rdata")

pelmed <- pelmed_obs %>% 
  filter(nom_fr == "Sterne caugek")

remove(pelmed_obs)

# Fermes pilotes 2011, 2013, 2016 - 2018
load("~/stage_M2/1.data/pilotebiotope.rdata")

ferme <- efglx %>% 
  filter(nomCite == "Thalasseus sandvicensis (Latham, 1787)")

# efglx %>% 
#   filter(nomCite == "Sterne caugek")

remove(efglx)

# Migralion lot 4
load("~/stage_M2/1.data/prenup2022.rdata")
load("~/stage_M2/1.data/postnup2022.rdata")
imp <- prenup22_obs
prenup <- imp

migralion <- imp %>%
  select(Espece, Effectif, geometry, Date_UTC) %>% 
  bind_rows(postnup2022_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC)) %>% 
  filter(Espece == "Sterne caugek")


# birdsamm %>% group_by(nom_fr) %>% count(sort = TRUE) 

remove(imp, postnup2022_obs)

# PNM Golfe du Lion
load("~/stage_M2/1.data/megaobs.rdata")

pnm <- obs_oiseaux %>% 
  filter(espece == "Sterne Caugek")

ggplot() + geom_sf(data = gdlmap)+
  geom_sf(data = sea.gdl2) + 
  geom_sf(data = pnm)

# SAMM 2011/2012 - 2018-2019
load("birdSAMM.rdata")

samm <- birdsamm %>% 
  filter(nom_fr == "Sterne caugek")


#---- Load grid cells  ----

load("~/stage_M2/1.data/gdlmap.rdata")
load("~/stage_M2/1.data/gdlhex.rdata")

# change CRS
effsamm <- st_transform(effsamm, crs = st_crs(sea.gdl2))
postnup2022_eff <- st_transform(postnup2022_eff, crs = st_crs(sea.gdl2))

# name grid cells
grid <- sea.gdl2 %>% 
  mutate(id = 1:nrow(sea.gdl2))

# ---- frame covariates ----
# dates of sampling
samp_dates(do.pel = T, do.sam = T, do.pnm = T, do.mig = T, do.count = F)

# Go to ("Data/Fond de carte/env_cov.R) and extract relevant covariates

# --- prepare data & check outliers ---- 

# store which cells of `grid` are sampled, for each dataset
sst <- sites.l(do.pel = T, do.sam = T, do.pnm = T, do.mig = T, grid = grid)
gridocc <- sst$grid


sitesocc <- sst
sitesocc$grid <- NULL
COV <- gridocc$depth.sc
xoords <- st_coordinates(st_centroid(gridocc))
ydt <- yndet(do.pel = T, do.sam = T,
             do.pnm = T, do.mig = T)

if(FALSE){
  ## test for incoherences btw yy and detcov
  # dimensions should be identical
  
  dy <- ydt$yy[[3]]
  dp <- ydt$pp[[3]]$det.cov.3
  dim(dy) == dim(dp)
  
  # diff in na's
  diff.na <-  length(which(is.na(dy)==T)) - length(which(is.na(dp)==T))
  diff.na
  # quel element est différent ? 
  if( diff.na < 0){
    id <- setdiff(which(is.na(dp)) ,which(is.na(dy)))
    dp[id] <-  0
  }
  if( diff.na > 0){
    id <- setdiff(which(is.na(dy)) ,which(is.na(dp)))
    dy[id] <-  0
  }
  
  
  # quel ligne est differente ? 
  dyt <- dy
  dpt <- dp
  dyt[!is.na(dy)] <- 9  
  dyt[is.na(dy)] <- 0  
  
  dpt[!is.na(dp)] <- 9  
  dpt[is.na(dp)] <- 0  
  
  which(rowSums(dyt) != rowSums(dpt))
  
  dy[100,]
  dp[100,]
  
  dim(which(is.na(ydt$yy[[1]]), arr.ind = T))
  dim(which(is.na(ydt$pp$det.cov.1[[1]]), arr.ind = T))
  
  # dy2 <- ydt$yy[[2]]
  # dy2[is.na(dy2)] <- 0
  # dp2 <- ydt$pp[[2]]
  # dp2[is.na(dp2)] <- 0
  # 
  
  which(is.na(ydt$pp[[1]]), arr.ind = T)
  dy1 <- ydt$yy[[1]]
  dp1 <- ydt$pp[[1]]
  class(dp1$det.cov.1)
  dim(which(is.na(dy1), arr.ind = T))
  dim(which(is.na(dp1$det.cov.1), arr.ind = T))
  
  pa[(!(unique(which(is.na(pa),arr.ind = T)[,1]) %in%
          unique(which(is.na(a),arr.ind = T)[,1]))==T),]
  a[(!(unique(which(is.na(pa),arr.ind = T)[,1]) %in%
         unique(which(is.na(a),arr.ind = T)[,1]))==T),]
  
  which((which(is.na(dy1)) %in% which(is.na(dp1))) ==T)
  # ydt$pp[[1]][442] <- 0
  
  dim(which(is.na(ydt$yy[[2]]), arr.ind = T))
  dim(which(is.na(ydt$pp[[2]]), arr.ind = T)) 
  which(is.na(ydt$yy[[1]])) == which(is.na(ydt$pp[[1]]))
  
  dim(which(is.na(ydt$yy[[3]]), arr.ind = T))
  dim(which(is.na(ydt$pp[[3]]), arr.ind = T))
  
  dim(which(is.na(ydt$yy[[4]]), arr.ind = T))
  dim(which(is.na(ydt$pp[[4]]), arr.ind = T)) 
}

# data list 

str(ydt)

data.list <- list(y =  ydt$yy,
                  occ.covs = COV,
                  det.covs = ydt$pp, 
                  sites = sitesocc)
str(data.list)
# Initial values
inits.list <- list(alpha = list(0, 0, 0, 0), 
                   beta = 0, 
                   z = rep(1, nrow(COV)))
# Priors
prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0, 0, 0), 
                                       var = list(2.72, 2.72, 2.72, 2.72 )))#) 2.72, 2.72)))
n.samples <- 8000

out.int <- intPGOcc(occ.formula = ~ COV + I(COV^2), 
                    det.formula = list(f.1 = ~ det.cov.1, 
                                       f.2 = ~ det.cov.2, 
                                       f.3 = ~ det.cov.3,
                                       f.4 = ~ det.cov.4), 
                    #f.4 = ~ pd), 
                    data = data.list,
                    inits = inits.list,
                    n.samples = n.samples, 
                    priors = prior.list, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    n.report = 2000, 
                    n.burn = 3000, 
                    n.thin = 25, 
                    n.chains = 3)

summary(out.int)
plot(out.int$beta.samples, density = FALSE)

waicOcc(out.int)

# posterior predictive check
ppc.int.out <- ppcOcc(out.int, 'freeman-tukey', group = 1)
summary(ppc.int.out)
# A Bayesian p-value that hovers around 0.5 indicates adequate model fit,
# while values less than 0.1 or greater than 0.9 suggest our model does not fit the data well
#   ppc.df <- data.frame(fit = ppc.int.out$fit.y, 
#                        fit.rep = ppc.int.out$fit.y.rep, 
#                        color = 'lightskyblue1')
#   ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
#   plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
#        ylab = 'Fit', xlab = 'True')
#   lines(ppc.df$fit, ppc.df$fit, col = 'black')


# Make sure to standardize using mean and sd from fitted model
batpred <- grid$depth.sc[,1]
X.0 <- cbind(1, batpred, batpred^2)
out.int.pred <- predict(out.int, X.0)

# Producing an SDM for HBEF alone (posterior mean)
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

# ---- test spIntPGOcc()----

spdata <- list(y = ydt$yy, 
               occ.covs = COV,
               det.covs = ydt$pp, 
               sites = sitesocc,
               coords = xoords)

# Initial values
spinits <- list(alpha = list(0, 0, 0, 0), 
                beta = 0, 
                z = rep(1, nrow(COV)),
                phi = 3/ .5,
                sigma.sq = 2,
                w = rep(0, nrow(COV)))
# Priors
spprior <- list(beta.normal = list(mean = 0, var = 2.72), 
                alpha.normal = list(mean = list(0, 0, 0, 0), 
                                    var = list(2.72, 2.72, 2.72, 2.72)),# 2.72)),
                phi.unif = c(3/1, 3/.1), 
                sigma.sq.ig = c(2, 2))

# Tuning
tuning.list <- list(phi = 0.3) 

# Number of batches
n.batch <- 3
# Batch length
batch.length <- 1000

n.samples <- 5000
outsp <- spIntPGOcc(occ.formula = ~ COV, 
                    det.formula = list(f.1 = ~ det.cov.1,
                                       f.2 = ~ det.cov.2,
                                       f.3 = ~ det.cov.3,
                                       f.4 = ~ det.cov.4), 
                    data = spdata,  
                    inits = spinits, 
                    n.batch = n.batch, 
                    batch.length = batch.length, 
                    accept.rate = 0.43, 
                    priors = spprior, 
                    cov.model = "exponential", 
                    tuning = tuning.list, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    NNGP = T, 
                    n.report = 10, 
                    n.burn = 50, 
                    n.thin = 1,
                    n.chains = 2)


summary(outsp)
waicOcc(outsp)
# Make sure to standardize using mean and sd from fitted model

batpred <- grid$depth.sc[,1]
xoord.0 <- st_coordinates(st_centroid(grid))
X.0 <- cbind(1, batpred)
outsp.pred <- predict(outsp, X.0, xoord.0)

# Producing an SDM for HBEF alone (posterior mean)
mean.psisp = apply(outsp.pred$psi.0.samples, 2, mean)
sd.psisp = apply(outsp.pred$psi.0.samples, 2, sd)

psisp <- ggplot() + 
  geom_sf(data = grid, aes(fill = mean.psisp), lwd = 0.1) +
  scale_fill_viridis_c() + 
  labs(title = 'Spatial Occupancy') +
  theme_bw()

sdpsisp <- ggplot() + 
  geom_sf(data = grid, aes(fill = sd.psisp), lwd = 0.1) +
  scale_fill_viridis_c(option = "B") + 
  labs(title = 'Spatial Occupancy',
       subtitle= "standard error") +
  theme_bw()


library(patchwork)

(psi | psisp) / (sdpsi | sdpsisp)
psisp + sdpsisp
