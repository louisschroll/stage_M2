# load packages

library(tidyverse)
library(sf)
library(stars)
library(amt)
library(nimble)
library(here)
library(patchwork)

# defin coolors palette  https://coolors.co/5f0f40-9a031e-fb8b24-e36414-0f4c5c
blue <-  "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <-  "#81b29a"

setwd("~/stage_M2/2.code/code_valentin")
#### Nmix #####

# Data ----
# Pelmed 2017-2021
load("~/stage_M2/1.data/pelmed2017_2021.rdata")
# peleff <- st_read(here("Data/Historiques/PELMED/PELMED2021/DonneesPelagis_Effort_2023-06-22.shp"))
# pelobs <- st_read(here("Data/Historiques/PELMED/PELMED2021/DonneesPelagis_Observations_2023-06-22.shp"))

load("contour_golfe_du_lion.Rdata")
load("~/stage_M2/1.data/gdlhex.rdata")
load("grid.rdata")

peleff <- peleff %>%
  select(seaState, date, geometry) %>%
  mutate(
    day = lubridate::day(date),
    month = lubridate::month(date),
    year = lubridate::year(date)
  ) %>%
  st_transform(crs = st_crs(sea.gdl2))

pelobs <- pelobs %>%
  dplyr::select(date, species, nom_fr, podSize, hhmmss, geometry) %>%
  filter(nom_fr == "Sterne caugek") %>%
  mutate(
    day = lubridate::day(date),
    month = lubridate::month(date),
    year = lubridate::year(date)
  ) %>%
  st_transform(crs = st_crs(sea.gdl2))

# 281 obs en 5 ans
pelobs %>%
  group_by(year) %>%
  count()

pelobs %>%
  group_by(year) %>%
  summarise(sum(podSize))

pelobs %>%
  pull(month) %>%
  unique()

# check data
if (F) {
  ggplot() + geom_sf(data = grid, alpha = 0) +
    geom_sf(data = contour_golfe) +
    geom_sf(data = peleff, color = red) +
    geom_sf(data = pelobs, color = purple) +
    theme_minimal() +
    labs(title = "Sterne caugek - nidification",
         subtitle = "GPS 2022, PELMED 2017-2021") +
    theme(plot.title = element_text(family = "Helvetica", face = "bold"))
}

pelobs
peleff

# name grid cells
grid <- sea.gdl2 %>%
  mutate(id = 1:nrow(sea.gdl2))

intpel <- st_intersection(peleff, grid)
pelsites <- unique(intpel$id)

if (F) {
  ggplot() +
    geom_sf(data = grid[migsites, ]) +
    geom_sf(data = c) +
    geom_sf(data = gdlmap)
}

ypel <- grid[pelsites, ] %>%
  mutate(
    idpel = 1:length(pelsites),
    ytot = 0,
    y2017 = 0,
    y2018 = 0,
    y2019 = 0,
    y2020 = 0,
    y2021 = 0,
    efftot = NA,
    eff2017 = NA,
    eff2018 = NA,
    eff2019 = NA,
    eff2020 = NA,
    eff2021 = NA
  )
intypel <- st_intersection(pelobs, ypel)
inteffpel <- st_intersection(peleff %>%
                               mutate(year = lubridate::year(date),
                                      month = month(date)), ypel)
inteffpel$length <- st_length(inteffpel)

### fill seff
for (i in 1:nrow(ypel)) {
  id <- which(inteffpel$idpel == ypel$idpel[i])
  ypel$efftot[i] <- sum(inteffpel$length[id])
  ypel$eff2017[i] <-
    sum(inteffpel$length[id][year(inteffpel$date)[id] == "2017"])
  ypel$eff2018[i] <-
    sum(inteffpel$length[id][year(inteffpel$date)[id] == "2018"])
  ypel$eff2019[i] <-
    sum(inteffpel$length[id][year(inteffpel$date)[id] == "2019"])
  ypel$eff2020[i] <-
    sum(inteffpel$length[id][year(inteffpel$date)[id] == "2020"])
  ypel$eff2021[i] <-
    sum(inteffpel$length[id][year(inteffpel$date)[id] == "2021"])
  
}

### fill count

#### fill w/ 0
ypel$y2017[ypel$eff2017 > 0] <- 0
ypel$y2018[ypel$eff2018 > 0] <- 0
ypel$y2019[ypel$eff2019 > 0] <- 0
ypel$y2020[ypel$eff2020 > 0] <- 0
ypel$y2020[ypel$eff2020 > 0] <- 0
ypel$y2021[ypel$eff2021 > 0] <- 0

#### fill w/ count
for (i in 1:nrow(intypel)) {
  y <-  intypel$year[i]
  ypel[intypel$idpel[i], paste0("y", y)] <-
    ypel[intypel$idpel[i], paste0("y", y)] + intypel$podSize[i]
  ypel$ytot[intypel$idpel[i]] <-
    ypel$ytot[intypel$idpel[i]] + intypel$podSize[i]
}

# dist_to_coast
load("~/stage_M2/1.data/gdlmap.Rdata")
dist_to_coast <-
  unlist(map2(sea.gdl2$x, st_union(gdlmap), st_distance))

ypel <- ypel %>%
  left_join(
    y = grid %>%
      st_drop_geometry() %>%
      mutate(dist_coast = dist_to_coast),
    by = "id"
  )

##### N-mix model model ####

# fit in Bayesian to be consistent

#code
nmix.ni <- nimbleCode({
  # priors
  for (i in 1:4) {
    a[i] ~ dnorm(0, 1)
  }
  
  for (i in 1:2) {
    b[i] ~ dnorm(0, 1)
  }
  
  # latent occu
  log(lambda[1:nsites]) <- a[1] + a[2] * prof[1:nsites] +
    a[3] * prof[1:nsites] * prof[1:nsites] +
    a[4] * dcol[1:nsites]
  
  # observation process
  logit(p[1:nsites, 1:nocc]) <- b[1] + b[2] * seff[1:nsites, 1:nocc]
  
  for (i in 1:nsites) {
    # likelihood
    N[i] ~ dpois(lambda[i])
    
    for (j in 1:nocc) {
      nobs[i, j] ~ dbin(p[i, j], N[i])
    }
  }
  
})

tibble(
  a = ypel$eff2017,
  b = ypel$eff2018,
  c = ypel$eff2019,
  d = ypel$eff2020,
  e = ypel$eff2021,
  row = 1:length(ypel$eff2017)
) %>%
  pivot_longer(-row) %>%
  mutate(value = (value - mean(value)) / sd(value)) %>%
  # distinct() %>%
  pivot_wider(names_from = name, values_from = value)

# format data
constants.o <- list(
  prof = scale(ypel$depth.x)[, 1],
  dcol = scale(ypel$dist_coast)[, 1],
  seff = scale(
    cbind(
      ypel$eff2017,
      ypel$eff2018,
      ypel$eff2019,
      ypel$eff2020,
      ypel$eff2021
    )
  ),
  nsites = nrow(ypel),
  nocc = 5
)

yy <- ypel %>%
  dplyr::select(starts_with("y20")) %>%
  st_drop_geometry()

data.o <- list(nobs = yy)

Ninit <- apply(yy, 1, sum) + 1
inits.o <- list(a = rnorm(4, 0, 1),
                b = rnorm(2, 0, 1),
                N = Ninit)

Rmodelo <- nimbleModel(
  code = nmix.ni,
  constants = constants.o,
  data = data.o,
  inits = inits.o
)

# start the nimble process
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
resNmix <-
  runMCMC(
    Cmcmco,
    niter = 510000,
    nburnin = 50000,
    nchains = 2,
    samplesAsCodaMCMC = TRUE
  )

# check convergence
mcmcplots::denplot(resNmix)
coda::effectiveSize(resNmix)

MCMCvis::MCMCsummary(object = resNmix, round = 2)

MCMCvis::MCMCtrace(
  object = resNmix,
  pdf = FALSE,
  ind = TRUE,
  Rhat = TRUE,
  n.eff = TRUE,
  params = "b"
)

# store results
res_b0ipp <- rbind(resNmix$chain1, resNmix$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[1]"))

res_profipp <- rbind(resNmix$chain1, resNmix$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[2]"))

res_qprofipp <- rbind(resNmix$chain1, resNmix$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[3]"))

res_dcolipp <- rbind(resNmix$chain1, resNmix$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("a[4]"))

#### RSF #####
# note that you can go directly to paragraph 'format for RSF 2' to save time
# raw data
load(here("Work/SSFOccu/Caugek/data/caugek06062023.rdata"))

# data description
head(imp)
names(imp)

imp2 %>%
  group_by(`individual-local-identifier`) %>%
  count() %>%
  print(n = 22)
# 22 individuals with several 1000"s of points

# crop to Golfe du lion
imp3 <- imp2 %>%
  st_crop(st_bbox(gdlmap))


if (F) {
  #plot tracks
  tracks <- imp3 %>%
    mutate(id = `individual-local-identifier`) %>%
    group_by(id) %>%
    summarise(do_union = F) %>%
    st_cast(to = "LINESTRING")
  
  ggplot() + geom_sf(data = gdlmap) +
    geom_sf(data = imp3[1, ])
}

# crop to Golfe du lion
imp3 <- imp2 %>%
  st_crop(st_bbox(gdlmap))

# rename individual id and select the columns
dat <- imp3 %>%
  mutate(id = `individual-local-identifier`) %>%
  dplyr::select(timestamp, id, geometry)


# data description: dates of tracks
dat %>%
  pull(timestamp) %>%
  lubridate::year() %>%
  unique()

dat %>%
  filter(lubridate::year(timestamp) == 2022) %>%
  pull(timestamp) %>%
  lubridate::month() %>%
  unique()

dat %>%
  filter(lubridate::year(timestamp) == 2022 &
           lubridate::month(timestamp) %in% c(6, 7, 8)) %>%
  group_by(id) %>%
  count()

# filter for summer 2022
dat <- dat %>%
  filter(lubridate::year(timestamp) == 2022 &
           lubridate::month(timestamp) %in% 6:8)

# filter for nesting individuals in 2021
# dat <- dat  %>%
#   filter(id %in% c("2021_Thau_central_02",
#                    "2021_Thau_central_04",
#                    "2021_Thau_central_05",
#                    "2021_Thau_central_06"))


# remove points on land
datt <- st_intersects(dat, st_union(gdlmap))
t <- apply(datt, 1, any)

dat2 <- dat[which(t == F), ]
remove(datt, t)

xy <- st_coordinates(dat2)

dat2 <- cbind(dat2, xy)

# remove land locations and make_track()
trk <- dat2 %>%
  st_as_sf(coords = c("X", "Y"), crs = st_crs(gdlmap)) %>%
  st_difference(st_union(gdlmap)) %>%
  st_drop_geometry() %>%
  make_track(X, Y, timestamp, id, crs = st_crs(dat))

if (F) {
  trkplot <- trk %>%
    st_as_sf(coords = c("x_", "y_"), crs = st_crs(gdlmap))
  
  ggplot() +
    geom_sf(data = gdlmap) +
    geom_sf(data = trkplot, aes(color = id)) +
    theme(legend.position = "none")
  # nidif n3, 4, 5,6,7,8, 9 thau, 10,12, 13, 15
}

trk1 <- trk %>%
  nest(data = -"id")

# resample data every hour
trk2 <- trk1 %>%
  mutate(steps = map(data, function(x)
    x %>% track_resample(rate = minutes(60), tolerance = minutes(5)) %>%  steps_by_burst()))

trk3 <- trk2 %>%
  dplyr::select(id, steps) %>%
  unnest(cols = steps)

# data description
trk3 %>%
  group_by(id) %>%
  count()

trk3 %>% nrow()

trk3 %>%
  pull(dt_) %>%
  min()

# plot step length (sl_) distribution
trk3 %>%
  ggplot() + geom_density(aes(x = sl_))

if (F) {
  ggplot() +
    geom_sf(data = gdlmap) +
    geom_sf(data = trk3 %>%
              st_as_sf(coords = c("x1_", "y1_"), crs = st_crs(gdlmap)), aes(color = id)) +
    theme(legend.position = "none")
}

data_ni <- trk3 %>%
  nest(data = -"id")

# save(data_ni,file = here("Work/RSFIPP/Caugek/ete/dataRSF_v2.rdata"))

# ---- format for RSF 2 ----
# generate available data
load("dataRSF_v2.rdata")

# if not loaded
# load(here("Seabirds/gdlmap.rdata"))
# load(here("Seabirds/gdlhex.rdata"))

trk3 <- data_ni %>%
  unnest(data)

# create mask to generate available points.
mask <- sea.gdl2 %>%
  st_union()

# place colony
load(here("Work/birdcol.Rdata"))

col7 <- birdcol[which(birdcol$code_site == "THAU_7")[1], ]

# number of locations
nptTrue <- data_ni %>%
  unnest(cols = c(data)) %>%
  nrow()

# create raster
sigma = as.numeric(max(st_distance(
  trk3 %>%
    st_as_sf(coords = c("x2_", "y2_"), crs = st_crs(col7)), col7
)) / 2.74)
map = 1.5 / (pi * sigma ^ 2) * exp(-sqrt(3) * as.numeric(st_distance(sea.gdl2, col7)) /
                                     sigma) #exponential negative
map = map * (1 / 0.95)

# distance to coast
dist_to_coast <-
  unlist(map2(sea.gdl2$x, st_union(gdlmap), st_distance))

# simulate only when distance to coast  20 km
rpraster <- sea.gdl2 %>%
  mutate(value = map,
         dist = as.numeric(st_distance(sea.gdl2, col7)))

# generate random points all in one
nptRand <- nptTrue * 11

nullCoords2 <- rpraster %>%
  sample_n(size = nptRand,
           replace = T,
           weight = rpraster$value) %>%
  st_centroid()

# plot check raster and available points
if (F) {
  ggplot() + #geom_sf(data= rpraster, aes(fill = value), lwd = 0) +
    geom_sf(data = nullCoords2[1:1000, ]) +
    geom_sf(
      data = trk3 %>%
        mutate(prof = prof$depth) %>%
        filter(prof > -200) %>%
        st_as_sf(coords = c("x2_", "y2_"), crs = st_crs(col7)),
      color = "gold"
    ) +
    geom_sf(data = sea.gdl2 %>% filter(depth > -200), lwd = 0)
}

nullCoord3 <- nullCoords2 %>%
  as_tibble() %>%
  arrange(dist)

nullCoord3 <- nullCoord3[1:(0.95 * nrow(nullCoord3)), ]

rpts <- nullCoord3 %>%
  mutate(
    case = 0,
    id = "0",
    x = st_coordinates(nullCoord3 %>% st_as_sf())[, "X"],
    y = st_coordinates(nullCoord3 %>% st_as_sf())[, "Y"]
  )


dfRSF <- trk3 %>% mutate(case = 1,
                         x = x2_,
                         y = y2_) %>%
  dplyr::select(id, case, x, y) %>%
  st_drop_geometry()

# Create covariates

# distance to colony
distcol <- dfRSF %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(gdlmap)) %>%
  st_distance(col7)

# bathymetry
prof <- dfRSF %>%
  mutate(locid = 1:nrow(dfRSF)) %>%
  st_as_sf(coords = c("x", "y"), crs = st_crs(gdlmap)) %>%
  st_intersection(sea.gdl2 %>%
                    dplyr::select(depth)) %>%
  dplyr::select(depth, locid)

data3 <- dfRSF %>%
  mutate(locid = 1:nrow(dfRSF))  %>%
  left_join(prof, by = "locid")

# trick to limite that max depth = 200m
data3 <- data3 %>%
  bind_rows(rpts) %>%
  mutate(depth.trunc = case_when(depth < -200 ~ -200,
                                 depth >= -200 ~ depth)) %>%
  mutate(depth.sc = scale(depth.trunc)[, 1],
         qprof = depth.sc * depth.sc)

# dist to coastline
distcoast <-
  st_distance(data3  %>% st_as_sf(coords = c("x", "y"), crs = st_crs(gdlmap)) , st_union(gdlmap))
length(which(is.na(distcoast)))

datan <- data3 %>%
  mutate(distcoast = distcoast,
         dcoast.sc = scale(as.numeric(distcoast))[, 1]) %>%
  dplyr::select(id, case, x, y, locid, depth, depth.sc, distcoast, dcoast.sc) %>%
  as_tibble()


# save(datan, file = here("Work/RSFIPP/Caugek/ete/dataRSF_v2.rdata"))
# load(here("Work/RSFIPP/Caugek/ete/nidif22_datrsf60min_v2.rdata"))
# plot data RSF + occupancy

# number of step for each individual
steperind <- datan %>%
  group_by(id) %>%
  summarise(nstep = sum(case)) %>%
  filter(nstep > 0)

steperind <- steperind %>%
  mutate(nrpt = nstep * 10,
         idd = 1:nrow(steperind))

# assign id to case == 0 , useful for random effect
steperind %>% pull(nrpt) %>% sum()
datan %>% filter(case == 0) %>% nrow()

dsamp <-  datan %>%
  mutate(rofull = 1:nrow(datan))  %>%
  filter(id == "0")


for (i in 1:nrow(steperind)) {
  dsamp <-  dsamp %>%
    filter(id == "0") %>%
    mutate(rosamp = 1:nrow(dsamp))
  
  rsamp <- dsamp %>%
    slice_sample(n = steperind$nrpt[i], replace = F) %>%
    mutate(id = steperind$id[i])
  
  datan$id[rsamp$rofull] <- rsamp$id[1]
  dsamp <- dsamp[-rsamp$rosamp, ]
}

# data description
datan %>%
  group_by(id, case) %>%
  count()

datan %>%
  pull(id) %>%
  unique()
datan %>%
  filter(id == "0") %>%
  nrow()

# run RSF nimble -----

dataNimb <- datan %>%
  filter(id != "0") %>%
  mutate(idind = as.numeric(as_factor(id)))

nindividual = nrow(steperind)

npts = dataNimb %>% nrow()

##### code ####
rsf.caugek <- nimbleCode({
  betaprof_pop ~ dnorm(0, 1)
  tauprof_pop ~ dunif(0, 1e2)
  betaqprof_pop ~ dnorm(0, 1)
  tauqprof_pop ~ dunif(0, 1e2)
  betadcol_pop ~ dnorm(0, 1)
  taudcol_pop ~ dunif(0, 1e2)
  
  beta0_pop ~ dnorm(0, 1)
  
  for (i in 1:nindividual) {
    ### PRIORS ###
    ## habitat cov
    beta_prof[i] ~ dnorm(betaprof_pop, sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop, sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop, sd = taudcol_pop)
    # movement cov
    beta_0[i] ~ dnorm(beta0_pop, 1e2)
  }
  
  # likelihood
  for (t in 1:npts) {
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * prof[t] +
      beta_qprof[idind[t]] * prof [t] * prof[t] +
      beta_dcol[idind[t]] * dcol[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
  
})
#
#
# --- bundle and run ----
# w <- dataNimb$case
# w[w==0] <- 1000
# # constants
# constants.ni <-  list(npts = npts,
#                       idind = dataNimb$idind,
#                       prof = dataNimb$depth.sc,
#                       dcol = dataNimb$dcoast.sc,
#                       nindividual = nindividual,
#                       w = w)
#
# # data
# data.ni <- list(kase = dataNimb$case)
#
# # Inits
# inits.ni <-  list(beta_0 = rep(0, nindividual),
#                   beta_prof =  rep(0, nindividual),
#                   beta_qprof=  rep(0, nindividual),
#                   beta_dcol =  rep(0, nindividual),
#                   betaprof_pop = 0 ,
#                   beta0_pop = 0 ,
#                   tauprof_pop = runif(1,0,5),
#                   betaqprof_pop = 0 ,
#                   tauqprof_pop = runif(1,0,5),
#                   betadcol_pop = 0 ,
#                   taudcol_pop = runif(1,0,5))
#
#
# # Nimble pre run
# Rmodel2 <- nimbleModel(code= rsf.caugek, constants = constants.ni, data = data.ni, inits = inits.ni)
# Rmodel2$initializeInfo()
# Rmodel2$calculate() # - 58854147
#
# # configure model
# conf2 <- configureMCMC(Rmodel2)
# conf2$printMonitors()
# #
# ## Build and compile MCMC
# Rmcmc2 <- buildMCMC(conf2)
# Cmodel2 <- compileNimble(Rmodel2)
# Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)
#
#
# # Run
# samplesRSF <- runMCMC(Cmcmc2, niter = 110000, nburnin = 10000, thin = 5, nchains = 1, samplesAsCodaMCMC = TRUE)  ## DT: use runMCMC
#
# # check convergence
# mcmcplots::denplot(samplesRSF)
# coda::effectiveSize(samplesRSF)
#
# # store results
# # when only 1 chain
# res_prof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
#   as_tibble() %>%
#   dplyr::select(starts_with("betaprof_pop"))
#
# res_tauprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
#   as_tibble() %>%
#   dplyr::select(starts_with("tauprof_pop"))
#
# res_qprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
#   as_tibble() %>%
#   dplyr::select(starts_with("betaqprof_pop"))
#
# res_dcol <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
#   as_tibble() %>%
#   dplyr::select(starts_with("betadcol_pop"))
#
# # make prediction maps
# if(F){
# load(here("Seabirds/gdlmap.rdata"))
# load(here("Seabirds/gdlhex.rdata"))
# gridpred <- sea.gdl2 %>% filter(depth > -250)
# dist_to_coast <- unlist(map2(gridpred$x, st_union(gdlmap), st_distance))
#
#
# pred_bathy <- tibble(prof.sc =  gridpred$depth.sc , prof = gridpred$depth, dcoast = scale(dist_to_coast)[,1]) %>%
#   mutate(pred= prof.sc* mean(res_prof$betaprof_pop) +
#            prof.sc*prof.sc* mean(res_qprof$betaqprof_pop) +
#             0*dcoast* mean(res_dcol$betadcol_pop))
#
# ggplot(data = pred_bathy) +
#   geom_point(aes( x = prof, y = pred))
# }


#### INTEGRATED MODEL #####

#####  code ----
# with NIMBLE package
caugekRSFIPP <- nimbleCode({
  # IPP
  # priors
  intIPP ~ dnorm(0, 1)
  
  for (i in 1:2) {
    b[i] ~ dnorm(0, 1)
  }
  
  # latent
  log(lambda[1:nsites]) <- intIPP + betaprof_pop * prof[1:nsites] +
    betaqprof_pop * prof[1:nsites] * prof[1:nsites] +
    betadcol_pop * dcol[1:nsites]
  
  # observation process
  logit(p[1:nsites, 1:nocc]) <- b[1] + b[2] * seff[1:nsites, 1:nocc]
  for (i in 1:nsites) {
    # likelihood
    N[i] ~ dpois(lambda[i])
    
    for (j in 1:nocc) {
      nobs[i, j] ~ dbin(p[i, j], N[i])
    }
  }
  
  # integrated coefficients
  betaprof_pop ~ dnorm(0, 1)
  tauprof_pop ~ dunif(0, 1e2)
  betaqprof_pop ~ dnorm(0, 1)
  tauqprof_pop ~ dunif(0, 1e2)
  betadcol_pop ~ dnorm(0, 1)
  taudcol_pop ~ dunif(0, 1e2)
  
  beta0_pop ~ dnorm(0, 1)
  # intercept fixe
  #intRSF ~ dnorm(0,1)
  
  for (i in 1:nindividual) {
    ### PRIORS ###
    ## habitat cov
    beta_prof[i] ~ dnorm(betaprof_pop, sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop, sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop, sd = taudcol_pop)
    # intercept cov
    beta_0[i] ~ dnorm(beta0_pop, 1e6)
  }
  
  # ll
  for (t in 1:npts) {
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * profRSF[t] +
      beta_qprof[idind[t]] * profRSF [t] * profRSF[t] +
      beta_dcol[idind[t]] * dcolRSF[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
  
})
#

# data for integrated model
w <- dataNimb$case
w[w == 0] <- 1000


constants.int <-  list(
  #RSF
  npts = npts,
  idind = dataNimb$idind,
  profRSF = dataNimb$depth.sc,
  dcolRSF = dataNimb$dcoast.sc,
  nindividual = nindividual,
  w = w,
  # IPP
  prof = scale(ypel$depth.x)[, 1],
  dcol = scale(ypel$dist_coast)[, 1],
  # RV
  seff = scale(
    cbind(
      ypel$eff2017,
      ypel$eff2018,
      ypel$eff2019,
      ypel$eff2020,
      ypel$eff2021
    )
  ),
  nsites = nrow(ypel),
  nocc = 5
)

# data
yy <- ypel %>%
  dplyr::select(starts_with("y20")) %>%
  st_drop_geometry()

data.int <- list(kase = dataNimb$case,
                 nobs = yy)

# Inits
Ninit <- apply(yy, 1, sum) + 1
inits.int <-  list(
  beta_0 =     rep(0, nindividual),
  beta_prof =  rep(0, nindividual),
  beta_qprof =  rep(0, nindividual),
  beta_dcol =  rep(0, nindividual),
  betaprof_pop = 0 ,
  tauprof_pop = runif(1, 0, 5),
  betaqprof_pop = 0 ,
  tauqprof_pop = runif(1, 0, 5),
  betadcol_pop = 0 ,
  taudcol_pop = runif(1, 0, 5),
  beta0_pop = 0,
  intIPP = rnorm(1, 0, 1),
  b = rnorm(2, 0, 1),
  N = Ninit
)

# Nimble pre run
Rmodel3 <-
  nimbleModel(
    code = caugekRSFIPP,
    constants = constants.int,
    data = data.int,
    inits = inits.int
  )
Rmodel3$initializeInfo()
Rmodel3$calculate() # - 29993
# configure model
conf3 <- configureMCMC(Rmodel3)
## Build and compile MCMC
Rmcmc3 <- buildMCMC(conf3)

Cmodel3 <- compileNimble(Rmodel3)
Cmcmc3 <- compileNimble(Rmcmc3, project = Cmodel3)

# Run
samplesint <-
  runMCMC(
    Cmcmc3,
    niter = 11000,
    nburnin = 1000,
    nchains = 1,
    samplesAsCodaMCMC = TRUE
  ) ## DT: use runMCMC

# check convergence
mcmcplots::denplot(samplesint)
coda::effectiveSize(samplesint)

# add iterations in Nimble
# niter_ad <- 100000
# Cmcmc3$run(niter_ad, reset = FALSE)
# more_samples <- as.matrix(Cmcmc3$mvSamples)
# mcmcplots::denplot(more_samples)

# store results with 1 chain
res_b0int <-
  samplesint %>%  #rbind(samplesintrand$chain1, samplesintrand$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("intIPP"))

res_profint <-
  samplesint %>%  #rbind(samplesintrand$chain1, samplesintrand$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaprof_pop"))

res_tauprof21 <-
  samplesint %>%  #rbind(samplesintrand$chain1, samplesintrand$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("tauprof_pop"))


res_qprofint <-
  samplesint %>%  #rbind(samplesintrand$chain1, samplesintrand$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaqprof_pop"))


res_dcolint <-
  samplesint %>%  #rbind(samplesintrand$chain1, samplesintrand$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betadcol_pop"))


#THE PLOT #####
# plot intercept ----
#
plot_int <- res_b0ipp %>%
  mutate(b0 = `a[1]`,
         model = "counts")  %>%
  bind_rows(res_b0int %>% mutate(b0 = intIPP,
                             model = "integrated mod")) %>%
  ggplot(aes(model, b0, color = model)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)) +
  geom_boxplot(
    width = .25,
    outlier.shape = NA
  ) +
  labs(title = "Intercept",
       x= "",
       y = "") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(colour = 'black', size = 12),
        plot.title = element_text(size = 16, face = "bold", family  ="Helvetica", 
                                  hjust = 0.5),
        plot.subtitle = element_text(family  ="Helvetica", hjust = 0.5),
        legend.position = "right",
        legend.title = element_blank()) +
  guides(color = guide_legend(label.position = "top", keywidth = unit(1, "pt")),
         fill = "none")


 plot_int

# ggsave("intercept_caugek.png", plot = plot_int, path = here("Work/RSFIPP/Caugek/ete"), dpi = 200)


# ---- plot effect size ----
# load(here("Work/RSFIPP/Caugek/ete/nidif_resipp.rdata"))
# load(here("Work/RSFIPP/Caugek/ete/res_RSFIPPv30min.rdata"))
plot_prof <- tibble(model = "RSF",
                    mu_prof = res_prof$betaprof_pop) %>%
  bind_rows(
    tibble(
      mu_prof = res_profipp %>%
        slice_sample(n = 10000) %>%
        pull(`a[2]`),
      model = "Poisson GLM"
    ),
    #  tibble(mu_prof = res_profint$betaprof_pop, model = "Integrated model"),
    # tibble(mu_prof = res_profRSF21$betaprof_pop, model = "RSF21"),
    tibble(
      mu_prof = res_profint %>%
        slice_sample(n = 10000) %>%
        pull(betaprof_pop),
      model = "Integrated model"
    )
  ) %>%
  ggplot(aes(mu_prof, model, color = model)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  # geom_point(
  #   size = 1.3,
  #   alpha = .01,
  #   aes(color = model),
  #   position = position_jitter(
  #     seed = 1, width = .1
  #   )) +
  labs(
    title = "Linear bathymetry effect" ,
    subtitle = "on Sandwich tern space-use",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  scale_fill_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.text.y = element_blank(),
    plot.title = element_text(
      size = 16,
      face = "bold",
      family  = "Helvetica",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family  = "Helvetica", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank(),
    plot.title.position = "plot"
  ) +
  guides(color = guide_legend(
    label.position = "top",
    keywidth = unit(1, "pt"),
    nrow = 1
  ),
  fill = "none")


plot_prof

# ggsave(here("Work/SSFOccu/Caugek/ete/2indVS15ind.png"), plot = last_plot(), dpi  =300, width  = 30, height = 20, unit = "cm", bg = "white")

#

plot_qprof <- tibble(model = "RSF",
                     mu_qprof = res_qprof$betaqprof_pop) %>%
  bind_rows(
    tibble(mu_qprof = res_qprofipp$`a[3]`, model = "Poisson GLM"),
    tibble(mu_qprof = res_qprofint$betaqprof_pop, model = "Integrated model")
  ) %>%
  ggplot(aes(mu_qprof, model, color = model)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  # geom_point(
  #   size = 1.3,
  #   alpha = .01,
  #   aes(color = model),
  #   position = position_jitter(
  #     seed = 1, width = .1
  #   )) +
  labs(
    title = "Quadratic bathymetry effect",
    subtitle = "on Sandwich tern space-use",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  scale_fill_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.text.y = element_blank(),
    plot.title = element_text(
      size = 16,
      face = "bold",
      family  = "Helvetica",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family  = "Helvetica", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(
    label.position = "top",
    keywidth = unit(1, "pt"),
    nrow = 1
  ),
  fill = "none")

plot_qprof

plot_dcol <- tibble(model = "RSF",
                    mu_dcol = res_dcol$betadcol_pop) %>%
  bind_rows(
    tibble(mu_dcol = res_dcolipp$`a[4]`, model = "Poisson GLM"),
    tibble(mu_dcol = res_dcolint$betadcol_pop, model = "Integrated model")
  ) %>%
  ggplot(aes(mu_dcol, model, color = model)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  # geom_point(
  #   size = 1.3,
  #   alpha = .01,
  #   aes(color = model),
  #   position = position_jitter(
  #     seed = 1, width = .1
  #   )) +
  labs(
    title = "Distance to coast effect",
    subtitle = "on Sandwich tern space-use",
    x = "",
    y = ""
  ) +
  theme_minimal() +
  scale_color_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  scale_fill_manual(values = c(
    "RSF" = red,
    "Integrated model" = purple,
    "Poisson GLM" = blue
  )) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.text.y = element_blank(),
    plot.title = element_text(
      size = 16,
      face = "bold",
      family  = "Helvetica",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family  = "Helvetica", hjust = 0.5),
    legend.position = "bottom",
    legend.title = element_blank()
  ) +
  guides(color = guide_legend(
    label.position = "top",
    keywidth = unit(1, "pt"),
    nrow = 1
  ),
  fill = "none")
plot_dcol
#library(patchwork)
pp <-  (plot_prof / plot_qprof / plot_dcol) +
  plot_layout(guides = "collect", tag_level = "new") &
  theme(
    legend.position = "bottom",
    legend.text = element_text(
      family  = "Helvetica",
      face = "bold",
      size = 12
    )
  )
pp

## the maps ----
# plot
# grid predictions
load(here("Seabirds/gdlmap.rdata"))
load(here("Seabirds/gdlhex.rdata"))

gridpred <- sea.gdl2 %>%
  mutate(dcoast = scale(as.numeric(st_distance(
    sea.gdl2, st_union(gdlmap)
  )))[, 1],
  dist = dcoast) %>%
  st_crop(peleff %>% st_transform(crs = st_crs(sea.gdl2)) %>% st_buffer(dist = 10000)) %>%
  filter(depth > -200)

gdlpred <- gdlmap %>% st_crop(gridpred)

ggplot() + geom_sf(data = gdlpred) +
  geom_sf(data = gridpred, aes(fill = depth))

cor(gridpred$depth, gridpred$dcoast)

# IPP/NMix model
# omipp <- apply(cbind(res_profipp,res_qprofipp, res_dcolipp) %>%
#                  sample_n(size = 50000, replace = F),1, function(x = double(1)){
#   omega <-  x[1] * gridpred$depth.sc + x[2] * gridpred$depth.sc *gridpred$depth.sc + x[3] * gridpred$dcoast
#   return(omega)
# } )
#
# ripp_sd <- apply(omipp, 1, function(x) { r <- sd(exp(x))
#  return(r)})
# ripp <- apply(omipp, 1, function(x) { r <- mean(exp(x))
# return(r)})
# ripp_cv <- ripp_sd /ripp

pipp_m <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = ripp), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + 
  theme_minimal() +
  scale_fill_continuous(labels = c("Low", "High"),
                        breaks = c(54.83, 140.14)) +
  #scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", 30, direction = -1,
  #                       labels = c("Low", "High"), breaks = c(54.83,138.14))+
  labs(title = "Poisson GLM") +
  theme(
    legend.direction = "horizontal",
    plot.title.position = 'plot',
    legend.position = c(0.5, 1.5),
    plot.margin = margin(t = 4, 1, 1, 1, "lines"),
    plot.title = element_text(
      family = "Helvetica",
      face = "bold",
      hjust = 0.5
    )
  ) +
  guides(
    fill = guide_colourbar(
      title = "Relative space-use",
      title.theme = element_text(
        family = "Helvetica",
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      title.position = 'top',
      title.hjust = .5,
      barwidth = unit(10, 'lines'),
      barheight = unit(.5, 'lines')
    ),
    colour = "none"
  )


pipp_sd <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = ripp_sd), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
  theme(
    legend.direction = "horizontal",
    plot.title.position = 'plot',
    legend.position = c(0.5, 1.5),
    plot.margin = margin(t = 4, 1, 1, 1, "lines"),
    plot.title = element_text(
      family = "Helvetica",
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family = "Helvetica", hjust = 0.5)
  ) +
  labs(title = "Poisson GLM") +
  scale_fill_viridis_c(labels = c("Low", "High"),
                       breaks = c(54.5, 163)) +
  guides(
    fill = guide_colourbar(
      title = "Standard deviation",
      title.position = 'top',
      title.hjust = .5,
      title.theme = element_text(
        family = "Helvetica",
        face = "bold",
        size = 16,
        hjust = 0.5
      ),
      barwidth = unit(10, 'lines'),
      barheight = unit(.5, 'lines')
    ),
    colour = "none"
  )

# pipp_cv <-  ggplot() +
#   geom_sf(data = gridpred , aes(fill = ripp_cv), lwd = 0) +
#   geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
#   theme(legend.direction = "horizontal",
#         plot.title.position = 'plot',
#         legend.position =c(0.5,1.5),
#         plot.margin = margin(t=4,1,1,1, "lines"),
#         plot.title = element_text(family = "Helvetica", face = "bold", hjust = 0.5)) +
#   labs(title = "Poisson GLM" )+
#   scale_fill_viridis_c(labels = c("Low", "High"), breaks = c(0.842,1.19)) +
#   guides(fill = guide_colourbar(title = "Coefficient of variation",
#                                 title.position = 'top', title.hjust = .5,
#                                 title.theme = element_text(family = "Helvetica", face = "bold", size = 16, hjust = 0.5),
#                                 barwidth = unit(10, 'lines'), barheight = unit(.5, 'lines')),
#          colour = "none")
#
# pipp <- pipp_m | pipp_sd # | pipp_cv
#
# # RSF
# omrsf <- apply(cbind(res_prof %>%
#                        pull(betaprof_pop),
#                      res_qprof %>%
#                        pull(betaqprof_pop),
#                      res_dcol%>%
#                        pull(betadcol_pop)),1, function(x = double(1)){
#   omega <-  x[1] * gridpred$depth.sc + x[2] * gridpred$depth.sc *gridpred$depth.sc + x[3] * gridpred$dcoast
#   return(omega)
# } )
#
# rrsf_sd <- apply(omrsf, 1, function(x) { r <- sd(exp(x))
# return(r)})
# rrsf <- apply(omrsf, 1 ,function(x) { r <- mean(exp(x))
# return(r)})
# rrsf_cv <- rrsf_sd / rrsf
#
prsf_m <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = rrsf), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + 
  theme_minimal() +
  labs(title = "RSF") +
  # scale_fill_continuous(labels = c("Low", "High"), breaks = c(54.83,140.14))+
  #scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", 30, direction = -1)+
  theme(
    legend.position = "top",
    plot.title = element_text(
      family = "Helvetica",
      hjust = 0.5,
      face = "bold"
    )
  ) +
  guides(fill = "none",
         colour = "none")

prsf_sd <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = rrsf_sd), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
  labs(title = "RSF") +
  theme(
    legend.position = "top",
    plot.title = element_text(
      family = "Helvetica",
      hjust = 0.5,
      face = "bold"
    )
  ) +
  scale_fill_viridis_c() +
  guides(fill = "none",
         colour = "none")
#
# prsf_cv <-  ggplot() +
#   geom_sf(data = gridpred , aes(fill = rrsf_cv), lwd = 0) +
#   geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
#   labs(title = "RSF")+
#   theme(legend.position = "top",
#         plot.title = element_text(family = "Helvetica",hjust = 0.5, face = "bold")) +
#   scale_fill_viridis_c() +
#   guides(fill = "none",
#          colour = "none")
#
# prsf <- prsf_m  | prsf_sd # | prsf_cv #+ plot_annotation(title = "SSF model",
# # theme = theme(plot.title = element_text(family = "Helvetica", hjust = 0.5, face = "bold"))))


# integrated model
omint <- apply(cbind(res_profint, res_qprofint, res_dcolint) %>%
                 sample_n(size = 50000, replace = F), 1, function(x = double(1)) {
                   omega <-
                     x[1] * gridpred$depth.sc + x[2] * gridpred$depth.sc * gridpred$depth.sc + x[3] * gridpred$dcoast
                   return(omega)
                 })

rint_sd <- apply(omint, 1, function(x) {
  r <- sd(exp(x))
  return(r)
})
rint <- apply(omint, 1 , function(x) {
  r <- mean(exp(x))
  return(r)
})
rint_cv <- rint_sd / rint

pint_m <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = rint), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal() +
  labs(title = "Integrated model") +
  #scale_fill_paletteer_c("ggthemes::Orange-Blue-White Diverging", 30, direction = -1)+
  theme(
    legend.position = "top",
    plot.title = element_text(
      family = "Helvetica",
      hjust = 0.5,
      face = "bold"
    )
  ) +
  guides(fill = "none",
         colour = "none")

pint_sd <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = rint_sd), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
  labs(title = "Integrated model") +
  theme(
    legend.position = "top",
    plot.title = element_text(
      family = "Helvetica",
      hjust = 0.5,
      face = "bold"
    )
  ) +
  scale_fill_viridis_c() +
  guides(fill = "none",
         colour = "none")

pint_cv <-  ggplot() +
  geom_sf(data = gridpred , aes(fill = rint_cv), lwd = 0) +
  geom_sf(data = gdlpred, lwd = 0.1) + theme_minimal()  +
  labs(title = "Integrated model") +
  theme(
    legend.position = "top",
    plot.title = element_text(
      family = "Helvetica",
      hjust = 0.5,
      face = "bold"
    )
  ) +
  scale_fill_viridis_c() +
  guides(fill = "none",
         colour = "none")

pint <-
  pint_m |
  pint_sd #  | pint_cv + plot_annotation(title = "SSF - Occupancy model")
#theme = theme(plot.title = element_text(family = "Helvetica", hjust = 0.5, face = "bold"))))



# global plot maps
# pgmaps <- (pipp / prsf / pint)
#
# ggsave(here("Work/RSFIPP/Caugek/ete/resmaps_v5.png"), plot = pgmaps, dpi  =200, width  = 20, height = 20, unit = "cm")

# -----SD/CV PRECISION  distribution of space-use standard deviation ----
pred_precision <-
  tibble(value = ripp_cv, model = "Poisson GLM") %>%
  bind_rows(tibble(value = rrsf_cv, model = "RSF")) %>%
  bind_rows(tibble(value = rint_cv, model = "Integrated model")) %>%
  ggplot(aes(value, model, color = model)) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  # geom_point(
  #   size = 1.3,
  #   alpha = .01,
  #   aes(color = model),
  #   position = position_jitter(
  #     seed = 1, width = .1
  #   )) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Space-use precision",
    subtitle = "Predicted CV among grid-cells",
    x = "Coefficient of variation (CV)",
    y = ""
  ) +
  scale_color_manual(
    values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "Poisson GLM" = blue
    )
  ) +
  scale_fill_manual(
    values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "Poisson GLM" = blue
    )
  ) +
  theme(
    axis.text.x = element_text(colour = 'black', size = 12),
    axis.text.y = element_text(colour = 'black', size = 8),
    plot.title = element_text(
      size = 16,
      face = "bold",
      family  = "Helvetica",
      hjust = 0.5
    ),
    plot.subtitle = element_text(family  = "Helvetica", hjust = 0.5),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  guides(color = "none", #guide_legend(label.position = "top", keywidth = unit(1, "pt"), nrow = 1),
         fill = "none")

pred_precision
ggsave(
  here("Work/RSFIPP/Caugek/ete/predprecision_v1.png"),
  plot = pred_precision,
  dpi  = 200,
  width  = 20,
  height = 15,
  unit = "cm",
  bg = "white"
)


# mega plot final ----
library(patchwork)

pleft <- pgmaps

pright <- (pp / pred_precision)

pult <- (pleft | pright) +
  plot_layout(
    widths = c(2, 1),
    heights = c(1),
    guide = "keep"
  ) +
  plot_annotation(
    title = '',
    caption = 'Source: Migralion project & PELMED 2017-2021',
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      legend.position = 'top'
    )
  )


pult
ggsave(
  here("Work/RSFIPP/Caugek/ete/pglobal_v4.png"),
  plot = pult,
  dpi  = 300,
  width  = 50,
  height = 30,
  unit = "cm",
  bg = "white"
)
