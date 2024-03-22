#### RSF #####

# note that you can go directly to paragraph 'format for RSF 2' to save time
# raw data
load("1.data/caugek06062023.rdata")

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


if(F){
  #plot tracks
  tracks <- imp3 %>% 
    mutate(id = `individual-local-identifier`) %>% 
    group_by(id) %>% 
    summarise(do_union = F) %>% 
    st_cast(to = "LINESTRING")
  
  ggplot() + geom_sf(data = gdlmap) + 
    geom_sf(data = imp3[1,])
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
  filter(lubridate::year(timestamp)== 2022) %>% 
  pull(timestamp) %>% 
  lubridate::month() %>% 
  unique()

dat %>% 
  filter(lubridate::year(timestamp) ==2022 & 
           lubridate::month(timestamp) %in% c(6,7,8)) %>% 
  group_by(id) %>% 
  count()

# filter for summer 2022
dat <- dat %>% 
  filter(lubridate::year(timestamp) ==2022 & 
           lubridate::month(timestamp) %in% 6:8) 

# filter for nesting individuals in 2021
# dat <- dat  %>%
#   filter(id %in% c("2021_Thau_central_02",
#                    "2021_Thau_central_04",
#                    "2021_Thau_central_05",
#                    "2021_Thau_central_06"))


# remove points on land
datt <- st_intersects(dat,st_union(gdlmap))
t <- apply(datt, 1, any)

dat2 <- dat[which(t ==F),]
remove(datt, t)

xy <- st_coordinates(dat2)

dat2 <- cbind(dat2, xy)

# remove land locations and make_track()
trk <- dat2 %>% 
  st_as_sf(coords = c("X", "Y"), crs = st_crs(gdlmap)) %>% 
  st_difference(st_union(gdlmap)) %>% 
  st_drop_geometry() %>% 
  make_track(X, Y, timestamp, id, crs = st_crs(dat))

if(F){
  trkplot <- trk %>%  
    st_as_sf(coords = c("x_", "y_"), crs = st_crs(gdlmap))
  
  ggplot() + 
    geom_sf(data = gdlmap) +
    geom_sf(data = trkplot, aes(color = id)) + 
    theme(legend.position = "none")
  # nidif n3, 4, 5,6,7,8, 9 thau, 10,12, 13, 15
}

trk1 <- trk %>% 
  nest(data = - "id")

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
  ggplot() + geom_density(aes(x= sl_))

if(F){
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

col7 <- birdcol[which(birdcol$code_site == "THAU_7")[1],]

# number of locations
nptTrue <- data_ni %>% 
  unnest(cols = c(data)) %>% 
  nrow()

# create raster
sigma = as.numeric(max(st_distance(trk3 %>% 
                                   st_as_sf(coords = c("x2_","y2_"),crs= st_crs(col7)), col7))/2.74)
map = 1.5/(pi*sigma^2)*exp(-sqrt(3)*as.numeric(st_distance(sea.gdl2, col7))/sigma) #exponential negative
map = map*(1/0.95)

# distance to coast
dist_to_coast <- unlist(map2(sea.gdl2$x, st_union(gdlmap), st_distance))

# simulate only when distance to coast  20 km
rpraster <- sea.gdl2 %>% 
  mutate(value = map,
         dist = as.numeric(st_distance(sea.gdl2, col7)))

# generate random points all in one 
nptRand <- nptTrue * 11

nullCoords2 <- rpraster %>% 
  sample_n(size = nptRand, replace = T, weight = rpraster$value) %>% 
  st_centroid()

# plot check raster and available points
if(F){
  ggplot() + #geom_sf(data= rpraster, aes(fill = value), lwd = 0) + 
    geom_sf(data = nullCoords2[1:1000,])+ 
    geom_sf(data = trk3 %>%
              mutate(prof = prof$depth) %>% 
              filter(prof > -200) %>% 
              st_as_sf(coords = c("x2_","y2_"),crs= st_crs(col7)), color ="gold") + 
    geom_sf(data = sea.gdl2 %>% filter(depth > -200), lwd = 0)
}

nullCoord3 <- nullCoords2 %>% 
  as_tibble() %>% 
  arrange(dist) 

nullCoord3 <- nullCoord3[1:(0.95*nrow(nullCoord3)),]

rpts <- nullCoord3 %>%
  mutate(case = 0,
         id = "0",
         x = st_coordinates(nullCoord3 %>% st_as_sf())[,"X"],
         y = st_coordinates(nullCoord3 %>% st_as_sf())[,"Y"])


dfRSF <- trk3 %>% mutate(case = 1,
                         x=x2_,
                         y= y2_) %>% 
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
  mutate(depth.trunc = case_when(depth < - 200 ~ -200,
                                 depth >= -200 ~ depth)) %>% 
  mutate(depth.sc = scale(depth.trunc)[,1],
         qprof = depth.sc * depth.sc)

# dist to coastline
distcoast <-  st_distance(data3  %>% st_as_sf(coords = c("x","y"), crs = st_crs(gdlmap)) , st_union(gdlmap))
length(which(is.na(distcoast)))

datan <- data3 %>%
  mutate(distcoast = distcoast,
         dcoast.sc = scale(as.numeric(distcoast))[,1]) %>% 
  dplyr::select(id, case, x,y,locid, depth, depth.sc, distcoast, dcoast.sc) %>% 
  as_tibble()


# save(datan, file = here("Work/RSFIPP/Caugek/ete/dataRSF_v2.rdata"))
# load(here("Work/RSFIPP/Caugek/ete/nidif22_datrsf60min_v2.rdata"))
# plot data RSF + occupancy

# number of step for each individual
steperind <- datan %>%
  group_by(id) %>% 
  summarise(nstep = sum(case)) %>% 
  filter(nstep >0)

steperind <- steperind %>% 
  mutate(nrpt = nstep * 10,
         idd = 1:nrow(steperind))

# assign id to case == 0 , useful for random effect
steperind %>% pull(nrpt) %>% sum()
datan %>% filter(case ==0) %>% nrow()

dsamp <-  datan %>%
  mutate(rofull = 1:nrow(datan))  %>% 
  filter(id == "0") 


for(i in 1:nrow(steperind)){
  
  dsamp <-  dsamp %>%
    filter(id == "0") %>% 
    mutate(rosamp = 1:nrow(dsamp))
  
  rsamp <- dsamp %>% 
    slice_sample(n = steperind$nrpt[i], replace = F) %>% 
    mutate(id = steperind$id[i])
  
  datan$id[rsamp$rofull] <- rsamp$id[1]
  dsamp <- dsamp[- rsamp$rosamp,]
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
  
  betaprof_pop ~ dnorm(0,1)
  tauprof_pop ~ dunif(0,1e2)
  betaqprof_pop ~ dnorm(0,1)
  tauqprof_pop ~ dunif(0,1e2)
  betadcol_pop ~ dnorm(0,1)
  taudcol_pop ~ dunif(0,1e2)
  
  beta0_pop ~ dnorm(0,1)
  
  for( i in 1:nindividual){
    
    ### PRIORS ###
    ## habitat cov
    beta_prof[i] ~ dnorm(betaprof_pop, sd = tauprof_pop)
    beta_qprof[i] ~ dnorm(betaqprof_pop, sd = tauqprof_pop)
    beta_dcol[i] ~ dnorm(betadcol_pop, sd = taudcol_pop)
    # movement cov
    beta_0[i] ~ dnorm(beta0_pop,1e2)
  }
  
  # likelihood
  for(t in 1:npts){
    
    logit(omega[t]) <- beta_0[idind[t]] +
      beta_prof[idind[t]] * prof[t] + 
      beta_qprof[idind[t]] *prof [t] * prof[t] +
      beta_dcol[idind[t]] * dcol[t]
    
    kase[t] ~ dbinom(omega[t], w[t])
  }
})


# --- bundle and run ----
w <- dataNimb$case
w[w==0] <- 1000
# constants
constants.ni <-  list(npts = npts,
                      idind = dataNimb$idind,
                      prof = dataNimb$depth.sc,
                      dcol = dataNimb$dcoast.sc,
                      nindividual = nindividual,
                      w = w)

# data
data.ni <- list(kase = dataNimb$case)

# Inits
inits.ni <-  list(beta_0 = rep(0, nindividual),
                  beta_prof =  rep(0, nindividual),
                  beta_qprof=  rep(0, nindividual),
                  beta_dcol =  rep(0, nindividual),
                  betaprof_pop = 0 ,
                  beta0_pop = 0 ,
                  tauprof_pop = runif(1,0,5),
                  betaqprof_pop = 0 ,
                  tauqprof_pop = runif(1,0,5),
                  betadcol_pop = 0 ,
                  taudcol_pop = runif(1,0,5))


# Nimble pre run
Rmodel2 <- nimbleModel(code= rsf.caugek, constants = constants.ni, data = data.ni, inits = inits.ni)
Rmodel2$initializeInfo()
Rmodel2$calculate() # - 58854147

# configure model
conf2 <- configureMCMC(Rmodel2)
conf2$printMonitors()
#
## Build and compile MCMC
Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)


# Run
samplesRSF <- runMCMC(Cmcmc2, niter = 110000, nburnin = 10000, thin = 5, nchains = 1, samplesAsCodaMCMC = TRUE)  ## DT: use runMCMC

# check convergence
mcmcplots::denplot(samplesRSF)
coda::effectiveSize(samplesRSF)

# store results
# when only 1 chain
res_prof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaprof_pop"))

res_tauprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("tauprof_pop"))

res_qprof <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betaqprof_pop"))

res_dcol <- samplesRSF %>%  #rbind(samplesRSF$chain1, samplesRSF$chain2) %>%
  as_tibble() %>%
  dplyr::select(starts_with("betadcol_pop"))

# make prediction maps
if(F){
load(here("Seabirds/gdlmap.rdata"))
load(here("Seabirds/gdlhex.rdata"))
gridpred <- sea.gdl2 %>% filter(depth > -250)
dist_to_coast <- unlist(map2(gridpred$x, st_union(gdlmap), st_distance))


pred_bathy <- tibble(prof.sc =  gridpred$depth.sc , prof = gridpred$depth, dcoast = scale(dist_to_coast)[,1]) %>%
  mutate(pred= prof.sc* mean(res_prof$betaprof_pop) +
           prof.sc*prof.sc* mean(res_qprof$betaqprof_pop) +
            0*dcoast* mean(res_dcol$betadcol_pop))

ggplot(data = pred_bathy) +
  geom_point(aes( x = prof, y = pred))
}
