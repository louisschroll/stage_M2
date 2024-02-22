rm(list = ls())

# load necessary package
library(INLA)
library(RISDM)
library(terra)
library(sf)
library(tidyverse)
load("~/stage_M2/1.data/all_seabirds_counts.rdara")

# Get covars and simulated data
covars_raster <- rast("~/stage_M2/1.data/static_covariates_raster.tif")
#crs(covars_raster) <- crs(pelmed_obs)
plot(covars_raster)
names(covars_raster)
# simulated_data <- simulateData.isdm(rasterCovars = covars_raster[[c("bathymetry", "dist_to_shore")]],
#                          rasterBiasCovar = covars_raster[["slope"]],
#                          control = list(doPlot = F, set.random.seed=T, random.seed = 123456))

# Create a mesh
meshy <- makeMesh(covars_raster$slope, max.n = c(1000, 500),
                  dep.range = 2,
                  offset = 10,
                  expans.mult = 7.5,
                  doPlot = F)

checkMesh(meshy)
pelmed_obs2 <- pnm_obs %>% 
  rename(X = lat, Y = long) %>% 
  as_tibble() %>% 
  filter(nom_fr == "goeland leucophee") %>% 
  mutate(PA = ifelse(effectif>1, 1, 0),
         transect = seq(1,nrow(pelmed_obs2)))

migralion_obs2 <- st_as_sf(migralion_obs, coords = c("X", "Y")) %>% 
  st_transform(crs = "EPSG:4326") %>% 
  st_coordinates() %>% 
  cbind(migralion_obs) %>% 
  as_tibble() %>% 
  filter(nom_fr == "goeland leucophee") %>% 
  mutate(transect = seq(1,nrow(migralion_obs2)))


migralion_eff2 <- st_transform(migralion_eff, crs = "EPSG:4326")

df <- covars_raster %>% as.data.frame() %>% bind_cols(crds(covars_raster))
ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = bathymetry)) +
  scale_fill_distiller(palette = "Spectral", name = "T°") +
  geom_point(data = migralion_obs2, aes(x = X, y = Y, size = Effectif)) +
  geom_sf(data = migralion_eff2)


# Fit model
fm <- isdm(observationList = list(AAdat = migralion_obs2,
                                  PAdat = pelmed_obs2),
           covars = covars_raster,
           mesh = meshy,
           responseNames = c(AA="Effectif", PA="PA"),
           sampleAreaNames = c(AA="transect", PA="transect"),
           distributionFormula = ~0+bathymetry+dist_to_shore,
           artefactFormulas = list(AA=~1, PA=~1),
           control = list(prior.range = c(0.5, 0.1),
                          prior.space.sigma = c(2, 0.1),
                          coord.names = c("X", "Y")))

fm <- isdm(observationList = list(AAdat = pelmed_obs2),
           covars = covars_raster,
           mesh = meshy,
           responseNames = c(AA="effectif"),
           sampleAreaNames = c(AA="transect"),
           distributionFormula = ~0+bathymetry+dist_to_shore,
           artefactFormulas = list(AA=~1),
           control = list(prior.range = c(0.5, 0.1),
                          prior.space.sigma = c(2, 0.1),
                          coord.names = c("X", "Y")))


summary(fm)

# check residuals
plot(fm, nFigRow = 2, ask = F)

# Prediction
fm$preds <- predict(object = fm, covars = covars_raster, S = 1000)

tmprast <- fm$preds$field[[c(2,1,3)]]

plot(tmprast, 
     range = range(values(tmprast), na.rm=TRUE),
     col = hcl.colors(25, "viridis", rev = TRUE),
     nc = 3)

