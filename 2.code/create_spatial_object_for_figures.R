
library(dplyr)
library(sf)
library(terra)

load("~/stage_M2/1.data/gdlmap.Rdata")
depth_raster <- terra::rast("1.data/static_covariates_raster.tif")

area_edge <- sf::st_crop(
  gdlmap %>% select(geometry, NAME),
  xmin = 640000,
  ymin = 6115601,
  xmax = 930000,
  ymax = 6350000
)

plot(area_edge)

# Get an elevation map around the Gulf of Lion (gdl=golfe du lion) ----
library(elevatr)

raster_elevation_gdl <- get_elev_raster(locations = area_edge,
                                   z = 12,
                                   clip = "locations",
                                   src="aws"
) %>%
  rast() %>%
  # ramener les valeurs négatives à zéro
  app(fun=function(x){ x[x < 0] <- 0; return(x)}) %>% 
  aggregate(fact=20)

# gdl_ombr_raster <- shade(terrain(raster_elevation_gdl, v='slope')/180*pi,
#                         terrain(raster_elevation_gdl, v='aspect')/180*pi,
#                         angle=40, direction=330
# )

# plot(gdl_ombr_raster, col=grey(0:100/100))

#gdl_ombr_raster_2 <- aggregate(gdl_ombr_raster, fact=2)

ggplot() +
  layer_spatial(raster_elevation_gdl) +
  scale_fill_gradientn(colors=terrain.colors(n=40), na.value = 0) +
  theme_bw()

writeRaster(raster_elevation_gdl, 
            file = "1.data/spatial_objects/raster_elevation_gdl.tif", 
            overwrite=TRUE)

crs(gdl_elev_raster) <- crs(gdlmap)

# Crop the map of bathymetry and save it ----
depth_raster <- project(depth_raster, y = crs(gdlmap)) %>% 
  crop(gdl_elev_raster)

writeRaster(depth_raster, 
            file = "1.data/spatial_objects/raster_depth.tif", 
            overwrite=TRUE)

# Save some cities in a sf object ----
library(RJSONIO)
locateCountry = function(nameCity, codeCountry) {
  cleanCityName = gsub(' ', '%20', nameCity)
  url = paste(
    "http://nominatim.openstreetmap.org/search?city="
    , cleanCityName
    , "&countrycodes="
    , codeCountry
    , "&limit=9&format=json"
    , sep="")
  resOSM = fromJSON(url)
  if(length(resOSM) > 0) {
    return(c(resOSM[[1]]$lon, resOSM[[1]]$lat))
  } else return(rep(NA,2)) 
}

cities = data.frame(nom=c("Montpellier", "Perpignan", "Narbonne", "Marseille"), pays=rep("FR",4))
print(cities)
coord = t(apply(cities, 1, function(aRow) locateCountry(aRow[1], aRow[2])))
cities_sf <- as_tibble(coord) %>% 
  rename(lat = V1, lon = V2) %>% 
  bind_cols(cities) %>% 
  st_as_sf(coords = c("lat", "lon")) %>% 
  st_set_crs(value=st_crs(depth_raster)) %>% 
  st_transform(crs = st_crs(area_edge))

save(cities_sf, file = "1.data/spatial_objects/cities_sf.RData")


# Save a sf of France and surrounding coutries
library(rnaturalearth)
sf_use_s2(FALSE)
sf_france <- st_crop(
  x = world,
  xmin = -4.5,
  ymin = 41,
  xmax = 8,
  ymax = 51
) %>% 
  filter(sovereignt != "United Kingdom") %>% 
  select(name_en)

plot(sf_france)
save(sf_france, file = "1.data/spatial_objects/sf_france.RData")
