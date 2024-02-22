


cities = data.frame(nom=c("Montpellier", "Perpignan", "Narbonne", "Marseille"), pays=rep("FR",4))
print(cities)

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

library(tidyverse)
coord = t(apply(cities, 1, function(aRow) locateCountry(aRow[1], aRow[2])))
cities_sf <- as_tibble(coord) %>% 
  rename(lat = V1, lon = V2) %>% 
  bind_cols(cities)# %>% 
  #st_as_sf(coords = c("lat", "lon"))

library("sf")
library("rnaturalearth")
library(ggspatial)

# https://r-spatial.org/r/2018/10/25/ggplot2-sf.html
world <- ne_countries(scale = "medium", returnclass = "sf")
st_crs(cities_sf) <- st_crs(world)
class(world)

area <- st_crop(x=world, xmin=2, ymin=40, xmax=8, ymax=45)

grid2 <- st_transform(grid, st_crs(world))
ggplot() + 
  # geom_text(data=cities_sf, aes(x=lat, y=lon, label = nom),
  #           color = "darkblue", fontface = "bold", check_overlap = FALSE) +
  geom_sf(data = grid2, aes(fill = mean.psi), lwd = 0.1) +
  scale_fill_viridis_c() + 
  #coord_sf(xlim = c(2, 8), ylim = c(40, 45), expand = TRUE) +
  geom_sf(data = world, fill= "antiquewhite") +
  # select an area
  # add scale (ggspatial package)
  annotation_scale(location = "bl", width_hint = 0.5) +
  # add noth arrow (ggspatial package)
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) + 
  annotate(geom = "text", x = 4.5, y = 42.7, label = "Lion Gulf", 
           fontface = "italic", color = "grey22", size = 6) +
  theme(panel.background = element_rect(fill = "aliceblue")) + 
  labs(title = 'Occupancy') 




## ---- osmdata ----

#first option (not recommended in this case)
library(osmdata)
library(tidyverse)
library(sf)

town <- 'London' 
location <- town %>% opq()

#second option (recommended)
# xmin, ymin, xmax, ymax
coords <- matrix(c(3, 42, 6, 43.7), byrow = TRUE, nrow = 2, ncol = 2, dimnames = list(c('x','y'),c('min','max'))) 
location <- coords %>% opq()
location <- opq(c(3, 42, 6, 43.7))
location <- getbb("mediterranean sea") %>% opq()
available_features()
water <- location %>%
  add_osm_feature(key = "natural", 
                  value = c("water")) %>%
  osmdata_sf()

ggplot() + geom_sf(data = water$osm_multipolygons, fill = 'light blue') + theme_minimal()


#build different types of streets
main_st <- data.frame(type = c("motorway","trunk","primary","motorway_junction","trunk_link","primary_link","motorway_link"))
st <- data.frame(type = available_tags('highway'))
st <- subset(st, !type %in% main_st$type)



