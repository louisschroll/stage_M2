
# laod data
load("gdlmap.rdata")
load("raw_data/prenup22.rdata")

# load packages
library(lubridate)
library(tidyverse)
library(sf)
# ---- EFFORT from c1/c2/c3 ----

head(c1)

# linestring with all points : nights and days
prenup_eff <- c1 %>% 
  bind_rows(c2, c3) %>% 
  group_by(day) %>%
  summarise() %>%
  st_cast("LINESTRING") 

# keep all track locations
pt_eff <- c1 %>% bind_rows(c2, c3)

# ---- FILTER EFFORT to visual obs ----

# they did not provide start or end of effort. 
# Then we can decide to start/end the effort at the time of the first/last detection of each day
imp %>% names()

# prepare a mask with obs data, i.e. collected when they where "en effort"
# we round the obs time to minute to compare with track locations
mask <- imp %>% 
  mutate(day = lubridate::day(Date_UTC),
         month = lubridate::month(Date_UTC),
         time = round_date(ymd_hms(Time_UTC), "minute")) %>% 
  select(day, month, geometry, time) 


# filter all track locations when obs have been recorded
# i.e. remove night locations
tr_eff <- pt_eff %>% 
  mutate(rtime = round_date(Time, "minute")) %>% 
  filter(rtime %in% mask$time) %>% st_transform(crs = st_crs(gdlmap)) %>% 
  group_by(day, .drop = F) %>% 
  summarise(do_union = FALSE) %>% 
  st_cast(to = "LINESTRING")

# voilou
if(F){
  ggplot() + 
    geom_sf(data = gdlmap) +
    geom_sf(data = imp)+
    geom_sf(data = tr_eff, color ="red") 
}


prenup22_obs <- imp 
prenup22_eff <- tr_eff
save(prenup22_eff, prenup22_obs , file= "prenup22_v2.rdata")
  