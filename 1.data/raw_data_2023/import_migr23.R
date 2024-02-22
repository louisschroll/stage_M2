# # install.packages("archive")
setwd("~/stage_M2/1.data/raw_data")
# library(archive)
# archive_extract(
#   "Donnees_MIGRALION_2023.7z",
#   dir = "data2023",
#   files = NULL
# )

# Process Lot 4 sea surveys data
library(here)
library(tidyverse)
library(sf)
library(readxl)
load(here("~/stage_M2/1.data/gdlmap.rdata"))
#load(here("medmap.rdata"))

# ---- OBSERVATIONS ----
postnup_file = "Export_DATA_bateau_suivi_visuel_MigraLion_2023_postnuptiale.xlsx"
# prenup_file = "Export_DATA_bateau_suivi_visuel_MigraLion_2023_prenuptiale.xlsx"
# postnup_file = prenup_file
imp <- read_xlsx(postnup_file,
                 sheet= "DATA distance sampling") %>% 
  filter(!is.na(X_project)) %>% 
  mutate(x = X_project, y = Y_project) %>% 
  st_as_sf(coords = c("x", "y"), crs= 2154) %>% 
  st_transform(crs= st_crs(gdlmap))


# ---- EFFORT from data ---- 

effort.xls <- read_xlsx(postnup_file,
                        sheet= "DATA Effort")

# start effort
debut.eff <- effort.xls %>% 
  select(WKT_debut_effort, Time_UTC_debut_effort, Transect_effort) %>% 
  st_as_sf(wkt = "WKT_debut_effort", crs= 2154) %>% 
  st_transform(crs= st_crs(gdlmap)) %>% 
  mutate(position = "start") %>% 
  rename(geometry = "WKT_debut_effort",
         Time = Time_UTC_debut_effort)%>% 
  mutate(day = lubridate::date(Time))

# end effort
fin.eff <- effort.xls %>% 
  select(WKT_fin_effort, Time_UTC_fin_effort, Transect_effort) %>% 
  st_as_sf(wkt = "WKT_fin_effort", crs= 2154) %>% 
  st_transform(crs= st_crs(gdlmap)) %>% 
  mutate(position = "end") %>% 
  rename(geometry = "WKT_fin_effort",
         Time = Time_UTC_fin_effort) %>% 
  mutate(day = lubridate::date(Time))

# create linestring
eff <- bind_rows(debut.eff,fin.eff) %>% 
  group_by(Transect_effort, day, .drop = F) %>%
  summarise() %>%
  st_cast("LINESTRING") 

  ggplot()+
    #geom_sf(data = eff, aes(color = as_factor(Transect_effort)))+
    geom_sf(data = eff[which(eff$Transect_effort =="T3"),], aes(color = as_factor(as.character(day))), lwd = 2) +
    geom_sf(data = imp %>% filter(Date_UTC == as.Date("2023/09/25")))

  ggplot()+
   # geom_sf(data = eff[which(eff$Transect_effort =="T3"),], aes(color = as_factor(as.character(day))), lwd = 2) +
    geom_sf(data = eff %>% filter(day == as.Date("2023/09/25")),col='darkcyan', lwd = 2) +
    geom_sf(data = imp %>% filter(Date_UTC == as.Date("2023/09/25")))
    
# ---- WEATHER conditions ----
  
meteo.xls <- read_xlsx(postnup_file,
                          sheet= "DATA Meteo")

meteo <- meteo.xls %>% 
  st_as_sf(coords = c("X_tablette","Y_tablette"), crs= 2154) %>% 
  st_transform(crs= st_crs(gdlmap))

names(meteo)
names(imp)
unique(imp$Espece)
# 

if(FALSE){
  ggplot()+
    geom_sf(data =gdlmap) +
    geom_sf(data = eff, color ="red") +
    geom_sf(data = meteo, color = "darkgreen", alpha = 0.2, lwd = 1 ) +
    theme_minimal() +
    labs(title = "Observations visuelles",
         subtitle = "Lot 4 - postnuptiale 2023",
         caption = "Data : Biotope, Migralion")+
    theme(plot.title = element_text(face= "bold"),
          plot.subtitle = element_text(face= "italic"))
  
  ggsave("lot4obs.png", path = getwd(), dpi =300)
}

# save obs + effort + weather from datatable (No GPX trace)
if(FALSE){
  postnup2023_obs <- imp
  postnup2023_eff <- eff
  postnup2023_weather <- meteo
  
  save(postnup2023_weather,
       postnup2023_obs,
       postnup2023_eff,
       file = here("postnup2023.rdata"))
}

# if(FALSE){
  postnup2023_obs <- imp
  postnup2023_eff <- eff
  postnup2023_weather <- meteo

  save(postnup2023_weather,
       postnup2023_obs,
       postnup2023_eff,
       file = here("postnup2023.rdata"))
# }


#load(here("Lot3_Telemetrie/Data/caugek/imp_caugek.rdata"))

head(imp)
names(imp)

imp %>% 
  as.data.frame() %>% 
  group_by(Espece) %>% 
  count(sort =T) %>% 
  print(n=51)


# --- EFFORT FROM GPX ----
common_path = "data2023/Campagnes_postnuptiale_2023/Traces GPS/"
campagne1 = "2023_postnuptiale_Campagne1/"
campagne2 = "2023_postnuptiale_Campagne2/"

# campagne 1
p1 <- gpx::read_gpx(paste0(common_path, campagne1, "2023-09-06 08.20.33 Jour.gpx"))
p2 <- gpx::read_gpx(paste0(common_path, campagne1, "2023-09-07 00.00.36 Jour.gpx"))
p3 <- gpx::read_gpx(paste0(common_path, campagne1, "2023-09-08 00.00.01 Jour.gpx"))
p4 <- gpx::read_gpx(paste0(common_path, campagne1, "2023-09-09 00.00.16 Jour.gpx"))

# campagne 2 
d1 <- gpx::read_gpx(paste0(common_path, campagne2, "2023-09-24 17.01.33 Jour.gpx"))
d2 <- gpx::read_gpx(paste0(common_path, campagne2, "2023-09-25 00.00.20 Jour.gpx"))
d3 <- gpx::read_gpx(paste0(common_path, campagne2, "2023-09-26 00.00.28 Jour.gpx"))
d4 <- gpx::read_gpx(paste0(common_path, campagne2, "2023-09-26 14.31.49 Day.gpx"))



c1 <- bind_rows(p1$tracks$`2023-09-06 08:20:33 Jour`, 
                p2$tracks$`2023-09-07 00:00:36 Jour`,
                p3$tracks$`2023-09-08 00:00:01 Jour`,
                p4$tracks$`2023-09-09 00:00:16 Jour`) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = st_crs(gdlmap)) %>% 
  mutate(day = lubridate::date(Time)) 

c2 <- bind_rows(d1$tracks$`2023-09-24 17:01:33 Jour`,
                d2$tracks$`2023-09-25 00:00:20 Jour`,
                d3$tracks$`2023-09-26 00:00:28 Jour`,
                d4$tracks$`2023-09-26 14:31:49 Jour`) %>%
  st_as_sf(coords = c("Longitude","Latitude"), crs = st_crs(gdlmap)) %>% 
  mutate(day = lubridate::date(Time))


if(FALSE){
  ggplot() + 
    geom_sf(data = gdlmap) +
    geom_sf(data = ttp, aes(color = date))
}

# 8 days 
unique(c(c1$day, c2$day))

# create linestring
ttp <- bind_rows(c1, c2) %>% 
  group_by(day) %>%
  summarise(do_union = FALSE) %>%
  st_cast("LINESTRING")

# ---- EXAMPLE EXTRACT CAUGEK TERNS ----
caugek_sp2022 <- imp %>% 
  filter(Espece =="Sterne caugek" )

if(FALSE){ # plot terns obs
  ggplot()+
    geom_sf(data =gdlmap) +
    geom_sf(data = ttp, aes(color = as_factor(as.numeric(day)))) +
    geom_sf(data = caugek_sp2022 , color = "darkgreen", alpha = 0.9, lwd = 1 ) +
    theme_minimal() +
    scale_color_discrete() + 
    labs(title = "Observations visuelles - Sternes caugek",
         subtitle = "Lot 4 - prénuptiale 2022",
         caption = "Data : Biotope, Migralion")+
    theme(plot.title = element_text(face= "bold"),
          plot.subtitle = element_text(face= "italic"))
}


# ---- EFFORT from c1/c2 ----

head(c1)

# linestring with all points : nights and days
prenup_eff <- c1 %>% 
  bind_rows(c2) %>% 
  group_by(day) %>%
  summarise() %>%
  st_cast("LINESTRING") 

# keep all track locations
pt_eff <- c1 %>% bind_rows(c2, c3)

# ---- FILTER EFFORT to visual obs ----

# for 2023 they provided start or end of effort. 
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


metadata <- "postnup2023: every observation detected during campagins 2023.
                ."
prenup22_obs <- imp 
prenup22_eff <- tr_eff
save(prenup22_eff, prenup22_obs, metadata, file= "postnup2023.rdata")


