# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/import_and_format_data/prepare_obs_data.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-02-06
#
# Script description:
# Prepare the data containing the observations and effort for the analyses.
# Filter the species of interest, select the relevant colomns, add year and 
# session columns and homogenize columns and birds name. 
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# load packages
library(tidyverse)
library(sf)
library(stringi)

format_name <- function(dataset, crs_ref){
  selected_name <- c("goeland leucophee", "mouette pygmee", "sterne caugek", 
                     "puffin yelkouan","puffin ind", "mouette melanocephale",
                     "fou de bassan","oceanite tempete", "mouette tridactyle", 
                     "pingouin torda", "puffin de scopoli", "sterne pierregarin",
                     "sterne ind", "puffin des baleares", "grand cormoran",
                     "mouette rieuse","puffin yelkouan/baleares","plongeon arctique", 
                     "fou de bassan","macareux moine", "labbe parasite",
                     "puffin cendre", "laride ind", "mouette/sterne ind", 
                     "cormoran ind", "labbe parasite/pomarin", "labbe pomarin",
                     "goeland brun", "sterne/guifette ind", "cormoran huppe", 
                     "mouette ind", "alcide", "alcide ind", "goeland", 
                     "goeland gris ind", "grand goeland ind", "grand grebe ind",
                     "grand puffin ind", "mouette", "oceanite ind", 
                     "petit puffin", "petit puffin ind", "pingouin ou guillemot", 
                     "plongeon ind", "puffin", "puffin cendre / de scopoli", 
                     "sterne moyenne ind", "grand labbe", "sterne naine",
                     "goeland d'audouin", "puffin fuligineux", "sterne caindenne")
  
  dataset2 <- dataset %>% 
    mutate(nom_fr = stri_trans_general(str = nom_fr, id = "Latin-ASCII")) %>%  # remove accents
    mutate(nom_fr = tolower(nom_fr)) %>% # in lower case
    mutate(nom_fr = str_replace_all(nom_fr, pattern = "sp", replacement = "ind")) %>% 
    mutate(nom_fr = str_replace_all(nom_fr, pattern = "sp.", replacement = "ind")) %>% 
    mutate(nom_fr = str_replace_all(nom_fr, pattern = "ind.", replacement = "ind")) %>% 
    filter(nom_fr %in% selected_name) %>% 
    mutate(year = lubridate::year(date),
           month = lubridate::month(date))
  
  if (!(st_crs(dataset2) == crs_ref)){
    dataset2 <- dataset2 %>% st_transform(crs = crs_ref)
  }
  
  return(dataset2)
}


add_session_column2 <- function(data_obs, data_eff) {
  asso_session_date <- unique(data_eff %>% as_tibble %>% select(session, date))
  data_obs$session <- NA
  for (sessionId in unique(asso_session_date$session)){
    date_session <- asso_session_date %>% filter(session == sessionId) %>% pull(date) %>% as.Date()
    data_obs <- data_obs %>% mutate(session = ifelse(date %in% date_session, sessionId, session))
  }
  data_obs <- data_obs %>% mutate(session = as.factor(session))
  return(data_obs)
}


# Pelmed 2017 -> 2020
load("1.data/pelmed2017_2021.rdata")
load("1.data/pelmed.rdata")

ref_coordinate_system <- st_crs(pelmed_obs) 

pelmed_obs <- pelobs %>%  
  select(lat, lon, podSize, transect, routeType, sighting, date, 
         hhmmss, famille_fr, nom_fr, geometry) %>% 
  format_name(ref_coordinate_system) %>%
  mutate(nom_fr = ifelse(nom_fr=="puffin cendre / de scopoli", "puffin de scopoli", nom_fr)) %>%
  mutate(nom_fr = ifelse(nom_fr=="grand puffin ind", "puffin de scopoli", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="goeland gris ind", "goeland leucophee", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="grand goeland ind", "goeland leucophee", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="petit puffin ind", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin des baleares", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin yelkouan", "petit puffin", nom_fr)) %>%
  mutate(nom_fr = ifelse(nom_fr=="oceanite ind", "oceanite tempete", nom_fr)) %>%
  rename(effectif = podSize) %>% 
  mutate(session = as.factor(year))

pelmed_eff <- peleff %>% 
  select(effort, date, hhmmss, seaState, lat, lon, legLengKm, geometry) %>% 
  mutate(year = year(date),
         session = as.factor(year),
         transect_name = 1:nrow(.)) %>% 
  st_transform(crs = ref_coordinate_system)


# SAMM 2011/2012 - 2018-2019
load("~/stage_M2/1.data/birdSAMM.rdata")


samm_eff <- effsamm %>% 
  select(nom_suivi, flight, date, hhmmss, DATE_TIME1, seaState, swell, turbidity, skyGlint, 
         glareFrom, glareTo, glareSever, cloudCover, lat, lon, speed, altitude,
         aircraft, seaStIndex, geometry) %>% 
  mutate(year = lubridate::year(date),
         session = as.factor(nom_suivi),
         transect_name = 1:nrow(.),
         Time = scale(as.Date(DATE_TIME1))) %>% 
  st_transform(crs = ref_coordinate_system)

samm_obs <- birdsamm %>% 
  select(transect, passage, flight, date, hhmmss, podSize, observer, 
         lat, lon, speed, altitude, famille_fr, nom_fr, geometry) %>% 
  rename(effectif = podSize) %>% 
  format_name(ref_coordinate_system) %>%
  mutate(nom_fr = ifelse(nom_fr=="puffin cendre / de scopoli", "puffin de scopoli", nom_fr)) %>%
  mutate(nom_fr = ifelse(nom_fr=="grand puffin ind", "puffin de scopoli", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="grand goeland ind", "goeland leucophee", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="goeland gris ind", "goeland leucophee", nom_fr)) %>%
  mutate(nom_fr = ifelse(nom_fr=="petit puffin ind", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="oceanite ind", "oceanite tempete", nom_fr)) %>%
  add_session_column2(samm_eff)

# PNM Golfe du Lion
load("~/stage_M2/1.data/megaobs.rdata")

pnm_eff <- transect %>% 
  select(date_tr, long_m, camp, geometry) %>% 
  rename(date = date_tr) %>% 
  mutate(year = year(date),
         session = as.factor(camp),
         transect_name = 1:nrow(.)) 


pnm_obs <- obs_oiseaux %>% 
  select(espece, famille, date, heure, lat, long, nb, geometry) %>% 
  mutate(date = dmy(date)) %>% 
  rename(nom_fr = espece,
         effectif = nb) %>% 
  mutate(nom_fr = ifelse(nom_fr=="Puffin Yelkouan de Mediterranee", "petit puffin", nom_fr)) %>%
  format_name(ref_coordinate_system) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin des baleares", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="goeland", "goeland leucophee", nom_fr)) %>% 
  add_session_column2(pnm_eff)


# Migralion 2022, 2023
load("~/stage_M2/1.data/prenup2022_v3.rdata")
load("~/stage_M2/1.data/postnup2022.rdata")
load("~/stage_M2/1.data/postnup2023.rdata")
load("~/stage_M2/1.data/prenup2023.rdata")

session_migralion <- c(rep("prenup_2022_A", 3), rep("prenup_2022_B", 4), 
                       rep("postnup_2022_A", 5), rep("postnup_2022_B", 6), 
                       rep("prenup_2023_A", 6), rep("prenup_2023_B", 7), 
                       rep("postnup_2023_A", 7), rep("prenup_2023_B", 5))

session_migralion <- c(rep("prenup_2022", 7), 
                       rep("postnup_2022", 11), 
                       rep("prenup_2023", 13), 
                       rep("postnup_2023", 12))

migralion_eff <- prenup22_eff %>% 
  bind_rows(postnup2022_eff, prenup2023_eff, postnup2023_eff) %>% 
  mutate(date = day) %>% 
  mutate(year = year(date),
         session = as.factor(session_migralion),
         transect_name = 1:nrow(.)) 

migralion_obs <- prenup22_obs %>%
  select(Espece, Effectif, geometry, Date_UTC, Time_UTC) %>% 
  bind_rows(postnup2022_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC)) %>% 
  bind_rows(postnup2023_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC)) %>% 
  bind_rows(prenup2023_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC)) %>% 
  rename(nom_fr = Espece,
         date = Date_UTC,
         effectif = Effectif) %>% 
  mutate(date = as.Date(date)) %>% 
  format_name(ref_coordinate_system) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin yelkouan", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin des baleares", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin yelkouan/baleares", "petit puffin", nom_fr)) %>%
  mutate(nom_fr = ifelse(nom_fr=="puffin yelkouan/baleare", "petit puffin", nom_fr)) %>% 
  add_session_column2(migralion_eff) %>% 
  drop_na(session)

# migralion_obs %>% ggplot() +
#   geom_sf(aes(color = as.factor(year))) +
#   geom_sf(data = migralion_eff) +
#   facet_wrap(~session)

# fermes pilotes
load("~/stage_M2/1.data/pilotebiotope.rdata")

scientific_name <- c("Larus michahellis Naumann, 1840", 
                     "Hydrocoloeus minutus (Pallas, 1776)",
                     "Thalasseus sandvicensis (Latham, 1787)",
                     "Puffinus yelkouan (Acerbi, 1827)",
                     "Puffinus Brisson, 1760 sp.",
                     "Ichthyaetus melanocephalus (Temminck, 1820)",
                     "Morus bassanus (Linnaeus, 1758)",
                     "Hydrobates pelagicus (Linnaeus, 1758)",
                     "Rissa tridactyla (Linnaeus, 1758)",
                     "Alca torda Linnaeus, 1758",
                     "Calonectris diomedea (Scopoli, 1769)",
                     "Sterna hirundo Linnaeus, 1758",
                     "Sterninae",
                     "Puffinus mauretanicus Lowe, 1921",
                     "Phalacrocorax carbo (Linnaeus, 1758)",
                     "Chroicocephalus ridibundus (Linnaeus, 1766)",
                     "Puffin yelkouan / Baléares",
                     "Gavia arctica (Linnaeus, 1758)",
                     "Sula bassana (Linnaeus, 1758)",
                     "Fratercula arctica (Linnaeus, 1758)",
                     "Stercorarius parasiticus (Linnaeus, 1758)")

french_name <- c("goeland leucophee", 
                 "mouette pygmee",
                 "sterne caugek",
                 "petit puffin",
                 "puffin ind",
                 "mouette melanocephale",
                 "fou de bassan",
                 "oceanite tempete",
                 "mouette tridactyle",
                 "pingouin torda",
                 "puffin de scopoli",
                 "sterne pierregarin",
                 "sterne ind",
                 "puffin des baleares",
                 "grand cormoran",
                 "mouette rieuse",
                 "petit puffin",
                 "plongeon arctique",
                 "fou de bassan",
                 "macareux moine",
                 "labbe parasite")


ferme_obs <- efglx %>% 
  mutate(nom_fr = case_when(
    nomCite %in% scientific_name ~ french_name[match(nomCite, scientific_name)],
    TRUE ~ nomCite)) %>%
  select("nom_fr", "denombrementMax", "denombrementMin", "jourDateDebut") %>% 
  rename(date = jourDateDebut) %>%
  format_name(ref_coordinate_system) %>% 
  #add_session_column() %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin cendre", "puffin de scopoli", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin yelkouan", "petit puffin", nom_fr)) %>% 
  mutate(nom_fr = ifelse(nom_fr=="puffin des baleares", "petit puffin", nom_fr)) 
  

ferme_eff <- NULL


# Save the data in a .radta
save(pelmed_obs, pelmed_eff, 
     samm_obs, samm_eff, 
     pnm_obs, pnm_eff,
     migralion_obs, migralion_eff,
     ferme_obs, ferme_eff,
     file = "~/stage_M2/1.data/all_seabirds_counts.rdata")

# comptage / espece
# count_bird_obs <- function(data_obs){
#   count_df <- data_obs %>% 
#     as.data.frame() %>% 
#     group_by(nom_fr) %>% 
#     count(sort = TRUE) %>% 
#     filter(n>1)
#   return(count_df)
# }
# 
# b <- count_bird_obs(pelmed_obs) %>% rename(pelmed_count = n) 
# a <- count_bird_obs(migralion_obs) %>% rename(migralion_count = n)
# c <- count_bird_obs(samm_obs) %>% rename(samm_count = n)
# e <- count_bird_obs(ferme_obs) %>% rename(ferme_count = n)
# d <- count_bird_obs(pnm_obs) %>% rename(pnm_count = n)
# 
# count_df <- merge(a, b, by = "nom_fr", all = TRUE) %>% 
#   merge(c, by = "nom_fr", all = TRUE) %>% 
#   merge(d, by = "nom_fr", all = TRUE) %>% 
#   merge(e, by = "nom_fr", all = TRUE) %>% 
#   arrange(desc(migralion_count), desc(ferme_count))
# 
# library(writexl)
# write_xlsx(count_df, path = "comptage_par_espece.xlsx")