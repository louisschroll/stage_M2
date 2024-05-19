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


add_session_column <- function(data_obs, data_eff) {
  asso_session_date <- unique(data_eff %>% as_tibble %>% select(session, date))
  data_obs$session <- NA
  for (sessionId in unique(asso_session_date$session)){
    date_session <- asso_session_date %>% filter(session == sessionId) %>% pull(date) %>% as.Date()
    data_obs <- data_obs %>% mutate(session = ifelse(date %in% date_session, sessionId, session))
  }
  data_obs <- data_obs %>% mutate(session = as.factor(session))
  return(data_obs)
}


split_periods <- function(data_obs, session_summer, session_winter){
  # Create a new column indicating the species and the season
  data_obs <- data_obs %>% 
    mutate(species_name = case_when(
      session %in% session_summer ~ paste0(str_replace_all(nom_fr, " ", "_"), "_summer"),
      session %in% session_winter ~ paste0(str_replace_all(nom_fr, " ", "_"), "_winter")
    ))
  # get species viewed less than 5 times or viewed at only one survey (no repicate)
  species_to_exclude <- data_obs %>% 
    as_tibble() %>%
    count(species_name) %>%
    filter(n<=5) %>% 
    pull(species_name)
  
  species_to_exclude2 <- data_obs %>% 
    as_tibble() %>%
    group_by(species_name) %>% 
    summarise(n_replicate = length(unique(session))) %>% 
    filter(n_replicate == 1) %>% 
    pull(species_name)
  
  species_to_exclude <- unique(species_to_exclude, species_to_exclude2)
  
  # Filter NA, indeterminate species and species without enough data
  data_obs <- data_obs %>% 
    filter(!(species_name %in% species_to_exclude)) %>% 
    filter(!is.na(species_name)) %>% 
    filter(str_detect(string = species_name, pattern = "ind", negate = T))
  
  return(data_obs)
}

# rename_species <- function(data_obs, tibble_name){
#   data_obs %>% 
#     mutate(nom_fr = case_when(
#       nom_fr %in% tibble_name$old_name ~ tibble_name$new_name[match(nomCite, tibble_name$old_name)],
#       TRUE ~ nom_fr
#     )) 
# }


# Pelmed 2017 -> 2020
load("1.data/pelmed2017_2021.rdata")
load("1.data/pelmed.rdata")

ref_coordinate_system <- st_crs(pelmed_obs) 

pelmed_obs <- pelobs %>%  
  select(lat, lon, podSize, transect, routeType, sighting, date, 
         hhmmss, famille_fr, nom_fr, geometry) %>% 
  format_name(ref_coordinate_system) %>%
  mutate(nom_fr = case_when(
    nom_fr == "puffin cendre / de scopoli" ~ "puffin de scopoli",
    nom_fr == "grand puffin ind" ~ "puffin de scopoli",
    nom_fr == "goeland gris ind" ~ "goeland leucophee",
    nom_fr == "grand goeland ind" ~ "goeland leucophee",
    nom_fr == "petit puffin ind" ~ "petit puffin",
    nom_fr == "puffin des baleares" ~ "petit puffin",
    nom_fr == "puffin yelkouan" ~ "petit puffin",
    nom_fr == "oceanite ind" ~ "oceanite tempete",
    TRUE ~ nom_fr  
  )) %>%
  rename(effectif = podSize) %>% 
  mutate(session = as.factor(year)) %>% 
  split_periods(session_summer = c("2017","2018", "2019", "2021"), 
                session_winter = c(""))


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
         session = as.factor(nom_suivi %>% str_remove("_")),
         transect_name = 1:nrow(.),
         Time = scale(as.Date(DATE_TIME1))) %>% 
  st_transform(crs = ref_coordinate_system) %>% 
  filter(session != "SAMM_1(2summer)")

samm_obs <- birdsamm %>% 
  select(transect, passage, flight, date, hhmmss, podSize, observer, 
         lat, lon, speed, altitude, famille_fr, nom_fr, geometry) %>% 
  rename(effectif = podSize) %>% 
  format_name(ref_coordinate_system) %>%
  mutate(nom_fr = case_when(
    nom_fr == "puffin cendre / de scopoli" ~ "puffin de scopoli",
    nom_fr == "grand puffin ind" ~ "puffin de scopoli",
    nom_fr == "grand goeland ind" ~ "goeland leucophee",
    nom_fr == "goeland gris ind" ~ "goeland leucophee",
    nom_fr == "petit puffin ind" ~ "petit puffin",
    nom_fr == "oceanite ind" ~ "oceanite tempete",
    nom_fr == "pingouin ou guillemot" ~ "alcide ind",
    TRUE ~ nom_fr  
  )) %>%
  add_session_column(samm_eff) %>% 
  filter(session != "SAMM_1(2summer)") %>% 
  split_periods(session_summer = c(""), 
                session_winter = c("SAMM1(1winter)", "SAMM2(1winter)"))

# PNM Golfe du Lion ----
load("~/stage_M2/1.data/megaobs.rdata")

pnm_eff <- transect %>% 
  select(date_tr, long_m, camp, geometry) %>% 
  rename(date = date_tr) %>% 
  mutate(year = year(date),
         session = as.factor(camp %>% str_remove("_")),
         transect_name = 1:nrow(.)) 

pnm_obs <- obs_oiseaux %>% 
  select(espece, famille, date, heure, lat, long, nb, geometry) %>% 
  mutate(date = dmy(date)) %>% 
  rename(nom_fr = espece,
         effectif = nb) %>% 
  format_name(ref_coordinate_system) %>% 
  mutate(nom_fr = case_when(
    nom_fr == "puffin yelkouan de mediterranee" ~ "petit puffin",
    nom_fr == "puffin des baleares" ~ "petit puffin",
    nom_fr == "goeland" ~ "goeland leucophee",
    nom_fr == "puffin" ~ "puffin ind",
    TRUE ~ nom_fr
  )) %>% 
  add_session_column(pnm_eff) %>% 
  split_periods(session_summer = c("2019printemps", "2020printemps", "2021printemps"), 
                session_winter = c("2019automne", "2020automne"))


# Migralion 2022, 2023 ----
load("~/stage_M2/1.data/prenup2022_v3.rdata")
load("~/stage_M2/1.data/postnup2022.rdata")
load("~/stage_M2/1.data/postnup2023.rdata")
load("~/stage_M2/1.data/prenup2023.rdata")

session_migralion <- c(rep("prenup2022", 7), 
                       rep("postnup2022", 11), 
                       rep("prenup2023", 13), 
                       rep("postnup2023", 12))

migralion_eff <- prenup22_eff %>% 
  bind_rows(postnup2022_eff, prenup2023_eff, postnup2023_eff) %>% 
  mutate(date = day) %>% 
  mutate(year = year(date),
         session = as.factor(session_migralion),
         transect_name = 1:nrow(.)) 

migralion_obs <- prenup22_obs %>%
  select(Espece, Effectif, geometry, Date_UTC, Time_UTC, Heure_UTC) %>% 
  bind_rows(postnup2022_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC, Heure_UTC)) %>% 
  bind_rows(postnup2023_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC, Heure_UTC)) %>% 
  bind_rows(prenup2023_obs %>% 
              select(Espece, Effectif, geometry, Date_UTC, Time_UTC, Heure_UTC)) %>% 
  rename(nom_fr = Espece,
         date = Date_UTC,
         effectif = Effectif) %>% 
  mutate(date = as.Date(date)) %>% 
  format_name(ref_coordinate_system) %>% 
  mutate(nom_fr = case_when(
    nom_fr == "puffin yelkouan" ~ "petit puffin",
    nom_fr == "puffin des baleares" ~ "petit puffin",
    nom_fr == "puffin yelkouan/baleares" ~ "petit puffin",
    nom_fr == "puffin yelkouan/baleare" ~ "petit puffin",
    nom_fr %in% c("labbe pomarin", "labbe parasite", "labbe parasite/pomarin") ~ "labbe",
    TRUE ~ nom_fr  
  )) %>%
  add_session_column(migralion_eff) %>% 
  drop_na(session) %>% 
  split_periods(session_summer = c("prenup2022", "prenup2023"), 
                session_winter = c("postnup2022", "postnup2023")) %>% 
  mutate(species_name = case_when(
    species_name == "mouette_pygmee_summer" ~ "mouette_pygmee_winter",
    species_name == "oceanite_tempete_winter" ~ "oceanite_tempete_summer",
    species_name == "puffin_de_scopoli_winter" ~ "puffin_de_scopoli_summer",
    species_name == "macareux_moine_winter" ~ "macareux_moine_summer",
    TRUE ~ species_name
  ))

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
    TRUE ~ nomCite
  )) %>%
  select("nom_fr", "denombrementMax", "denombrementMin", "jourDateDebut") %>% 
  rename(date = jourDateDebut) %>%
  format_name(ref_coordinate_system) %>% 
  #add_session_column() %>% 
  mutate(nom_fr = case_when(
    nom_fr == "puffin cendre" ~ "puffin de scopoli",
    nom_fr == "puffin yelkouan" ~ "petit puffin",
    nom_fr == "puffin des baleares" ~ "petit puffin",
    TRUE ~ nom_fr
  )) %>% 
  mutate(date = as.Date(date),
         month = month(date)) %>% 
  mutate(species_name = case_when(
    month %in% c(1:3, 9:12) ~ paste0(str_replace_all(nom_fr, " ", "_"), "_winter"),
    month %in% c(4:8) ~ paste0(str_replace_all(nom_fr, " ", "_"), "_summer"),
  )) 

ferme_eff <- NULL


# Filter the species with enough data across the data sets and save ----
species_to_keep <- as_tibble(pelmed_obs) %>% count(species_name) %>% 
  full_join(as_tibble(samm_obs) %>% count(species_name), 
            by = "species_name") %>% 
  full_join(as_tibble(pnm_obs) %>% count(species_name), 
            by = "species_name") %>% 
  full_join(as_tibble(migralion_obs) %>% count(species_name), 
            by = "species_name") %>% 
  # full_join(as_tibble(ferme_obs) %>% count(species_name), 
  #           by = "species_name") %>% 
  mutate(row_sum = rowSums(select(., -species_name), na.rm = TRUE)) %>% arrange(row_sum) %>% print(n = 50)
  filter(row_sum > 95) %>% 
  pull(species_name)

migralion_obs <- migralion_obs %>% 
  filter(species_name %in% species_to_keep)

pnm_obs <- pnm_obs %>% 
  filter(species_name %in% species_to_keep)

samm_obs <- samm_obs %>% 
  filter(species_name %in% species_to_keep)

pelmed_obs <- pelmed_obs %>% 
  filter(species_name %in% species_to_keep)

# Save the data in a .rdata
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