# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2 _Copie/2.code/pt2_telemetry_and_count/prepare_data_for_RSF.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-03-28
#
# Script description:
# From .csv file download on movebank, this script convert the data to a dataframe
# with presence and pseudo-absence points needed for performing RSF. 
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(sf)

# Load data
load("grid.rdata")

prepare_data_for_RSF <- function(input_file, output_file){
  #' @input_file: path to the csv file with movebank data
  #' @output_file: path to the file in which the data will be saved
  
  convert_sf_to_df <- function(sf_object){
    sf_object %>% st_coordinates() %>% 
      bind_cols(st_drop_geometry(sf_object)) 
  }
  
  ## create a data frame with presence points (case=1) from the input file 
  df_RSF_case1 <- read.csv(input_file) %>% 
    as_tibble() %>% 
    rename(individual_id = individual.local.identifier,
           X = location.long,
           Y = location.lat,
           time = timestamp) %>% 
    select(X, Y, individual_id, time) %>% 
    # convert to sf to use st_intersection
    filter(!is.na(X) & !is.na(Y)) %>% 
    st_as_sf(dim = "XY", coords = c("X", "Y"), remove = TRUE,
             na.fail = TRUE, crs = st_crs(grid$geometry)) %>% 
    st_transform(crs = st_crs(grid)) %>% 
    # keep only points in the study area
    st_intersection(grid %>% select(-geometry)) %>% 
    # reconvert to tibble for more efficiency
    convert_sf_to_df() %>% 
    # add columns link to time 
    mutate(month = month(time),
           year = year(time),
           time = as_datetime(time)) %>% 
    # resample data every hour
    mutate(hourly_time = floor_date(time, "hour")) %>% 
    group_by(individual_id, hourly_time) %>%                          
    slice(1) %>% 
    ungroup() %>% 
    select(-c(hourly_time, time)) %>% 
    # add case
    mutate(case = 1)  
  
  ## Create a data frame for pseudo-absence (case = 0)
  # generate available points for pseudo absence (evenly accross the study area)
  n_random_points <- nrow(df_RSF_case1) * 10
  
  sf_available_points <- grid %>% 
    sample_n(size = n_random_points, replace = TRUE) %>% 
    st_centroid() 
  
  df_available_points <- st_coordinates(sf_available_points) %>% 
    as_tibble() %>% 
    bind_cols(sf_available_points %>% st_drop_geometry() %>% select(-geometry)) %>% 
    mutate(case = 0) %>% 
    relocate(case)
  
  # Assign id to pseudo absence (useful for random effect)
  # number of step for each individual for each month
  step_per_ind <- df_RSF_case1 %>% 
    group_by(individual_id, month, year) %>% 
    summarise(nstep = sum(case)) %>% 
    ungroup() %>% 
    filter(nstep > 0) %>% 
    mutate(n_random_points = nstep * 10) %>% 
    uncount(n_random_points) %>% 
    select(-nstep)
  
  shuffled_step_per_ind <- step_per_ind[sample(nrow(step_per_ind)), ]
  
  df_RSF_case0 <- bind_cols(shuffled_step_per_ind, df_available_points) %>% 
    relocate(X, Y)
  
  ## Bind the two df
  df_RSF <- df_RSF_case1 %>% bind_rows(df_RSF_case0) %>% 
    mutate(individual_id = as.numeric(as_factor(individual_id)))
  
  # ggplot(df_RSF) + 
  #   geom_sf(data = grid %>% st_union()) +
  #   geom_point(aes(x=X, y=Y, color = as.factor(case)))
  
  # ggplot() + 
  #   geom_sf(data = grid %>% st_union()) +
  #   geom_sf(data = sf_available_points[1000:2000,], color = "blue")
  
  # ggplot(df_RSF_case0) +
  #   geom_point(aes(x = X, y = Y, col = as.factor(individual_id))) +
  #   theme(legend.position = "none")
  
  save(df_RSF, file = output_file)
}

# prepare_data_for_RSF(input_file = "0.raw_data/Yellow-Legged Gull - France - Espagne -ID_PROG 990.csv", 
#                      output_file = "1.data/RSF_data_yellow_legged_gull.rdata")
# 
# prepare_data_for_RSF(input_file = "0.raw_data/Thalasseus sandvicensis Med - CRBPO 1190.csv", 
#                      output_file = "1.data/RSF_data_sandwich_tern.rdata")
# 
# prepare_data_for_RSF(input_file = "0.raw_data/Puffinus yelkouan - Yelkouan shearwater - Port-Cros France - ID_PROG 1190.csv", 
#                      output_file = "1.data/RSF_puffin_yelkouan.rdata")
# 
# prepare_data_for_RSF(input_file = "0.raw_data/Navigation in Scopoli's shearwaters (data from Pollonara et al. 2015).csv", 
#                      output_file = "1.data/RSF_puffin_scopoli_pollonara.rdata")
# 
# prepare_data_for_RSF(input_file = "0.raw_data/Calonectris diomedea - Scopoli's shearwater - Riou Marseille France - ID_PROG 1190.csv", 
#                      output_file = "1.data/RSF_puffin_scopoli_marseille.rdata")


