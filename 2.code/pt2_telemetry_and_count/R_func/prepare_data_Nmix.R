#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/prepare_data_nmix.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-04
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------


prepare_data_Nmix <- function(data_list, grid, species, selected_cov){
  n.occ.cov <- length(selected_cov) + 1
  n.det.cov <- 2
  
  nmix_tibble_list <- map(data_list, 
                     ~{filter_one_species(.x, species) %>% 
                         get_count_and_effort(grid)
                       }) %>% 
    Filter(f=function(x) nrow(x)!=0) # remove dataset with 0 obs
  
  result_list <- map(nmix_tibble_list, ~ {
    .x %>%
      select(starts_with("transectLength"), starts_with("effectif"), cells_id) %>%
      relocate(cells_id) %>%
      pivot_longer(-c(cells_id), names_to = c(".value", "session"), names_pattern = "(.+)_(.+)") %>%
      # Keep only cells that have been sampled
      filter(!is.na(transectLength)) %>% 
      # Scale transect length (mean = 0 and sd = 1)
      mutate(transectLength = (transectLength - mean(transectLength))/sd(transectLength))
  }) %>% 
    bind_rows(.id = "dataset_nb") %>% 
    mutate(dataset_nb = as.numeric(as.factor(dataset_nb)))
  
  # Get the observation data
  data <- list(nobs = result_list$effectif)
  
  # Set initial values
  N0 <- result_list %>% group_by(cells_id) %>% 
    summarise(N0 = sum(effectif)+1) %>% 
    arrange(by = "cells_id") %>% 
    pull(N0)
  
  ndatasets <- max(result_list$dataset_nb)
  
  initial.values <- list(beta = rnorm(n.occ.cov,0,1), 
                         alpha = matrix(rnorm(n.det.cov*ndatasets), nrow = n.det.cov),
                         N = N0)
  
  # Get the constants
  constants <- list(XN = get_occ_cov(nmix_tibble_list, selected_cov),
                    dataset_nb = result_list$dataset_nb,
                    transect_length = result_list$transectLength,
                    site_id = get_new_id(result_list, id_column="cells_id"),
                    n.occ.cov = n.occ.cov,
                    nsites_total = n_distinct(result_list$cells_id),
                    nsampled_points = nrow(result_list),
                    ndatasets = ndatasets)
  
  return(list(data = data, constants = constants, inits = initial.values))
}


filter_one_species <- function(data, species){
  obs <- data$obs %>% filter(species_name == species)
  session_w_obs <- obs %>% pull(session) %>% unique()
  effort_data = data$eff %>% filter(session %in% session_w_obs)
  return(list(obs = obs, eff = effort_data))
}

# filter_all_dataset <- function(dataset_list, species){
#   filtered_dataset_list <- map(dataset_list, ~ filter_one_dataset(.x, species)) 
#   
#   removed_dataset = names(dataset_list)[unlist(map(filtered_dataset_list, is.null))]
#   if (length(removed_dataset) > 0){
#     print(paste0(removed_dataset, 
#                  " contain no observation for this species (",
#                  species, ")"))
#   }
#   # remove NULL values with compact (when no observation are made in a dataset)
#   return(compact(filtered_dataset_list))
# }



get_count_and_effort <- function(data, grid){
  #' Get the count data and compute transect length for all sampled cells
  #' @data: a list containing two sf data frame, obs with the observed counts and
  #' eff which contains the transects. Both data frame should have a column
  #' "session" indicating the corresponding session for each transect or observation
  #' @grid: a sf data frame containing the grid, with a column "id_cells"
  #' @Output: same format as grid with additionnal columns (transect_length and 
  #' effectif for each session)
  
  intersect_grid_obs <- st_intersection(data$obs, grid) 
  intersect_grid_eff <- st_intersection(data$eff, grid) 
  
  id_sampled_cells <- intersect_grid_eff %>% pull(cells_id) %>% unique() %>% sort()
  
  sampled_gridcells <- grid %>% 
    filter(cells_id %in% id_sampled_cells) %>% 
    as_tibble() %>% 
    select(-grid_c)
  
  # Get session names (only when the species is observed)
  session_names <- data$obs %>% pull(session) %>% unique()
  
  for (k in seq_along(session_names)){
    sessionK <- session_names[k]
    # transect length
    intersect_eff_sessionK <- intersect_grid_eff %>% 
      filter(session == sessionK) 
    
    df_transect_length <- tibble(cells_id = intersect_eff_sessionK$cells_id, 
                                 len = st_length(intersect_eff_sessionK)) %>% 
      group_by(cells_id) %>% 
      summarise(transectLength = sum(len))
    
    # Retrieve effectif 
    effectif_df <- intersect_grid_obs %>%
      as_tibble() %>% 
      filter(session == sessionK) %>% 
      select(cells_id, effectif) %>% 
      group_by(cells_id) %>% 
      summarise(effectif = sum(effectif)) %>% 
      arrange(cells_id)
    
    # Add in the grid
    sampled_gridcells <- sampled_gridcells %>% 
      left_join(df_transect_length, by = join_by(cells_id)) %>% 
      left_join(effectif_df, by = join_by(cells_id)) %>% 
      mutate(effectif = ifelse(is.na(effectif) & !is.na(transectLength), 0, effectif)) %>%
      rename_with(~ paste0(.x, "_", as.character(sessionK), recycle0 = TRUE), 
                  .cols = c(effectif, transectLength))
  }
  return(sampled_gridcells)
}


get_occ_cov <- function(nmix_tibble, selected_cov){
  map(nmix_tibble, function(x) select(x, all_of(c(selected_cov, "cells_id")))) %>% 
    bind_rows() %>% 
    unique() %>% 
    arrange(by = cells_id) %>% 
    select(all_of(selected_cov)) %>% 
    mutate(intersect = 1) %>% 
    relocate(intersect) %>% 
    as.matrix()
}





