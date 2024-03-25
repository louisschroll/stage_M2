


# ---- Prepare data for N-mixture ----
prepare_data_for_int_nmix <- function(data_list, grid, species, selected_cov){
  nmix_tibble_list <- map(data_list, ~{prepare_tibble_for_nmix(.x, grid, species)})
  
  common_sites <- map(nmix_tibble_list, function(x) select(x, all_of(c(selected_cov, "id")))) %>% 
    bind_rows() %>% 
    unique() %>% 
    arrange(by = id)
  
  occurence_cov_df <- common_sites %>% 
    select(-id) %>% 
    mutate(intersect = 1) %>% 
    relocate(intersect)
  
  id_correspondance <- tibble(old_id = common_sites$id, new_id = 1:nrow(common_sites))
  
  id_list <- map(nmix_tibble_list, 
                 ~{left_join(.x %>% select(id), id_correspondance, by = c("id" = "old_id")) %>% pull(new_id)})
  
  effectif_list <- map(nmix_tibble_list, ~{.x %>% select(starts_with("effectif"))})
  
  # effort_list <- map(nmix_tibble_list, ~{list(
  #   transect_length = .x %>% select(starts_with("transect_length")) %>% scale(),
  #   session = .x %>% select(starts_with("session")))})
  
  effort_list <- map(nmix_tibble_list, ~{.x %>% select(starts_with("transect_length"))})
  
   nmix_data_list <- list(
    effectif = Nmix_data_tibble %>% select(starts_with("effectif")),
    occurence.cov = occurence_cov_df,
    
    )
   
   nmix_constants_list <- list(
     sites_id = id_list,
     XN = occurence_cov_df,
     transect_length = effort_list,
     # detection.cov = list(
     #   transect_length = Nmix_data_tibble %>% select(starts_with("transect_length")),
     #   session = Nmix_data_tibble %>% select(starts_with("session")))
     nsites = sapply(effectif_list, nrow),
     nreplicates = sapply(effectif_list, ncol),
     n.occ.cov = ncol(occurence_cov_df)
   )
}
prepare_tibble_for_nmix <- function(data, grid, species){
  filtered_data <- filter_one_species(data, species)

  Nmix_data_tibble <- get_count_and_effort(filtered_data, grid) 

  # Nmix_data_list <- list(
  #   sites_id = id_correspondance$new_id,
  #   effectif = Nmix_data_tibble %>% select(starts_with("effectif")),
  #   occurence.cov = occurence_cov_df,
  #   detection.cov = list(
  #     transect_length = Nmix_data_tibble %>% select(starts_with("transect_length")),
  #     session = Nmix_data_tibble %>% select(starts_with("session")))
  # )
  
  return(Nmix_data_tibble)
}

filter_one_species <- function(data, species){
    obs = data$obs %>% filter(species_name == species)
    session_w_obs <- obs %>% pull(session) %>% unique()
    effort_data = data$eff %>% filter(session %in% session_w_obs)
    return(list(obs = obs, eff = effort_data))
}


get_count_and_effort <- function(data, grid){
  #' Get the count data and compute transect length for all sampled cells
  #' @data: a list cntaining two sf data frame, obs with the observed counts and
  #' eff which contains the transects. Both data frame should have a column
  #' "session" indicating the corresponding session for each transect or observation
  #' @grid: a sf data frame containing the grid, with a column "id_cells"
  #' @Output: same format as grid with additionnal columns (transect_length and 
  #' effectif for each session)
  intersect_grid_obs <- st_intersection(data$obs, grid) 
  intersect_grid_eff <- st_intersection(data$eff, grid) 
  
  id_sampled_cells <- intersect_grid_eff %>% pull(id) %>% unique() %>% sort()
  
  sampled_gridcells <- grid %>% 
    filter(id %in% id_sampled_cells) %>% 
    as_tibble() %>% 
    select(-grid_c)
  
  # Get session names (only when the species is observed)
  session_names <- data$obs %>% pull(session) %>% unique()
  
  for (k in seq_along(session_names)){
    sessionK <- session_names[k]
    # transect length
    intersect_eff_sessionK <- intersect_grid_eff %>% 
      filter(session == sessionK) 
    
    df_transect_length <- tibble(id = intersect_eff_sessionK$id, 
                                 len = st_length(intersect_eff_sessionK)) %>% 
      group_by(id) %>% 
      summarise(transect_length = sum(len))
    
    # effectif df
    effectif_df <- intersect_grid_obs %>%
      as_tibble() %>% 
      filter(session == sessionK) %>% 
      select(id, effectif) %>% 
      group_by(id) %>% 
      summarise(effectif = sum(effectif)) %>% 
      arrange(id)
    
    # Add in the grid
    sampled_gridcells <- sampled_gridcells %>% 
      left_join(df_transect_length, by = join_by(id)) %>% 
      left_join(effectif_df, by = join_by(id)) %>% 
      mutate(effectif = ifelse(is.na(effectif), 0, effectif),
             transect_length = ifelse(is.na(transect_length), 0, transect_length),
             session = sessionK) %>% 
      rename_with(~ paste0(.x, as.character(sessionK), recycle0 = TRUE), 
                  .cols = c(effectif, transect_length, session))
  }
  return(sampled_gridcells)
}



prepare_data_for_Nmixture <- function(data, grid, species, selected_cov){
  filtered_data <- filter_one_species(data, species)
  
  Nmix_data_tibble <- get_count_and_effort(filtered_data, grid) 
  
  Nmix_data_list <- list(
    effectif = Nmix_data_tibble %>% select(starts_with("effectif")),
    occurence.cov = Nmix_data_tibble %>% 
      select(all_of(selected_cov)) %>% 
      mutate(intersect = 1) %>% 
      relocate(intersect),
    detection.cov = list(
      transect_length = Nmix_data_tibble %>% select(starts_with("transect_length")),
      session = Nmix_data_tibble %>% select(starts_with("session")))
  )
  
  return(Nmix_data_list)
}


# 
# start_time <- Sys.time()
# Nmix_data_tibble <- get_count_and_effort(filtered_data, grid)
# end_time <- Sys.time()
# execution_time <- end_time - start_time
# print(execution_time)
# 
# start_time <- Sys.time()
# Nmix_data_tibble2 <- get_count_and_effort2(filtered_data, grid)
# end_time <- Sys.time()
# execution_time <- end_time - start_time
# print(execution_time)
# 
# all(Nmix_data_tibble == Nmix_data_tibble2)