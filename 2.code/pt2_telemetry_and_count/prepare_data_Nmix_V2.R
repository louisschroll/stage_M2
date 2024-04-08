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


prepare_data_Nmix_V2 <- function(data_list, grid, species, selected_cov){
  ndataset <- length(data_list)
  if (ndataset == 1){
    data_nmix <- filter_one_species(data_list[[1]], species) %>% 
      get_count_and_effort(grid) %>% 
      prepare_data_for_1_dataset(selected_cov)
  } else {
    data_nmix <- map(data_list, 
                     ~{filter_one_species(.x, species) %>% 
                         get_count_and_effort(grid)}) %>% 
      prepare_data_for_several_datasets(selected_cov)
  }
  return(data_nmix)
}


filter_one_species <- function(data, species){
  obs <- data$obs %>% filter(species_name == species)
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
      summarise(transectLength = sum(len))
    
    # Retrieve effectif 
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
      mutate(effectif = ifelse(is.na(effectif) & !is.na(transectLength), 0, effectif),
             transectLength = ifelse(is.na(transectLength), 0, transectLength)) %>%
      rename_with(~ paste0(.x, "_", as.character(sessionK), recycle0 = TRUE), 
                  .cols = c(effectif, transectLength))
  }
  return(sampled_gridcells)
}


prepare_data_for_1_dataset <- function(nmix_tibble, selected_cov){
  n.occ.cov <- length(selected_cov) + 1
  n.det.cov <- 2
  # Get count data
  effectif_df <- get_effectif(nmix_tibble)
  data <- list(nobs = effectif_df)
  
  # Set initial values
  initial.values <- list(beta = rnorm(n.occ.cov, 0, 1), 
                         alpha = rnorm(n.det.cov, 0, 1), 
                         N = apply(effectif_df, 1, sum) + 1)
  
  # Get constants
  constants <- list(XN = get_occ_cov(nmix_tibble),
                    transect_length = get_and_scale_transect_length(nmix_tibble),
                    nsites = nrow(effectif_df),
                    nreplicates = ncol(effectif_df),
                    n.occ.cov = n.occ.cov)
  
  return(list(data = data, constants = constants, inits = initial.values))
}


prepare_data_for_several_datasets <- function(nmix_tibble_list, selected_cov){
  n.occ.cov <- length(selected_cov) + 1
  n.det.cov <- 2
  result_list <- map(nmix_tibble_list, ~ {
    .x %>%
      select(starts_with("transectLength"), starts_with("effectif"), id) %>%
      relocate(id) %>%
      pivot_longer(-c(id), names_to = c(".value", "session"), names_pattern = "(.+)_(.+)") %>%
      filter(!is.na(effectif)) %>% 
      mutate(transectLength = (transectLength - mean(transectLength))/sd(transectLength))
  }) %>% 
    bind_rows(.id = "dataset_nb") %>% 
    mutate(dataset_nb = as.numeric(as.factor(dataset_nb)))
  
  # Get the observation data
  data <- list(nobs = result_list$effectif)
  
  # Set initial values
  N0 <- result_list %>% group_by(id) %>% 
    summarise(N0 = sum(effectif)+1) %>% 
    arrange(by = "id") %>% 
    pull(N0)
  
  ndatasets <- max(result_list$dataset_nb)
  
  initial.values <- list(beta = rnorm(n.occ.cov,0,1), 
                        alpha = matrix(rnorm(n.det.cov*ndatasets), nrow = n.det.cov),
                        N = N0)
  
  # Get the constants
  sampled_sites <- map(nmix_tibble_list, function(x) select(x, all_of(c(selected_cov, "id")))) %>% 
    bind_rows() %>% 
    unique() %>% 
    arrange(by = id) %>% 
    get_occ_cov()
  
  constants <- list(XN = sampled_sites,
                    dataset_nb = result_list$dataset_nb,
                    transect_length = result_list$transectLength,
                    site_id = get_new_sites_id(result_list),
                    n.occ.cov = n.occ.cov,
                    nsites_total = nrow(sampled_sites),
                    nsampled_points = nrow(result_list),
                    ndatasets = ndatasets)
  
  return(list(data = data, constants = constants, inits = initial.values))
}


get_effectif <- function(nmix_tibble){
  nmix_tibble %>% select(starts_with("effectif"))
}


get_occ_cov <- function(nmix_tibble){
  nmix_tibble %>% 
    select(all_of(selected_cov)) %>% 
    mutate(intersect = 1) %>% 
    relocate(intersect) %>% 
    as.matrix()
}


get_and_scale_transect_length <- function(nmix_tibble){
  nmix_tibble %>% 
    select(starts_with("transect_length")) %>% 
    mutate(row_nb = 1:nrow(.)) %>% 
    pivot_longer(-row_nb) %>%
    # Scale
    mutate(value = (value - mean(value))/sd(value)) %>%
    pivot_wider(names_from = name, values_from = value) %>%
    select(-row_nb) %>% 
    as.matrix()
}


get_new_sites_id <- function(result_list){
  sampled_sites_id <- unique(result_list$id) %>% sort()
  id_correspondance <- tibble(old_id = sampled_sites_id, new_id = 1:length(sampled_sites_id))
  
  new_sites_id <- result_list %>% select(id) %>% 
    left_join(id_correspondance, by = c("id" = "old_id")) %>%
    pull(new_id) 
  return(new_sites_id)
}


