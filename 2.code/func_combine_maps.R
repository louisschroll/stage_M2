# HEADER ------------------------------------------------------------------------
#
# Script name:  
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-05-20
#
# Script description:
#
#
# -------------------------------------------------------------------------------


compute_risk_score <- function(file_criteria_value){
  read_excel(file_criteria_value, skip = 1) %>% 
    # Compute collision score and rescale it out of ten
    mutate(
      collision_score = Alt * (Man + Ptf + Noc) / 3,
      collision_summer = collision_score * CS_summer,
      collision_winter = collision_score * CS_winter,
      collision_summer = round(collision_summer / max(collision_summer, na.rm = T) * 10, 2),
      collision_winter = round(collision_winter / max(collision_winter, na.rm = T) * 10, 2)) %>% 
    # Compute displacement score and rescale it out of ten
    mutate(
      displacement_score = (Dis * Spe),
      displacement_summer = displacement_score * CS_summer,
      displacement_winter = displacement_score * CS_winter,
      displacement_summer = round(displacement_summer / max(displacement_summer, na.rm = T) * 10, 2),
      displacement_winter = round(displacement_winter / max(displacement_winter, na.rm = T) * 10, 2)
    ) %>% 
    # Compute vulnerability (displacement + collision) score and rescale it out of ten
    mutate(
      # using the mean btw displacement and collision
      vul_mean_summer = rowMeans(select(., displacement_summer, collision_summer)),
      vul_mean_winter = rowMeans(select(., displacement_winter, collision_winter)),
      # using the max btw displacement and collision
      vul_maxi_summer = pmax(collision_summer, displacement_summer),
      vul_maxi_winter = pmax(collision_winter, displacement_winter),
      # using the formula from Furness et al., 2013 (I keep only this score)
      vul_furn = (Alt + Man + Ptf + Noc) * (Dis + Spe) / 8,
      vul_furn_summer = vul_furn * CS_summer,
      vul_furn_winter = vul_furn * CS_winter,
      vulnerability_summer = vul_furn_summer / max(vul_furn_summer, na.rm = T) * 10,
      vulnerability_winter = vul_furn_winter / max(vul_furn_winter, na.rm = T) * 10
    ) %>% 
    select(nom_fr, 
           collision_summer, collision_winter, 
           displacement_summer, displacement_winter,
           vulnerability_summer, vulnerability_winter) 
}


create_grid_risk <- function(grid_result, score_df){
  grid_result %>%
    select(grid_c) %>% 
    mutate(
      displacement_summer = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "summer",
        risk = "displacement"
      ),
      displacement_winter = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "winter",
        risk = "displacement"
      ),
      collision_summer = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "summer",
        risk = "collision"
      ),
      collision_winter = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "winter",
        risk = "collision"
      ),
      vulnerability_summer = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "summer",
        risk = "vulnerability"
      ),
      vulnerability_winter = compute_risk_col(
        results_grid = grid_result,
        score_df = score_df,
        season = "winter",
        risk = "vulnerability"
      )) %>% 
    mutate(
      displacement_average = rowMeans(select(as_tibble(.), 
                                             displacement_summer, 
                                             displacement_winter)),
      collision_average = rowMeans(select(as_tibble(.), 
                                          collision_summer, 
                                          collision_winter)),
      vulnerability_average = rowMeans(select(as_tibble(.), 
                                              vulnerability_summer, 
                                              vulnerability_winter))
      )
}


compute_risk_col <- function(results_grid,
                             score_df,
                             season = "summer",
                             risk = "vul_furn_summer") {
  if (season == "summer") {
    species_names_end = "_R"
    extract_pattern = "(?<=mean_psi_)(.*?)(?=_R)"
  } else if (season == "winter") {
    species_names_end = "_HR"
    extract_pattern = "(?<=mean_psi_)(.*?)(?=_HR)"
  }
  mean_psi_df <- results_grid %>%
    as_tibble() %>%
    select(starts_with("mean_psi")) %>%
    select(ends_with(species_names_end)) %>%
    mutate(row_id = 1:nrow(.)) %>%
    pivot_longer(-row_id, values_to = "mean_psi", names_to = "species_name") %>%
    mutate(species_name = str_extract(string = species_name, pattern = extract_pattern))
  
  species_list <- unique(mean_psi_df$species_name)
  
  filtered_score <- score_df %>% filter(nom_fr %in% species_list)
  species_order <- filtered_score %>% pull(nom_fr)
  score_matrix <- filtered_score %>% pull(paste0(risk, "_", season))
  
  mean_psi_matrix <- mean_psi_df %>%
    pivot_wider(names_from = species_name, values_from = mean_psi) %>%
    select(-row_id) %>%
    relocate(all_of(species_order)) %>%
    as.matrix()
  
  risk_column <- c(mean_psi_matrix %*% score_matrix)
  
  return(risk_column)
}




plot_vulnerability_map <- function(grid_vulnerability, 
                                   fill_with,
                                   legend_position = "bottom",
                                   plot_title = "",
                                   plot_title_size = 8,
                                   legend_title_size = 9,
                                   axis_text_size = 4,
                                   plot_font = "Arial"){
  ggplot() +
    geom_sf(data = grid_vulnerability, aes(fill = .data[[fill_with]]), color = NA, lwd = 0) +
    scale_fill_distiller(palette = "Spectral", 
                         breaks = c(min(grid_vulnerability[[fill_with]]), max(grid_vulnerability[[fill_with]])),
                         labels = c("Low", "High")) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    geom_sf(data = wind_farm, fill = NA, col = "black") +
    #geom_sf(data = wind_farm, fill = NA, col = "black") +
    labs(title = plot_title) +
    theme_bw() +
    theme(
      #text = element_text(family = plot_font),
      legend.position = legend_position,
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5, 
                                face = "bold", 
                                size = plot_title_size),
      plot.title.position = 'plot',
      axis.text = element_text(size = axis_text_size)
    ) +
    guides(
      fill = guide_colourbar(
        title = "Seabirds vulnerability",
        title.theme = element_text(
          face = "bold",
          size = legend_title_size,
          # family = plot_font,
          hjust = 0.5),
        title.position = 'top',
        ticks = FALSE,
        title.hjust = .5,
        barwidth = unit(6, 'lines'),
        barheight = unit(.4, 'lines'),
        label.theme = element_text(size = legend_title_size * 0.65)),
      colour = "none")
}
