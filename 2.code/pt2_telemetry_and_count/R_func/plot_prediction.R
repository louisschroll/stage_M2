#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/plot_prediction.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-04-03
#'
#' Script description:
#' This function makes two plots of the study area representing the mean intensity
#' of space use and the standard deviation associated with the predicted values, and 
#' return them in a list.
#' Mandatory arguments:
#' @new_grid: a sf object, with at least two columns called mean_psi and sd_psi
#' Facultative arguments:
#' @add_colony: a boolean indicating whether colony stored in countLL file should 
#' be plot on the map
#' @species_colony: if add_colony is true, species_colony is the name of the species
#' for which we want to plot colony
#' -------------------------------------------------------------------------------

load("~/stage_M2/1.data/contour_golfe_du_lion.rdata")
load("~/stage_M2/1.data/countLL.rdata")
wind_farm <- st_read(dsn = "~/stage_M2/1.data/EMODnet_HA_Energy_WindFarms_20231124", layer = "EMODnet_HA_Energy_WindFarms_pg_20231124") %>% 
  st_transform(crs = st_crs(pelmed_obs)) %>% 
  st_crop(st_bbox(grid)) %>% 
  select(geometry)


plot_prediction <- function(new_grid,
                            add_colonies = F,
                            species_colony = NULL,
                            legend_position = "bottom",
                            plot_title = "",
                            plot_title_size = 10,
                            legend_title_size = 9,
                            axis_text_size = 6,
                            plot_font = "Arial") {
    
  # Plot the intensity of space use
  mean_psi_plot <- ggplot() +
    geom_sf(data = new_grid, aes(fill = mean_psi), color = NA, lwd = 0) +
    scale_fill_distiller(palette = "Spectral", 
                         breaks = c(min(new_grid$mean_psi), max(new_grid$mean_psi)),
                         labels = c("Low", "High")) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    geom_sf(data = wind_farm, fill = NA, col = "black") +
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
        title = "Relative space-use",
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
  
  if (add_colonies) {
    mean_psi_plot <- mean_psi_plot +
      geom_sf(data = countLL %>% filter(Species == species_colony) %>% st_crop(st_bbox(contour_golfe)),
              pch = 16, size = 2)
  }

  # Plot uncertainty in the prediction
  sd_psi_plot <- ggplot() +
    geom_sf(data = new_grid, aes(fill = sd_psi), color = NA, lwd = 0) +
    scale_fill_viridis_c(option = "inferno",
                         breaks = c(min(new_grid$sd_psi), max(new_grid$sd_psi)),
                         labels = c("Low", "High")
                         ) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = plot_title) +
    theme_bw() +
    theme(#text = element_text(family = plot_font),
          legend.position = legend_position,
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, 
                                    face = "bold", 
                                    size = plot_title_size),
          axis.text = element_text(size = axis_text_size)) +
    guides(
      fill = guide_colourbar(
        title = "Standard deviation",
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

  return(list(mean_psi_plot = mean_psi_plot, sd_psi_plot = sd_psi_plot))
}


plot_occupancy <- function(grid_occupancy, 
                           add_colonies = F, 
                           species_colony=NULL, 
                           legend_position="bottom", 
                           plot_title = "", 
                           plot_title_size = 10,
                           legend_title_size = 9,
                           axis_text_size = 6,
                           plot_font = "Arial"){
  # Plot the intensity of space use
  mean_psi_plot <- ggplot() +
    geom_sf(data = grid_occupancy, aes(fill = mean_psi), color = NA) +
    scale_fill_viridis_c(breaks = c(min(grid_occupancy$mean_psi), max(grid_occupancy$mean_psi)),
                         labels = c("Low", "High")) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    geom_sf(data = wind_farm, fill = NA, col = "black") +
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
        title = "Presence probability",
        limits = c(0,1),
        title.theme = element_text(
          face = "bold",
          size = legend_title_size,
          #family = plot_font,
          hjust = 0.5),
        title.position = 'top',
        ticks = FALSE,
        title.hjust = .5,
        barwidth = unit(6, 'lines'),
        barheight = unit(.4, 'lines'),
        label.theme = element_text(size = legend_title_size * 0.65)),
      colour = "none")
  
  if (add_colonies) {
    mean_psi_plot <- mean_psi_plot +
      geom_sf(data = countLL %>% filter(Species == species_colony) %>% st_crop(st_bbox(contour_golfe)),
              pch = 16, size = 2)
  }
  
  # Plot uncertainty in the prediction
  sd_psi_plot <- ggplot() +
    geom_sf(data = grid_occupancy, aes(fill = sd_psi), color = NA) +
    scale_fill_viridis_c(option = "inferno",
                         breaks = c(min(grid_occupancy$sd_psi), max(grid_occupancy$sd_psi)),
                         labels = c("Low", "High")
    ) +
    geom_sf(data = contour_golfe, color = "black", fill = "antiquewhite") +
    labs(title = plot_title) +
    theme_bw() +
    theme(
      #text = element_text(family = plot_font),
          legend.position = legend_position,
          legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5, 
                                    face = "bold", 
                                    size = plot_title_size),
          axis.text = element_text(size = axis_text_size)) +
    guides(
      fill = guide_colourbar(
        title = "Standard deviation",
        title.theme = element_text(
          face = "bold",
          size = legend_title_size,
          #family = plot_font,
          hjust = 0.5),
        title.position = 'top',
        ticks = FALSE,
        title.hjust = .5,
        barwidth = unit(6, 'lines'),
        barheight = unit(.4, 'lines'),
        label.theme = element_text(size = legend_title_size * 0.65)),
      colour = "none")
  
  return(list(mean_psi_plot = mean_psi_plot, sd_psi_plot = sd_psi_plot))
}

blue <-  "#3d405b"
red <- "#e07a5f"
purple <-  "#81b29a"
#orange <- "#f4d35e"
orange <- "#8AB17D"
### Coeff
homemade_rename <- function(df, new_names){
  colnames(df) <- new_names
  return(df)
}

#covariate_names <- c("Chlorophyl A", "Salinity", "Sea Surface Height (mean value)", "Sea Surface Height (standard deviation)")
  
select_coeff_cols <- function(mcmc_samples, covariate_names, model_name){
  mcmc_samples %>% map(~{.x %>% 
      as_tibble() %>% 
      janitor::clean_names() %>% 
      select(starts_with("beta")) }) %>% 
    bind_rows() %>% 
    homemade_rename(tolower(covariate_names)) %>% 
    select(-starts_with("intercept")) %>% 
    mutate(model = model_name)
}


gather_coeff_values <- function(sampleNmixture = NULL, samplesRSF = NULL, samplesint = NULL, samples_spOcc = NULL, selected_cov){
  bind_rows(
    samplesNmixture %>%
      select_coeff_cols(
        covariate_names = c("intercept", selected_cov),
        model_name = "N-mixture"
      ),
    
    samplesRSF %>%
      select_coeff_cols(
        covariate_names = c("intercept", selected_cov),
        model_name = "RSF"
      ),
    
    samplesint %>%
      select_coeff_cols(
        covariate_names = c("intercept_Nmix", "intercept_RSF", selected_cov),
        model_name = "Integrated model"
      ),
    
    samples_spOcc %>% as_tibble() %>%
      janitor::clean_names() %>%
      select(-starts_with("intercept")) %>%
      #homemade_rename(new_names = covariate_names) %>% 
      mutate(model = "spOccupancy"),
  ) %>% 
    pivot_longer(-model, names_to = "covariates")
}


plot_coeff_by_model <- function(coeff_df, 
                                covariate_name, 
                                plot_title = "", 
                                plot_subtitle = "",
                                plot_title_size = 10,
                                axis_text_size = 7,
                                plot_font = "Arial"){
  ggplot(coeff_df %>% filter(covariates == covariate_name), aes(value, model, color = model)) +
    ggdist::stat_halfeye(
      adjust = .5,
      width = .6,
      .width = 0,
      justification = -.3,
      point_colour = NA,
      aes(fill = model)
    ) +
    geom_boxplot(width = .25,
                 outlier.shape = NA) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "",
      y = ""
    ) +
    theme_minimal() +
    scale_color_manual(values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "N-mixture" = blue,
      "spOccupancy" = orange
    )) +
    scale_fill_manual(values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "N-mixture" = blue,
      "spOccupancy" = orange
    )) +
    theme(
      #text = element_text(family = plot_font),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_blank(),
      plot.title = element_text(
        size = plot_title_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(hjust = 0.5, 
                                   size = 8),
      legend.position = "bottom",
      legend.title = element_blank()
    ) +
    guides(color = guide_legend(
      label.position = "top",
      keywidth = unit(1, "pt"),
      nrow = 1
    ),
    fill = "none")
}

# coeff_df <- gather_coeff_values(sampleNmixture, samplesRSF, samplesint, samples_spOcc, selected_cov)
# 
# sal_plot <- plot_coeff_by_model(coeff_df, 
#                     covariate_name = "sd_sal", 
#                     plot_title = "Salinity variability effect", 
#                     plot_subtitle = "(standard deviation of salinity)")
# 
# sd_ssh_plot <- plot_coeff_by_model(coeff_df, 
#                     covariate_name = "sd_ssh", 
#                     plot_title = "SSH variability effect", 
#                     plot_subtitle = "(standard deviation of SSH)")
# 
# mean_ssh_plot <- plot_coeff_by_model(coeff_df, 
#                                 covariate_name = "mean_ssh", 
#                                 plot_title = "SSH effect", 
#                                 plot_subtitle = "(mean of SSH)")
# 
# chl_plot <- plot_coeff_by_model(coeff_df, 
#                                 covariate_name = "log_mean_chl", 
#                                 plot_title = "Logarithmic chlorophyl A effect", 
#                                 plot_subtitle = "(mean)")
# 
# ((sd_ssh_plot | sal_plot) / (mean_ssh_plot | chl_plot)) +
#   plot_layout(guides = "collect", tag_level = "new") &
#   theme(
#     legend.position = "bottom",
#     legend.text = element_text(
#       family  = "Helvetica",
#       face = "bold",
#       size = 12
#     )
#   )

plot_space_use_precision <- function(grid_nmix, 
                                     grid_RSF, 
                                     grid_int, 
                                     grid_spOcc,
                                     plot_title_size = 10,
                                     plot_subtitle_size = 7,
                                     axis_text_size = 8,
                                     plot_font = "Arial"){
  tibble(value = grid_nmix$sd_psi, model = "N-mixture") %>%
  bind_rows(tibble(value = grid_RSF$sd_psi, model = "RSF"),
            tibble(value = grid_int$sd_psi, model = "Integrated model"),
            tibble(value = grid_spOcc$sd_psi, model = "spOccupancy")) %>%
  ggplot(aes(x = value, y = model, color = as.factor(model))) +
  ggdist::stat_halfeye(
    adjust = .5,
    width = .6,
    .width = 0,
    justification = -.3,
    point_colour = NA,
    aes(fill = model)
  ) +
  geom_boxplot(width = .25,
               outlier.shape = NA) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Space-use precision",
    subtitle = "Predicted CV among grid-cells",
    x = "Coefficient of variation (CV)",
    y = "") +
  scale_color_manual(
    values = c("RSF" = red, "Integrated model" = purple, "N-mixture" = blue)) +
  scale_fill_manual(
    values = c("RSF" = red, "Integrated model" = purple, "N-mixture" = blue)) +
  theme(
    #text = element_text(family = plot_font),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_text_size),
    plot.title = element_text(
      size = plot_title_size,
      face = "bold",
      hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, 
                                 size = plot_subtitle_size),
    legend.position = "top",
    legend.title = element_blank()) +
  guides(color = "none", fill = "none")
  }



plot_coeff_by_model2 <- function(coeff_df, 
                                 plot_title = "", 
                                 plot_subtitle = "",
                                 plot_title_size = 10,
                                 axis_text_size = 7,
                                 plot_font = "Arial"){
  ggplot(coeff_df, aes(value, model, color = model)) +
    facet_wrap(~covariates, scales = "fixed", ncol = 1) +
    ggdist::stat_halfeye(
      adjust = .5,
      width = .6,
      .width = 0,
      justification = -.3,
      point_colour = NA,
      aes(fill = model)
    ) +
    geom_boxplot(width = .25,
                 outlier.shape = NA) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = "",
      y = ""
    ) +
    theme_minimal() +
    scale_color_manual(values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "N-mixture" = blue,
      "spOccupancy" = orange
    )) +
    scale_fill_manual(values = c(
      "RSF" = red,
      "Integrated model" = purple,
      "N-mixture" = blue,
      "spOccupancy" = orange
    )) +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    theme(
      #text = element_text(family = plot_font),
      axis.text.x = element_text(size = axis_text_size),
      axis.text.y = element_blank(),
      plot.title = element_text(
        size = plot_title_size,
        face = "bold",
        hjust = 0.5
      ),
      plot.subtitle = element_text(hjust = 0.5, size = 8),
      legend.position = "bottom",
      legend.title = element_blank(),
      strip.text = element_text(face = "bold", size = plot_title_size),
      strip.background = element_rect(fill = "white", linetype = "solid",
                                      color = "black", linewidth = 0.9),
      panel.background = element_rect(fill = "transparent", color = "black", linewidth = 0.9)
    ) +
    guides(
      color = guide_legend(
        label.position = "top",
        keywidth = unit(1, "pt"),
        nrow = 1
      ),
      fill = "none"
    )
  
}
