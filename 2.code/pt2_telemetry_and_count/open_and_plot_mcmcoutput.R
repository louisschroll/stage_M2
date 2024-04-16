


load("3.results/mcmc_outputs/mcmc_output_puffin_de_scopoli_R.rdata")
load("1.data/all_seabirds_counts.rdata")
load("1.data/grid_cells.rdata")
load("1.data/RSF_data_scopoli_shearwaters.rdata")
# Load package
library(tidyverse, warn.conflicts = FALSE)
library(sf)
library(nimble, warn.conflicts = FALSE)
library(patchwork)

path_to_Rfunc <- "2.code/pt2_telemetry_and_count/R_func"
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

best_cov <- list(
  goeland_leucophee_HR = c(
    "log_dist_to_shore",
    "log_bathymetry",
    "mean_winter_SST",
    "mean_spring_SST",
    "mean_summer_SST",
    "mean_autumn_SST"
  ),
  goeland_leucophee_R = c(
    "log_dist_to_shore",
    "mean_winter_SST",
    "mean_summer_SST",
    "mean_autumn_SST"
  ),
  
  petit_puffin_HR = c(
    "log_dist_to_shore",
    "log_bathymetry",
    'mean_CHL',
    "sd_SAL",
    "mean_SSH"
  ),
  petit_puffin_R = c("mean_CHL", "mean_SSH"),
  
  puffin_de_scopoli_R = c("log_mean_CHL", "sd_SAL", "mean_SSH", "sd_SSH"),
  
  sterne_caugek_HR = c(
    "mean_CHL",
    "mean_SSH",
    "mean_winter_SST",
    "mean_spring_SST",
    "mean_summer_SST"
  ),
  sterne_caugek_R = c(
    "log_bathymetry",
    "mean_CHL",
    "mean_SSH",
    "sd_SSH",
    "log_sd_VEL"
  )
)

selected_cov <- best_cov[["sterne_caugek_R"]]
# check convergence
mcmcplots::traplot(samplesRSF)
mcmcplots::denplot(samplesRSF)

MCMCvis::MCMCsummary(object = samplesRSF, round = 2)

new_grid_RSF <-
  make_prediction(
    samplesRSF,
    grid,
    selected_cov,
    include_intercept = F,
    rsf_intercept = "beta_pop[1]"
  )

plots_rsf <- plot_prediction(new_grid_RSF, add_colonies = F)
#plots$mean_psi_plot + geom_point(data=df_RSF %>% filter(case==1), aes(x=X, y=Y))
plots_rsf$mean_psi_plot + plots_rsf$sd_psi_plot

# check convergence
mcmcplots::traplot(samplesNmixture)
mcmcplots::denplot(samplesNmixture)


plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot

# check convergence
mcmcplots::traplot(samplesint)
mcmcplots::denplot(samplesint)

MCMCvis::MCMCsummary(object = samplesint, round = 2)



plots_int$mean_psi_plot + plots_int$sd_psi_plot

# 1/3 - N-mixture
new_grid_nmix <- make_prediction(samplesNmixture, grid, selected_cov)
plots_nmix <- plot_prediction(new_grid_nmix, plot_title = "N-mixture")
# 2/3 - RSF
new_grid_RSF <- make_prediction(
  samplesRSF,
  grid,
  selected_cov,
  include_intercept = F,
  rsf_intercept = "beta_pop[1]"
)
  
plots_rsf <- plot_prediction(new_grid_RSF, plot_title = "RSF")
plots_rsf_mean <- plots_rsf$mean_psi_plot
# 3/3 - RSF & N-mixture
new_grid_int <- make_prediction(samplesint, 
                                grid, 
                                selected_cov, 
                                rsf_intercept = "beta_pop[1]")
plots_int <- plot_prediction(new_grid_int, 
                             add_colonies = F, 
                             plot_title = "Integrated model")

(plots_rsf_mean + plots_rsf$sd_psi_plot) /
  (plots_nmix$mean_psi_plot + plots_nmix$sd_psi_plot) /
  (plots_int$mean_psi_plot + plots_int$sd_psi_plot) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = '',
    #caption = 'Source: Migralion project & PELMED 2017-2021',
    theme = theme(
      plot.title = element_text(size = 20,
                                face = "bold",
                                hjust = 0.5),
      legend.position = "top"
    )
  )

