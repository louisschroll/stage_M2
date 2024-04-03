
# Load packages
library(dagR)
library(tidyverse)
library(ggplot2)
library(spatialEco)
library(raster)
library(sf)
library(amt)
library(here)
library(localGibbs) # Install from Theo Michelot's GitHub
library(tidyverse)
library(stars)

# Define coolors palette https://coolors.co/5f0f40-9a031e-fb8b24-e36414-0f4c5c
blue <- "#3d405b"
orange <- "#f2cc8f"
red <- "#e07a5f"
purple <- "#81b29a"

# Create landscape with 1 continuous covariate and 1 discrete covariate with spatial autocorrelation
# Landscape
covariate_list <- list()
covariate_list[[1]] <- simRaster(rho = 10, lim = c(-500, 500, -500, 500), res = 1)
covariate_list[[2]] <- simRaster(rho = 30, lim = c(-500, 500, -500, 500), res = 1)

raster1 <- st_as_stars(covariate_list[[1]])
raster2 <- st_as_stars(covariate_list[[2]])

# Discretize covariate 1
raster1$layer <- as_factor(ifelse(raster1$layer > 0.50, 1, 0))
# Scale covariate 2
raster2_scaled <- st_as_stars(scale(covariate_list[[2]]))

# Plot the landscape      
if (FALSE) {
  theme_set(theme_minimal())
  plot1 <- ggplot() +
    geom_stars(data = raster1) + 
    scale_fill_viridis_d(name = "Discrete", labels = c("A", "B")) +
    theme(legend.position = "top")
  
  plot2 <- ggplot() +
    geom_stars(data = raster2) +
    labs(fill = "Continuous") +
    theme(legend.position = "top")
  
  library(patchwork)
  plot1 + plot2
}

# Simulation parameters for biased Correlated Random Walk towards colony
n_individual <- 3    # Number of individuals
n_steps <- 300       # Number of steps per individual
preference <- 4      # Preference of habitat 2 over habitat 1
s0_param <- 0.8      # Parameter 1 for biased-CRW
k_param <- 0.4       # Parameter 2 for biased-CRW
n_tmp <- 10          # Number of potential locations at each step

# Simulation of tracks with a biased CRW towards colony
data_total <- data.frame(matrix(NA, nrow = 0, ncol = 10))
colnames(data_total) <- c('id', 'time', 'x', 'y', 'dist2col', 'step', 'angle', 'direction', 'veget', 'cov')

# Add individual heterogeneity
coef <- rnorm(n_individual, log(preference), 0.25) # Individual heterogeneity 

# Function to simulate tracks
for (j in 1:n_individual) {
  kcpt <- coef[j]
  
  data <- data.frame(matrix(NA, nrow = n_steps, ncol = 10))
  colnames(data) <- c('id', 'time', 'x', 'y', 'dist2col', 'step', 'angle', 'direction', 'veget', 'cov')
  data$id <- j
  data[1, 'time'] <- 0
  data[1, 'x'] <- 0
  data[1, 'y'] <- 0
  data[1, 'dist2col'] <- 0
  
  # First step from the colony
  tmp <- data.frame(matrix(NA, nrow = n_tmp, ncol = 7))
  colnames(tmp) <- c('new.coord.x', 'new.coord.y', 'veg.tmp', 'cov.tmp', 'p.tmp', 'step', 'angle')
  
  for (i in 1:n_tmp) {
    tmp$step[i] <- rgamma(n = 1, shape = 5, scale = 2)
    tmp$angle[i] <- runif(n = 1, 0, 2 * pi)
    tmp[i, 1:2] <- anglePoint(c(0, 0), angl = tmp$angle[i] + pi, len = tmp$step[i])
    tmp$veg.tmp[i] <- st_extract(raster1, cbind(tmp[i, 1], tmp[i, 2]))
    tmp$cov.tmp[i] <- st_extract(raster2_scaled, cbind(tmp[i, 1], tmp[i, 2]))
    tmp$p.tmp[i] <- exp(kcpt * (as.numeric(tmp$veg.tmp[i][[1]]) - 1))
  }
  
  loc_f <- which(rmultinom(1, 1, tmp$p.tmp) == 1)[1]
  
  data[2, 'time'] <- 1
  data[2, 'x'] <- tmp$new.coord.x[loc_f]
  data[2, 'y'] <- tmp$new.coord.y[loc_f]
  data[2, 'dist2col'] <- sqrt(tmp$new.coord.x[loc_f]^2 + tmp$new.coord.y[loc_f]^2)
  data[2, 'step'] <- tmp$step[loc_f]
  data[2, 'angle'] <- tmp$angle[loc_f]
  data[2, 'direction'] <- tmp$angle[loc_f]
  data[2, 'veget'] <- tmp$veg.tmp[loc_f][[1]]
  data[2, 'cov'] <- tmp$cov.tmp[loc_f][[1]]
  
  for (i in 3:n_steps) {
    tmp <- data.frame(matrix(NA, nrow = n_tmp, ncol = 9))
    colnames(tmp) <- c('new.coord.x', 'new.coord.y', 'veg.tmp', 'cov.tmp', 'p.tmp', 'step', 'angle', 'direction', 'dist2col')
    
    for (n in 1:n_tmp) {
      tmp$step[n] <- rgamma(n = 1, shape = 5, scale = 2)
      alpha <- atan2(data[i - 2, 'y'] - 0, data[i - 2, 'x'] - 0) - atan2(data[i - 2, 'y'] - data[i - 1, 'y'], data[i - 2, 'x'] - data[i - 1, 'x'])
      alpha <- 2 * ((abs(alpha) - 0) / (pi - 0)) - 1
      angle <- rnorm(n = 1, mean = 0, sd = s0_param * (1 + k_param * alpha))
      angle <- ifelse(angle > pi, angle - pi, ifelse(angle < (-pi), angle + pi, angle))
      tmp$angle[n] <- angle
      tmp$direction[n] <- data[i - 1, 'direction'] + angle
      tmp$direction[n] <- ifelse(tmp$direction[n] > 2 * pi, tmp$direction[n] - 2 * pi, ifelse(tmp$direction[n] < 0, tmp$direction[n] + 2 * pi, tmp$direction[n]))
      tmp[n, 1:2] <- anglePoint(c(data[i - 1, 'x'], data[i - 1, 'y']), angl = tmp$direction[n] + pi, len = tmp$step[n])
      tmp$dist2col[n] <- sqrt(tmp[n, 1]^2 + tmp[n, 2]^2)
      tmp$veg.tmp[n] <- st_extract(raster1, cbind(tmp[n, 1], tmp[n, 2]))
      tmp$cov.tmp[n] <- st_extract(raster2_scaled, cbind(tmp[n, 1], tmp[n, 2]))
      tmp$p.tmp[n] <- exp(kcpt * (as.numeric(tmp$veg.tmp[n][[1]]) - 1))
    }
    
    if (TRUE %in% is.na(tmp$p.tmp)) {
      loc_f <- which(is.na(tmp$p.tmp) == F)[1]
    } else {
      loc_f <- which(rmultinom(1, 1, na.omit(tmp$p.tmp)) == 1)[1]
    }
    
    data[i, 'step'] <- tmp$step[loc_f]
    data[i, 'angle'] <- tmp$angle[loc_f]
    data[i, 'direction'] <- tmp$direction[loc_f]
    data[i, 'x'] <- tmp[loc_f, 1]
    data[i, 'y'] <- tmp[loc_f, 2]
    data[i, 'dist2col'] <- tmp$dist2col[loc_f]
    data[i, 'time'] <- i - 1
    data[i, 'veget'] <- tmp$veg.tmp[loc_f][[1]]
    data[i, 'cov'] <- tmp$cov.tmp[loc_f][[1]]
  }
  
  data_total <- rbind(data_total, data)
}

head(data_total)
dim(data_total)
data_total$id <- factor(data_total$id)

# Plot the track
if (FALSE) {
  plot_track <- ggplot() +
    geom_stars(data = raster1, alpha = 0.6) +
    scale_fill_viridis_d(name = "Discrete covariate") +
    geom_path(data = data_total, aes(x, y, group = id, color = id), linewidth = 1.2) +
    theme_classic() +
    theme(legend.position = 'top') +
    guides(color = "none") +
    labs(subtitle = paste0("Preference for habitat = 1: p.sel = ", preference))
  
  plot_track
}

# Calculate estimates
preference_estimate <- table(data_total$veget)[1] / table(data_total$veget)[2]
-log(preference_estimate)
log(preference)
