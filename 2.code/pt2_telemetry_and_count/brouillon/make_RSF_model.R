# HEADER ------------------------------------------------------------------------
#
# Script name:  ~/stage_M2/2.code/pt2_telemetry_and_count/make_RSF_model.R
# Author:       Louis Schroll
# Email:        louis.schroll@ens-lyon.fr
# Date:         2024-03-28
#
# Script description:
#
#
# -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

library(nimble)
library(tidyverse)

load("1.data/RSF_data_sandwich_tern.rdata")
df_RSF <- df_RSF %>% filter(individual_id %in% 1:3)
selected_cov <- c("mean_SST")
# Run
n.iter = 15000
n.burnin = 0.1 * n.iter
n.chains = 2

data_rsf <- format_rsf_data_for_nimble(df_RSF, selected_cov)
samplesRSF <- run_RSF(data_rsf, n.iter = n.iter, n.burnin=n.burnin, n.chains=n.chains)

# check convergence
mcmcplots::traplot(samplesRSF)
mcmcplots::denplot(samplesRSF)
coda::effectiveSize(samplesRSF)

# make prediction maps
new_grid <- make_prediction(mcmc.output = samplesRSF, 
                            grid = grid, 
                            selected_cov = selected_cov, 
                            include_intercept=T)

plot <- plot_prediction(new_grid)
plot$mean_psi_plot + geom_point(data = df_RSF %>% filter(case==1), aes(x=X, y=Y), alpha = 0.1)
plot$sd_psi_plot
