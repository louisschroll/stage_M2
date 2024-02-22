
library(spOccupancy)
library(coda)
library(stars)
library(ggplot2)
library(tidyverse)
set.seed(102)

data("hbef2015")
str(hbef2015)

sp.names = dimnames(hbef2015$y)[[1]]
ovenHBEF = hbef2015
ovenHBEF$y = ovenHBEF$y[sp.names == "OVEN", , ]
table(ovenHBF$y)

oven.occ.formula <- ~ scale(Elevation) + I(scale(Elevation)^2)
oven.det.formula <- ~ scale(day) + scale(tod) + I(scale(day)^2)
# Check out the format of ovenHBEF
str(ovenHBEF)

# Format with explicit specification of inits for alpha and beta
# with four detection parameters and three occurrence parameters 
# (including the intercept).
oven.inits <- list(alpha = c(0, 0, 0, 0), 
                   beta = c(0, 0, 0), 
                   z = apply(ovenHBEF$y, 1, max, na.rm = TRUE))

oven.priors <- list(alpha.normal = list(mean = 0, var = 2.72), 
                    beta.normal = list(mean = 0, var = 2.72))

n.samples <- 5000
n.burn <- 3000
n.thin <- 2
n.chains <- 3

out <- PGOcc(occ.formula = oven.occ.formula, 
             det.formula = oven.det.formula, 
             data = ovenHBEF, 
             inits = oven.inits, 
             n.samples = n.samples, 
             priors = oven.priors, 
             n.omp.threads = 1, 
             verbose = TRUE, # rport progress
             n.report = 1000, # every 1000 iterations
             n.burn = n.burn, 
             n.thin = n.thin, 
             n.chains = n.chains)

summary(out)
plogis(2.13)

plot(out$beta.samples)
plot(out$alpha.samples)

ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
summary(ppc.out)


ppc.df <- data.frame(fit = ppc.out$fit.y, 
                     fit.rep = ppc.out$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')


diff.fit <- ppc.out$fit.y.rep.group.quants[3, ] - ppc.out$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')

waicOcc(out)
#  elpd is the expected log point-wise predictive density
# PD is the effective number of parameters

# Prediction

data(hbefElev)
str(hbefElev)

elev.pred <- (hbefElev$val - mean(ovenHBEF$occ.covs[, 1])) / sd(ovenHBEF$occ.covs[, 1])
# These are the new intercept and covariate data.
X.0 <- cbind(1, elev.pred, elev.pred^2) 
out.pred <- predict(out, X.0)

plot.dat <- data.frame(x = hbefElev$Easting, 
                       y = hbefElev$Northing, 
                       mean.psi = apply(out.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.pred$psi.0.samples, 2, sd), 
                       stringsAsFactors = FALSE)
# Make a species distribution map showing the point estimates,
# or predictions (posterior means)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') +
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'Mean OVEN occurrence probability') +
  theme_bw()



#### Integrated model
library(spOccupancy)
library(coda)
library(stars)
library(ggplot2)
library(tidyverse)
set.seed(102)

data(hbef2015)
data(neon2015)

sp.names <- dimnames(hbef2015$y)[[1]]
ovenHBEF <- hbef2015
ovenHBEF$y <- ovenHBEF$y[sp.names == "OVEN", , ]
ovenNEON <- neon2015
ovenNEON$y <- ovenNEON$y[sp.names == "OVEN", , ]
table(ovenHBEF$y)
table(ovenNEON$y)
str(ovenHBEF)
occ.covs.int <- rbind(ovenHBEF$occ.covs, ovenNEON$occ.covs)
str(occ.covs.int)

sites.int <- list(hbef = 1:nrow(ovenHBEF$occ.covs), 
                  neon = 1:nrow(ovenNEON$occ.covs) + nrow(ovenHBEF$occ.covs))
str(sites.int)

y.int <- list(hbef = ovenHBEF$y, 
              neon = ovenNEON$y)
str(y.int)

det.covs.int <- list(hbef = ovenHBEF$det.covs, 
                     neon = ovenNEON$det.covs)

data.int <- list(y = y.int, 
                 occ.covs = occ.covs.int, 
                 det.covs = det.covs.int, 
                 sites = sites.int)
str(data.int)

occ.formula.int <- ~ scale(Elevation) + I(scale(Elevation)^2)
det.formula.int = list(hbef = ~ scale(day) + scale(tod) + I(scale(day)^2), 
                       neon = ~ scale(day) + scale(tod) + I(scale(day)^2))

# Total number of sites
J <- nrow(data.int$occ.covs)
inits.list <- list(alpha = list(0, 0),
                   beta = 0, 
                   z = rep(1, J))

prior.list <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = list(0, 0), 
                                       var = list(2.72, 2.72)))

n.samples <- 8000
n.burn <- 3000
n.thin <- 25

# Approx. run time: < 15 sec
out.int <- intPGOcc(occ.formula = occ.formula.int,
                    det.formula = det.formula.int, 
                    data = data.int,
                    inits = inits.list,
                    n.samples = n.samples, 
                    priors = prior.list, 
                    n.omp.threads = 1, 
                    verbose = TRUE, 
                    n.report = 2000, 
                    n.burn = n.burn, 
                    n.thin = n.thin, 
                    n.chains = 3) 

summary(out.int)

ppc.int.out <- ppcOcc(out.int, 'freeman-tukey', group = 2)
summary(ppc.int.out)

waicOcc(out.int)

out.int.k.fold <- intPGOcc(occ.formula = occ.formula.int,
                           det.formula = det.formula.int, 
                           data = data.int,
                           inits = inits.list,
                           n.samples = n.samples, 
                           priors = prior.list, 
                           n.omp.threads = 1, 
                           verbose = FALSE, 
                           n.report = 2000, 
                           n.burn = n.burn, 
                           n.thin = n.thin, 
                           n.chains = 1,
                           k.fold = 4) 
out.int.k.fold.small <- intPGOcc(occ.formula = ~ 1, 
                                 det.formula = list(hbef = ~ 1, neon = ~ 1), 
                                 data = data.int,
                                 inits = inits.list,
                                 n.samples = n.samples, 
                                 priors = prior.list, 
                                 n.omp.threads = 1, 
                                 verbose = FALSE, 
                                 n.report = 2000, 
                                 n.burn = n.burn, 
                                 n.thin = n.thin, 
                                 n.chains = 1,
                                 k.fold = 4) 
# Summarize the CV results
out.int.k.fold$k.fold.deviance
out.int.k.fold.small$k.fold.deviance


out.int.k.fold.hbef <- intPGOcc(occ.formula = occ.formula.int,
                                det.formula = det.formula.int, 
                                data = data.int,
                                inits = inits.list,
                                n.samples = n.samples, 
                                priors = prior.list, 
                                n.omp.threads = 1, 
                                verbose = TRUE, 
                                n.report = 2000, 
                                n.burn = n.burn, 
                                n.thin = n.thin, 
                                k.fold = 4, 
                                k.fold.data = 1) 

# Look at CV results again 
# Single data source model
out.k.fold$k.fold.deviance
# Integrated model
out.int.k.fold.hbef$k.fold.deviance


# Make sure to standardize using mean and sd from fitted model
elev.pred <- (hbefElev$val - mean(data.int$occ.covs[, 1])) / sd(data.int$occ.covs[, 1])
X.0 <- cbind(1, elev.pred, elev.pred^2)
out.int.pred <- predict(out.int, X.0)

# Producing an SDM for HBEF alone (posterior mean)
plot.dat <- data.frame(x = hbefElev$Easting, 
                       y = hbefElev$Northing, 
                       mean.psi = apply(out.int.pred$psi.0.samples, 2, mean), 
                       sd.psi = apply(out.int.pred$psi.0.samples, 2, sd))

dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') + 
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'Mean OVEN occurrence probability using intPGOcc') +
  theme_bw()

ggplot() + 
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = sd.psi)) +
  scale_fill_viridis_c(na.value = 'transparent') + 
  labs(x = 'Easting', y = 'Northing', fill = '', 
       title = 'SD OVEN occurrence probability using intPGOcc') +
  theme_bw()

