
# install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE, type = "binary")
# devtools::install_github(repo = "https://github.com/Scott-Foster/RISDM", build_vignettes = FALSE)

# load necessary package
library(INLA)
library(RISDM)
library(terra)

# Get covars and simulated data
covars <- rast(system.file("extdata", "ACT_DemoData.grd", package = "RISDM"))
names(covars) <- c("lACC", "SMRZ", "TEMP")

dat <- simulateData.isdm(rasterCovars = covars[[c("SMRZ", "TEMP")]],
                         rasterBiasCovar = covars[["lACC"]],
                         control = list(doPlot = F, set.random.seed=T, random.seed = 123456))

# Create a mesh
meshy <- makeMesh(dat$covarBrick$lACC, max.n = c(1000, 500),
                  dep.range = 0.5,
                  offset = 10,
                  expans.mult = 7.5,
                  doPlot = F)

checkMesh(meshy)

# Fit model
fm <- isdm(observationList = list(POdat = dat$PO,
                                  AAdat = dat$AA,
                                  PAdat = dat$PA),
           covars = dat$covarBrick,
           mesh = meshy,
           responseNames = c(PA="PA", AA="AA"),
           sampleAreaNames = c(PO=NULL, PA="transectArea", AA="transectArea"),
           distributionFormula = ~0+SMRZ+TEMP,
           biasFormula = ~1+lACC,
           artefactFormulas = list(PA=~1, AA=~1),
           control = list(prior.range = c(0.5, 0.1),
                          prior.space.sigma = c(2, 0.1),
                          coord.names = c("x", "y")))


summary(fm)

# check residuals
plot(fm, nFigRow = 3, ask = F)

# Prediction
fm$preds <- predict(object = fm, covars = dat$covarBrick, S = 1000)

tmprast <- fm$preds$field[[c(2,1,3)]]

plot(tmprast, 
     range = range(values(tmprast), na.rm=TRUE),
     col = hcl.colors(25, "viridis", rev = TRUE),
     nc = 3)
df <- as.data.frame(covars) %>% 
  bind_cols(crds(covars))

ggplot() +
  geom_raster(data = df, aes(x = x, y = y, fill = TEMP)) +
  scale_fill_distiller(palette = "Spectral", name = "T°") +
  geom_point(data = dat$AA, aes(x = x, y = y, size = AA))
  


