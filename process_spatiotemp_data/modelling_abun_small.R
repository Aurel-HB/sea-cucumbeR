# This script create a model spatio-temporal to model small sea cucumber abundance


#############
#load packages
#############
library(dplyr) 
library(ggplot2) 
library(sf)
library(sdmTMB)
library(INLA)
library(viridis)
library(DHARMa)
library(here)
library(tweedie)
library(future)

#############################################################
#PART 1 - Load data and prepare mesh and prediction grid
#############################################################

# load the study area ###
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

# load the data ###
data_holotv <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_abun_tot_small_cov.rds",
        sep = "")
)

data_holotv <- data_holotv %>%
  mutate(X=X/1000) %>%
  mutate(Y=Y/1000) # the lat and long are in km for modelling


# prepare the grid for prediction ###
grid_tot <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_tot.rds",
        sep = "")
)

# Create the mesh for inference with INLA ###
#using the grid of the area
bnd <- INLA::inla.nonconvex.hull(cbind(grid_tot$long, grid_tot$lat),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(grid_tot$long, grid_tot$lat),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(as.data.frame(grid_tot)[,c(1,2)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic
mesh <- make_mesh(as.data.frame(data_holotv), c("X", "Y"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)


#############################################################
#PART 2 - Modelling the abundance of small individual
#############################################################

# model 1 without covariate and tweedie distribution
data <- as.data.frame(data_holotv)[,seq(1,ncol(data_holotv)-1,1)]

fit_tot <- sdmTMB(
  intensity ~ 1,#<< fixed intercept ignoring time
  data = data,
  mesh = mesh,
  family = tweedie(link = "log"),
  spatial = "on",
  time = "year",
  spatiotemporal = "ar1",#<< setting an AR(1) spatiotemporal process
  extra_time = c(2024,2025),#<< our list of extra years to be included
)

#summary model
summary(fit_tot)

#sanity model
sanity(fit_tot)

#check residual
set.seed(seed = 14)
model_sim <- simulate(fit_tot, nsim = 1000, type = "mle-mvn")
simulationOutput <- dharma_residuals(model_sim, fit_tot, return_DHARMa = TRUE)
plotQQunif(simulationOutput, testUniformity = FALSE, testOutliers = FALSE,
           testDispersion = FALSE)

data$resids <- residuals(fit_tot)
ggplot(data, aes(X, Y, col = resids)) + scale_colour_gradient2() +
  geom_point(size=3) + facet_wrap(~year) + coord_fixed() +theme(aspect.ratio = 3)


# Prediction ####
grid <- as.data.frame(grid_tot)[,seq(1,ncol(grid_tot)-1,1)]
names(grid)[c(1,6,7)] <- c("temp","X","Y")
predictions <- predict(fit_tot, newdata = grid,
                       return_tmb_object = TRUE)

# small function to make maps
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year, nrow = 1) +
    coord_fixed()+
    theme(aspect.ratio = 3)
}

# total prediction (random and fix effects)
plot_map(predictions$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt") +
  ggtitle("Prediction (fixed effects + all random effects)")

# spatial random effect
plot_map(predictions$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

# spatio-temproral random effect
plot_map(predictions$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()


#  total abundance
index <- get_index(predictions, area = 250000, bias_correct = TRUE)
# be careful our prediction is in nb/m² so area have to be in m²

ggplot(index, aes(year, est*0.4/1000)) + geom_line() +
  #geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.4) +
  xlab('Year') + ylab('Abundance estimate (nombre)')
