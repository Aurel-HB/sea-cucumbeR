# This code create a first approach to sdmTMB model for PPP

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(tweedie)
library(ggspatial)
library(sdmTMB)
library(INLA)

shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))
predata_intensity <- readRDS(
  paste(here(), "/point_process/Output/holo_simu_intensity.rds", sep=""))
data_abun <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_abun.rds", sep=""))
data_position <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_position.rds",
                     sep=""))
win.area <- 750

ggplot(data_position, aes(long, lat)) +
  geom_point(col = "darkblue", alpha = 0.1) +
  coord_cartesian(expand = FALSE)


res <- 100 # Determines resolution: lower value will increase number of zeroes generated

# zeros is generated on a grid for this example, but other strategies could be used
zeros <- expand.grid(
  x = seq(min(data_position$long), max(data_position$long), by = res),
  y = seq(min(data_position$lat), max(data_position$lat), by = res)
)

#create a data frame with the position
dat <- data.frame(
  x = data_position$long,
  y = data_position$lat
)

dat$present <- 1
zeros$present <- 0
all_dat <- rbind(dat, zeros)

all_dat$fpres <- as.factor(all_dat$present)
ggplot() +
  geom_point(
    data = all_dat,
    aes(x = x, y = y, col = fpres), size = 0.1, alpha = 0.3
  ) +
  coord_equal() +
  guides(col = guide_legend(title = "Present"))

#### create the mesh ####
# working in UTM km instead of UTM m by divided by 1000 to help creation mesh
all_dat$x <- all_dat$x/1000
all_dat$y <- all_dat$y/1000

mesh <- make_mesh(
  all_dat,
  xy_cols = c("x", "y"),
  cutoff = 1 # min. distance between knots in X-Y units
)

# create mesh with inla without zero ##
bnd <- INLA::inla.nonconvex.hull(cbind(dat$x, dat$y),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(dat$x, dat$y),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(data.frame(dat$x, dat$y)),
  boundary = list(bnd,bnd2),
  cutoff = 100, # minimum triangle edge length
  max.edge = c(1852, 18520), # inner and outer max triangle lengths
)
mesh2 <- make_mesh(as.data.frame(dat), c("x", "y"), mesh = mesh_inla)

plot(mesh2$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$x, dat$y, pch = 19, col = "red",cex = 0.3)

# create mesh with inla with pseudo-zero ##
bnd <- INLA::inla.nonconvex.hull(cbind(all_dat$x, all_dat$y),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(all_dat$x, all_dat$y),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(data.frame(all_dat['x'], all_dat['y'])),
  boundary = list(bnd,bnd2),
  cutoff = 100, # minimum triangle edge length
  max.edge = c(1852, 18520), # inner and outer max triangle lengths
)
mesh2 <- make_mesh(as.data.frame(all_dat), c("x", "y"), mesh = mesh_inla)

plot(mesh2$mesh, main = NA, edge.color = "grey60", asp = 1)
points(all_dat$x, all_dat$y, pch = 19, col = "red",cex = 0.3)



##### Downweighted Poisson Regression (DWPR) #####

## first calculate the weights ##
# small values at presence locations
all_dat$wt <- 1e-6

# pseudo-absences: area per quadrature point
tot_area <- diff(range(dat$x)) * diff(range(dat$y))
n_zeros <- length(which(all_dat$present == 0))

all_dat$wt <- ifelse(all_dat$present == 1,
                     1e-6, tot_area / n_zeros
)

## Then fit the model with the new weights and a Poisson distribution ##
fit <- sdmTMB(
  present / wt ~ 1,
  data = all_dat,
  mesh = mesh2,
  family = poisson(link = "log"),
  weights = all_dat$wt
)

## check the result ##
summary(fit)
sanity(fit)
# plot the random spatial effects
p <- predict(fit, newdata = zeros)
ggplot(p, aes(x, y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE)
# predict spatial distribution both in link (log) space:
ggplot(p, aes(x, y, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
# Or natural space:
ggplot(p, aes(x, y, fill = exp(est))) +
  geom_raster() +
  labs(fill = "Intensity\n(average point density)") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)


##### Poisson point process by log-linked Poisson GLMM #####

## Model 1 with spatial components

fit_ppp <- sdmTMB(
  present ~ 1,
  data = dat,
  mesh = mesh2,
  family = poisson(link = "log"),
  spatial = "on",
  share_range = FALSE)


## check the result ##
summary(fit_ppp)
sanity(fit_ppp)
# plot the random spatial effects
p <- predict(fit_ppp, newdata = zeros)
ggplot(p, aes(x, y, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE)
# predict spatial distribution both in link (log) space:
ggplot(p, aes(x, y, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
# Or natural space:
ggplot(p, aes(x, y, fill = exp(est))) +
  geom_raster() +
  labs(fill = "Intensity\n(average point density)") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)
