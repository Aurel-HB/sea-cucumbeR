# This code create a first approach to sdmTMB model for PPP for one stat square

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
data_position <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_position.rds",
        sep=""))
win.area <- 750

data_used <- data_position[data_position$station == 98,]


ggplot(data_used, aes(long, lat)) +
  geom_point(col = "darkblue", alpha = 0.1) +
  coord_cartesian(expand = FALSE)


res <- 10 # Determines resolution: lower value will increase number of zeroes generated

stat_trait <- data.frame(
  long = c(min(data_used$long),max(data_used$long),
           max(data_used$long),min(data_used$long),min(data_used$long)),
  lat = c(max(data_used$lat),max(data_used$lat),
           min(data_used$lat),min(data_used$lat),max(data_used$lat))
)

# zeros is generated on a grid for this example, but other strategies could be used
zeros <- expand.grid(
  x = seq(min(stat_trait$long), max(stat_trait$long), by = res),
  y = seq(min(stat_trait$lat), max(stat_trait$lat), by = res)
)

#create a data frame with the position
dat <- data.frame(
  x = data_used$long,
  y = data_used$lat
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

all_dat$x <- all_dat$x/1000
all_dat$y <- all_dat$y/1000

#### create the mesh #####
mesh <- make_mesh(
  all_dat,
  xy_cols = c("x", "y"),
  cutoff = 0.01 # min. distance between knots in X-Y units
)
plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$x, dat$y, pch = 19, col = "red",cex = 0.3)


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
  mesh = mesh,
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
  coord_fixed(ratio = 0.01,
              expand = FALSE,
              xlim = c(min(data_used$long), max(data_used$long)),
              ylim = c(min(data_used$lat), max(data_used$lat)))




