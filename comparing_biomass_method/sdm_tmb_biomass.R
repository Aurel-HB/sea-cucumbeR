# This code create a first approach to sdmTMB model for biomass with abunx0.4

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(tweedie)
library(ggspatial)
library(sdmTMB)
library(INLA)
library(fitdistrplus)


#########################
#test with simulated data
#########################
shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))
predata_intensity <- readRDS(
  paste(here(), "/point_process/Output/holo_simu_intensity.rds", sep=""))
data_abun <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_abun.rds", sep=""))

mean_weight <- 0.4
win.area <- (1.5/1000)*0.5
surfsquare <- 1.852*1.852 # in kilometer square

## Prepare the dataframe
dat <- as.data.frame(
  cbind(data_abun,as.data.frame(predata_intensity[,c("X","Y")])[,1:2]))

# working in UTM km instead of UTM m by divided by 1000 to help creation mesh
dat$long <- dat$X/1000
dat$lat <- dat$Y/1000

# calculate biomass
dat$biomass <- dat$abun*mean_weight

# calculate density
dat$density <- dat$abun/win.area

# prepare surface sampled for offset
dat$surf <- win.area

## Prepare grid for prediction
res <- 1.852 # Determines resolution

# generate a grid
grid <- expand.grid(
  long = seq(min(dat$long), max(dat$long), by = res),
  lat = seq(min(dat$lat), max(dat$lat), by = res)
)

# use the sampling grid to generate a extrapolation grid
grid <- shp_grid %>% 
  dplyr::filter(Id > 85) %>%
  dplyr::filter(Id != 188) %>%
  dplyr::filter(Id != 189) %>%
  st_centroid() %>%
  mutate(long = st_coordinates(.)[,1]/1000) %>%
  mutate(lat = st_coordinates(.)[,2]/1000) %>%
  as.data.frame() %>%
  select(long,lat) 

plot(grid)

## Create mesh with inla
#bnd <- INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
#                                 convex = -0.04)
bnd <- INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.03)
bnd2 = INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(data.frame(dat$long, dat$lat)),
  boundary = list(bnd,bnd2),
  cutoff = 3, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
)
mesh <- make_mesh(as.data.frame(dat), c("long", "lat"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)

## Create mesh with sdmtmb
#mesh2 <- make_mesh(
#  dat,
#  xy_cols = c("long", "lat"),
#  cutoff = 1.852 # min. distance between knots in X-Y units
#)
#
#plot(mesh2$mesh, main = NA, edge.color = "grey60", asp = 1)
#points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)

## Then fit the model 
fit <- sdmTMB(
  density ~ 1,
  data = dat,
  mesh = mesh,
  #mesh = mesh2,
  family = gaussian(link = "log"),
  spatial = "on")

## Then fit the model with surf but it doesn't work
#fit <- sdmTMB(
#  abun ~ 1 + log(surf),
#  data = dat,
#  mesh = mesh,
#  #offset = log(dat$surf),
#  #mesh = mesh2,
#  family = gaussian(link = "log"),
#  spatial = "on")

## check the result ##
summary(fit)
sanity(fit)

# check quantile residuals
dat$resids <- residuals(fit) # randomized quantile residuals
qqnorm(dat$resids)
qqline(dat$resids)

# check residual pattern
ggplot(dat, aes(long, lat, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  coord_fixed()

# plot the random spatial effects
p <- predict(fit, newdata = grid)
ggplot(p, aes(long, lat, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE)
# predict spatial distribution both in link (log) space:
ggplot(p, aes(long, lat, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
# Or natural space:
ggplot(p, aes(long, lat, fill = exp(est))) +
  geom_raster() +
  labs(fill = "Intensity\n(average point density)") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

# Get the index of the area
p_sdm <- predict(fit, newdata = grid, return_tmb_object = TRUE)
index <- get_index(p_sdm, area = surfsquare, bias_correct = TRUE)
ggplot(index, aes(type, est)) +
  geom_point(colour = "grey30")+
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2)+
  labs(x = "2025", y = "Abundance (individual)")

ggplot(index, aes(type, est*mean_weight/1000)) +
  geom_point(colour = "grey30")+
  geom_errorbar(aes(ymin = lwr*mean_weight/1000, ymax = upr*mean_weight/1000),
                width = 0.2)+
  ylim(0,index$lwr*mean_weight/1000+index$upr*mean_weight/1000)+
  labs(x = "2025", y = "Biomass (t)")+
  scale_y_continuous(labels = scales::label_comma())

bio_stock_mod <- index$est*mean_weight/1000

# Get the biomass of the stock by hand
bio_stock_mod <- 0.4*sum(exp(p$est))/1000


#########################
#work on the 2025 data
#########################

## Load data
data_abun <- readRDS(
  paste(here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
        sep=""))

shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))

mean_weight <- 0.4 # weight of an adult sea cucumber in kg
surfsquare <- 1.852*1.852 # in kilometer square

## Prepare the dataframe
dat <-  as.data.frame(data_abun)
dat$density <- dat$abun*mean_weight/(dat$area/1e+06) #density in kg/km2
dat$intensity <- dat$abun/(dat$area/1e+06) #density in nb/km2

# working in UTM km instead of UTM m by divided by 1000 to help creation mesh
dat$long <- dat$X/1000
dat$lat <- dat$Y/1000

## Prepare grid for prediction

# use the sampling grid to generate a extrapolation grid
grid <- shp_grid %>% 
  dplyr::filter(Id > 85) %>%
  dplyr::filter(Id != 188) %>%
  dplyr::filter(Id != 189) %>%
  st_centroid() %>%
  mutate(long = st_coordinates(.)[,1]/1000) %>%
  mutate(lat = st_coordinates(.)[,2]/1000) %>%
  as.data.frame() %>%
  select(long,lat) 

plot(grid)


## Create mesh with inla

bnd <- INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.03)
bnd2 = INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(data.frame(dat$long, dat$lat)),
  boundary = list(bnd,bnd2),
  cutoff = 3, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic
mesh <- make_mesh(as.data.frame(dat), c("long", "lat"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)


## Then fit the model
#test distribution with fidistrplus
fw <- fitdist(data_abun$intensity, "weibull")
fg <- fitdist(data_abun$intensity, "gamma")
fln <- fitdist(data_abun$intensity, "lnorm")
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fw, fln, fg), legendtext = plot.legend)
qqcomp(list(fw, fln, fg), legendtext = plot.legend)
cdfcomp(list(fw, fln, fg), legendtext = plot.legend)
ppcomp(list(fw, fln, fg), legendtext = plot.legend)

fit1 <- sdmTMB(
  intensity ~ 1,
  data = dat,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on")

fit2 <- sdmTMB(
  intensity ~ 1,
  data = dat,
  mesh = mesh,
  family = Gamma(link = "inverse"),
  spatial = "on")

#abundance density
fit <- sdmTMB(
  intensity ~ 1,
  data = dat,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on")

#biomasse density
fit <- sdmTMB(
  density ~ 1,
  data = dat,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on")

## check the result ##
summary(fit1)
sanity(fit1)# the model have several issues so try with another  

summary(fit)
sanity(fit)

AIC(fit2)
AIC(fit)

# check quantile residuals
dat$resids <- residuals(fit) # randomized quantile residuals
qqnorm(dat$resids)
qqline(dat$resids)

# check residual pattern
ggplot(dat, aes(long, lat, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  coord_fixed()

# plot the random spatial effects
p <- predict(fit, newdata = grid)
ggplot(p, aes(long, lat, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE)
# predict spatial distribution both in link (log) space:
ggplot(p, aes(long, lat, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
# Or natural space:
ggplot(p, aes(long, lat, fill = exp(est))) +
  geom_raster() +
  labs(fill = "Intensity\n(average point density)") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

# Get the index of the area
p_sdm <- predict(fit, newdata = grid, return_tmb_object = TRUE)
index <- get_index(p_sdm, area = surfsquare, bias_correct = TRUE)
ggplot(index, aes(type, est)) +
  geom_point(colour = "grey30")+
  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2)+
  labs(x = "2025", y = "Abundance (individual)")

ggplot(index, aes(type, est*mean_weight/1000)) +
  geom_point(colour = "grey30")+
  geom_errorbar(aes(ymin = lwr*mean_weight/1000, ymax = upr*mean_weight/1000),
                width = 0.2)+
  ylim(0,index$lwr*mean_weight/1000+index$upr*mean_weight/1000)+
  labs(x = "2025", y = "Biomass (t)")+
  scale_y_continuous(labels = scales::label_comma())

bio_stock_mod <- index$est*mean_weight/1000

# TAC 1.5% of the 5e percentile biomass distribution
# the 5e percentile is replaced by lwr
TAC_mod <- (1.5/100)*(index$lwr*mean_weight/1000)



## try with grid with a finer scale
# generate a grid
res <- 0.05 # resolution in km

grid <- expand.grid(
  long = seq(min(dat$long), max(dat$long), by = res),
  lat = seq(min(dat$lat), max(dat$lat), by = res)
)

# plot the random spatial effects
p_2 <- predict(fit, newdata = grid)
ggplot(p_2, aes(long, lat, fill = omega_s)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed(expand = FALSE)
# predict spatial distribution both in link (log) space:
ggplot(p_2, aes(long, lat, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed(expand = FALSE)
# Or natural space:
ggplot(p_2, aes(long, lat, fill = exp(est))) +
  geom_raster() +
  labs(fill = "Intensity\n(average point density)") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

ggplot(p_2, aes(long, lat, fill = exp(est)*0.001^2)) +
  geom_raster() +
  labs(fill = "Abundance per m²") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

ggplot(p_2, aes(long, lat, fill = exp(est)*(mean_weight/1000))) +
  geom_raster() +
  labs(fill = "Biomass in t per km²") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

ggplot(p_2, aes(long, lat, fill = exp(est)*(mean_weight)*0.001^2)) +
  geom_raster() +
  labs(fill = "Biomass in kg per m²") +
  scale_fill_viridis_c(trans = "sqrt") +
  coord_fixed(expand = FALSE)

#map
zee <- readRDS(
  "C:/Users/ahebertb/Documents/sea-cucumbeR/SIG/SIG_Data/sf_sampling_grid.rds")


#initialization
limit_zone <- data.frame(
  long = c(min(zee[[6]][[1]][[1]][,1]),
           max(zee[[6]][[1]][[1]][,1])),
  lat = c(min(zee[[6]][[1]][[1]][,2]),
          max(zee[[6]][[1]][[1]][,2]))
)
for (indice in (1:nrow(zee))){
  min <- min(zee[[6]][[indice]][[1]][,1])
  max <- max(zee[[6]][[indice]][[1]][,1])
  if(min < limit_zone$long[1]){
    limit_zone$long[1] <- min
  }
  if(max > limit_zone$long[2]){
    limit_zone$long[2] <- max
  }
  
  min <- min(zee[[6]][[indice]][[1]][,2])
  max <- max(zee[[6]][[indice]][[1]][,2])
  if(min < limit_zone$lat[1]){
    limit_zone$lat[1] <- min
  }
  if(max > limit_zone$lat[2]){
    limit_zone$lat[2] <- max
  }
}

ggplot(p_2, aes(long*1000, lat*1000, fill = exp(est)*mean_weight/1000))+
  geom_raster()+
  scale_fill_viridis_c(trans = "sqrt") +
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0ddff"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  coord_sf(xlim = limit_zone$long, ylim = limit_zone$lat, expand = FALSE)+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet",
       fill = "Biomass in t per km²")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tl", 
                         height = unit(0.7, "cm"), width = unit(0.7, "cm")) + 
  theme()
