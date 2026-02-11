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

#########################
#work on the 2025 data
#########################

## Load data ####
data_annot <- read.csv(
  paste(here(),"/comparing_biomass_method/Data/01_annotation_summary.csv",
        sep=""), sep = ",", header = TRUE)

data_coordinates <- read.csv(
  paste(here(),"/comparing_biomass_method/Data/02_coordinates_summary.csv",
        sep=""), sep = ",", header = TRUE)

# import the area for calculate the total abundance
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))
surfarea <- as.numeric(st_area(calcul_area)/1e+6)

field_view <- 1.5 # length in meter of the field view of the camera
mean_weight <- 0.4 # mean weight observed of an adult sea cucumber in kg
surfsquare <- 1.852*1.852 # in kilometer square

## process data ####
distance <- function(x1,y1,x2,y2){
  sqrt((x2-x1)**2+(y2-y1)**2)
}

data_coordinates$haul_distance <- NA
for (row in 1:nrow(data_coordinates)){
  data_coordinates$haul_distance[row] <- distance(
    data_coordinates$X.haul_start[row],
    data_coordinates$Y.haul_start[row],
    data_coordinates$X.shoot_end[row],
    data_coordinates$Y.shoot_end[row]
  )
}

# transform the time of a haul in a number of seconds
temps <- as.POSIXct(data_annot$Time.haul_start,
                    format = ("%H:%M")) -
  as.POSIXct(data_annot$Time.shoot_end, format = ("%H:%M"))
temps <- as.numeric(temps)
temps <- temps*-60
data_annot <- data_annot %>%
  mutate(haul_duration = temps)

data_annot$annotation_rate <- data_annot$Seconds_Read/data_annot$haul_duration


## Prepare the dataframe for the model ####
dat <-  merge(data_annot,data_coordinates,by="STN")
#calculate surface sampled and density
dat$annotation_distance <- dat$annotation_rate*dat$haul_distance
dat$surface <- dat$annotation_distance*field_view

dat$density <- dat$Total_Number*mean_weight/(
  dat$surface)#/1e+06) #density in kg/km2

# format the seabed_type column
unique(dat$Seabed_Type)
dat$Seabed_Type[grep("dollar",dat$Seabed_Type)] <- c("Sand/gravel/pebble banks")
dat$Seabed_Type[grep("Sand ",dat$Seabed_Type)] <- c("Coarse pebbles")
dat$Seabed_Type[grep("sand",dat$Seabed_Type)] <- c("Coarse pebbles")
dat$Seabed_Type[grep("Sand",dat$Seabed_Type)] <- c("Sand/gravel/pebble banks")
dat$Seabed_Type[grep("Gravel",
                     dat$Seabed_Type)] <- c("Large pebbles/gravel")
dat$Seabed_Type[grep("stone",
                     dat$Seabed_Type)] <- c("Block fields/large pebbles")
dat$Seabed_Type[grep("large",
                     dat$Seabed_Type)] <- c("Block fields/large pebbles")
dat$Seabed_Type[grep("mixture",
                     dat$Seabed_Type)] <- c("Large pebbles/gravel")
dat$Seabed_Type[grep("Mixture",
                     dat$Seabed_Type)] <- c("Large pebbles/gravel")
dat$Seabed_Type[grep("Pink pebbles",
                     dat$Seabed_Type)] <- c("Stone fields")

dat$Seabed_Type <- as.factor(
  dat$Seabed_Type
)

# choose the interest coloumn only
dat <- dat %>% select(STN,
                      Total_Number,
                      surface,
                      density,
                      Seabed_Type,
                      X.centroid_track,
                      Y.centroid_track)

# working in UTM km instead of UTM m by divided by 1000 to help creation mesh
dat$long <- dat$X.centroid_track/1000
dat$lat <- dat$Y.centroid_track/1000

# working only in the "Tuyau"
dat <- dat %>% filter(!is.na(as.numeric(STN)))

## Prepare grid for prediction
# use the study area polygon to generate a extrapolation grid
grid <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_bathy.rds",
        sep = "")
)
plot(grid)

## Create mesh with inla ####
bnd <- INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.03)
bnd2 = INLA::inla.nonconvex.hull(cbind(dat$long, dat$lat),
                                 convex = -0.15)
mesh_inla <- INLA::inla.mesh.2d(
  #loc = as.matrix(as.data.frame(grid)[,c(1,2)]),
  loc = as.matrix(dat[,c(8,9)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic
mesh <- make_mesh(dat, c("long", "lat"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)

#using the polygon of the area
bnd <- INLA::inla.nonconvex.hull(cbind(grid$long, grid$lat),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(grid$long, grid$lat),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(as.data.frame(grid)[,c(1,2)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic
mesh <- make_mesh(as.data.frame(dat), c("long", "lat"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)

## Then fit the model

fit1 <- sdmTMB(
  density ~ 1,
  data = dat,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on")

fit2 <- sdmTMB(
  density ~ 0+Seabed_Type,
  data = dat,
  mesh = mesh,
  family = lognormal(link = "log"),
  spatial = "on")
#check the result
summary(fit1)
sanity(fit1)

summary(fit2)
sanity(fit2)

AIC(fit1)
AIC(fit2)

# check quantile residuals
dat$resids <- residuals(fit1) # randomized quantile residuals
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
index <- get_index(p_sdm, area = res^2, bias_correct = TRUE)
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
  long = seq(min(coordinate$X), max(coordinate$X), by = res),
  lat = seq(min(coordinate$Y), max(coordinate$Y), by = res)
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
