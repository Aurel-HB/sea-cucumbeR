# this script compile all the work tested to compare the Gamma and Lognormal
# distribution in the context of sea cucumber biomass densities distribution
# model. This code produce the analysis and result for the paper associated to
# this study.

# load packages ####
library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sdmTMB)
library(INLA)
library(RColorBrewer)
library(ape)
#library(cowlot)

# load data ####
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

data_holotv <- readRDS(
  paste(here(),"/process_spatiotemp_data/Data/processed/data_abun_tot_cov.rds",
        sep = "")
)

grid_bathy <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_bathy.rds",
        sep = "")
)

grid_tmp <- readRDS(
  paste(here(),
      "/environment_exploration/Environment_Data/processed/grid_bottom_tmp.rds",
      sep = "")
)

grid_chl <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_chl.rds",
        sep = "")
)

# prepare data ####
## Prepare the grid for prediction ####
# The survey is only in may so we keep the may data
grid_tmp <- grid_tmp[,grep(pattern = "May",names(grid_tmp))]
grid_chl <- grid_chl[,grep(pattern = "May",names(grid_chl))]

# The grid have to be replicated for each year with the environmental covariate values
grid_tot <- data.frame()
for (i in 1:5){
  data <- st_join(grid_tmp[,i],grid_bathy)
  data <- st_join(data,grid_chl[,i])
  data <- data %>% mutate(year = 2020+i)
  names(data) <- c("bottomT","long","lat","bathy","chla","geometry","year")
  grid_tot <- rbind(grid_tot,data[,c(1,2,3,4,5,7,6)])
} # the lat and long are in km for modelling
rm(data)
## Prepare the date for modelling ####
# synchronize variable name
names(grid_tot)[2:3] <- c("X","Y")
# keep only station in the study area
data_holotv <- 
  st_join(calcul_area,data_holotv) # to get intersection of points and poly

# convert abundance in biomass
data_holotv$biomass <- data_holotv$abun*0.4
data_holotv$density.t_km2 <- (data_holotv$biomass/1000)/(data_holotv$area/1e6)

# convert distance from meter to kilometer
data_holotv$X <- data_holotv$X/1000
data_holotv$Y <- data_holotv$Y/1000
data_holotv$area <- data_holotv$area/1e6

# keep only the column needed
data_estimate <- as.data.frame(data_holotv) %>% 
  select(station,X,Y,density.t_km2,year)
grid_proj <- as.data.frame(grid_tot)[,1:6]

# extract the surface of the study area
stock_surface <- as.numeric(st_area(calcul_area)/1e6)


###-###-###-###-###-###-###-###
# Mesh sensibility analysis ####
###-###-###-###-###-###-###-###

## creation of the mesh with R-INLA ####
#using the grid of the area
bnd <- INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(as.data.frame(grid_bathy)[,c(1,2)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic

mesh_regular <- make_mesh(data_estimate,
                       c("X", "Y"), mesh = mesh_inla)
plot(mesh_regular$mesh, main = NA, edge.color = "grey60", asp = 1)
points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(data_estimate[,c("X","Y")]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic

mesh_data_fit <- make_mesh(data_estimate,
                          c("X", "Y"), mesh = mesh_inla)
plot(mesh_data_fit$mesh, main = NA, edge.color = "grey60", asp = 1)
points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)

## fit model for the 2 meshes####
data_estimate$time <- data_estimate$year

fit_gamma_regular <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year),#<< fixed intercept ignoring time
  data = data_estimate,
  mesh = mesh_regular,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)

fit_gamma_data_fit <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year),#<< fixed intercept ignoring time
  data = data_estimate,
  mesh = mesh_data_fit,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)

fit_logn_regular <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year),#<< fixed intercept ignoring time
  data = data_estimate,
  mesh = mesh_regular,
  family = lognormal(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)

fit_logn_data_fit <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year),#<< fixed intercept ignoring time
  data = data_estimate,
  mesh = mesh_data_fit,
  family = lognormal(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)

## check the result ####
sanity(fit_gamma_regular)
sanity(fit_gamma_data_fit)
sanity(fit_logn_regular)
sanity(fit_logn_data_fit)

par(mfrow = c(2,2))
data_estimate$resids_gamma_regular <- residuals(fit_gamma_regular) 
# randomized quantile residuals
qqnorm(data_estimate$resids_gamma_regular, main = "Normal Q-Q Plot Gamma Regular Mesh")
qqline(data_estimate$resids_gamma_regular)

data_estimate$resids_gamma_data_fit <- residuals(fit_gamma_data_fit) 
# randomized quantile residuals
qqnorm(data_estimate$resids_gamma_data_fit, main = "Normal Q-Q Plot Gamma Data-fit Mesh")
qqline(data_estimate$resids_gamma_data_fit)

data_estimate$resids_logn_regular <- residuals(fit_logn_regular) 
# randomized quantile residuals
qqnorm(data_estimate$resids_logn_regular, main = "Normal Q-Q Plot Lognormal Regular Mesh")
qqline(data_estimate$resids_logn_regular)

data_estimate$resids_logn_data_fit <- residuals(fit_logn_data_fit) 
# randomized quantile residuals
qqnorm(data_estimate$resids_logn_data_fit, main = "Normal Q-Q Plot Lognormal Data-fit Mesh")
qqline(data_estimate$resids_logn_data_fit)
par(mfrow = c(1,1))

### Moran index per year ###
moran <- data.frame()
sim <- c()
for (i in 1:1000){
  data_estimate$resids_gamma_regular <- residuals(fit_gamma_regular) 
  data_estimate$resids_gamma_data_fit <- residuals(fit_gamma_data_fit)
  data_estimate$resids_logn_regular <- residuals(fit_logn_regular)
  data_estimate$resids_logn_data_fit <- residuals(fit_logn_data_fit)
  
  for (time in unique(data_estimate$year)){
    # calcul weight
    dists <- as.matrix(dist(data_estimate[
      grep(time,data_estimate$year), c("X", "Y")]))
    inv_dists <- 1 / dists
    diag(inv_dists) <- 0
    inv_dists[is.infinite(inv_dists)] <- 0
    
    #gamma regular
    data <- Moran.I(data_estimate$resids_gamma_regular[data_estimate$year==time],
                    inv_dists, scaled = TRUE)
    data$year <- time
    data$mesh <- "Regular"
    data$model <- "Gamma"
    moran <- rbind(moran, data)
    
    #gamma data-fit
    data <- Moran.I(data_estimate$resids_gamma_data_fit[data_estimate$year==time],
                    inv_dists, scaled = TRUE)
    data$year <- time
    data$mesh <- "Data-fit"
    data$model <- "Gamma"
    moran <- rbind(moran, data)
    
    #logn regular
    data <- Moran.I(data_estimate$resids_logn_regular[data_estimate$year==time],
                    inv_dists, scaled = TRUE)
    data$year <- time
    data$mesh <- "Regular"
    data$model <- "Lognormal"
    moran <- rbind(moran, data)
    
    #logn data-fit
    data <- Moran.I(data_estimate$resids_logn_data_fit[data_estimate$year==time],
                    inv_dists, scaled = TRUE)
    data$year <- time
    data$mesh <- "Data-fit"
    data$model <- "Lognormal"
    moran <- rbind(moran, data)
  }
  sim <- c(sim, rep(paste("sim",i,sep="_"),16))
}
moran$sim <- sim
moran$signif <- FALSE
moran$signif[moran$p.value<0.05] <- TRUE
summary_moran <- data.frame(moran[grep("sim_1000",moran$sim),
                                  c("year","mesh","model")],
                            percent_ns=0)
for (i in 1:nrow(summary_moran)){
  time = summary_moran$year[i]
  type = summary_moran$mesh[i]
  distrib = summary_moran$model[i]
  data <- moran %>% filter(year==time) %>% filter(mesh==type) %>%
    filter(model==distrib)
  summary_moran$percent_ns[i] <- length(grep(TRUE,data$signif))*100/nrow(data)
}
rm(time,type,data)

ggplot(data = summary_moran, aes(x = as.factor(year), y = as.factor(mesh)))+
  #geom_raster(aes(fill = percent_ns))+
  geom_point(aes(colour = percent_ns),shape=15, size=13)+
  #scale_colour_viridis_c(option = "plasma")+
  scale_color_gradient2(low="slateblue",mid = "#eeeeaa",high = "red3",
                        midpoint = 10)+
  geom_point(shape=22, size=15)+
  geom_text(aes(label = percent_ns))+
  facet_wrap(~model, nrow=1)+
  labs(fill = "% significant Moran's index", x = "Year", y = "Mesh type")
