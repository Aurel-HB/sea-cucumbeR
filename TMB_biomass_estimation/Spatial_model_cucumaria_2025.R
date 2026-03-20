# this code test a spatial model on the 2025 data using TMB

# install
install.packages("here")
install.packages("fmesher")
remotes::install_version("INLA", version = "25.06.13",
                         repos = c(getOption("repos"),
                                   INLA = "https://inla.r-inla-download.org/R/testing"),
                         dep = TRUE) 


# Load TMB
library( TMB )
library("here")
library("ggplot2")
library( dplyr )
library("sf")
library("fmesher")

# Data
setwd( paste(here(),"/TMB_biomass_estimation/Data",sep="") )
HOLOTVSPM2025 <- readRDS( "HOLOTVSPM2025_summary.rds" )
area <- readRDS("study_calcul_area.rds")
prediction_grid <- readRDS("grid_bathy.rds")

# function to optimize from the package TMBhelper
setwd( paste(here(),"/TMB_biomass_estimation",sep="") )
source("fit_tmb.R")
source("TMBAIC.R")

# working only in the "Tuyau"
dat <- HOLOTVSPM2025 %>% filter(!is.na(as.numeric(STN)))

# choose the interest coloumn only
dat <- dat %>% select(biomass_density,
                      X.centroid_track,
                      Y.centroid_track,
                      surface)
names(dat) <- c("density", "long", "lat","area")

# working in UTM km instead of UTM m by divided by 1000 to help creation mesh
dat$long <- dat$long/1000
dat$lat <- dat$lat/1000

# transform the density to be in kg/km2
dat$density <- dat$density*1000
ggplot()+
  geom_point(aes(dat$long,dat$lat,color=log(dat$density)))+
  scale_color_gradient2()+
  theme(aspect.ratio = 3)

## Create mesh with fmesher ####
bnd <- fm_nonconvex_hull(cbind(prediction_grid$long, prediction_grid$lat),
                         convex = -0.04)
bnd2 = fm_nonconvex_hull(cbind(prediction_grid$long, prediction_grid$lat),
                         convex = -0.15)
mesh <- fm_mesh_2d(
  #loc = as.matrix(as.data.frame(grid)[,c(1,2)]),
  loc = as.matrix(dat[,c(2,3)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic
# show mesh
plot(mesh, main = NA, edge.color = "grey60", asp = 1)
points(dat$long, dat$lat, pch = 19, col = "red",cex = 0.3)

# show data before model ####
#show knot
mesh_loc <- as.data.frame(mesh$loc)
mesh_loc[,3] <- seq(1,mesh$n,1)
names(mesh_loc) <- c("X","Y","idx")
ggplot(mesh_loc,aes(X,Y, label = idx))+geom_text()+theme(aspect.ratio = 3)

ggplot()+
  geom_point(aes(dat$long,dat$lat,color=log(dat$density)))+
  scale_colour_gradient2()+
  geom_text(aes(mesh_loc$X,mesh_loc$Y, label = mesh_loc$idx))+
  geom_text(aes(dat$long,dat$lat,label=mesh$idx$loc),colour="blue", size = 3)+
  theme(aspect.ratio = 3)

###### Method spatial model -- Optimize using TMB ####
# Create matrices in fmesher/INLA for spde
# spde <- fm_fem(mesh, order=2)
# M0 = κ⁴ * c0 + 2 * κ² * c1 + c2
# M1 = κ⁴ * g1 + 2 * κ² * g2
# M2 = τ² * c0
spde <- INLA::inla.spde2.matern(mesh, alpha=2)

# create projection matrix from vertices to samples 
A_is = fm_evaluator( mesh , loc=as.matrix(dat[,c(2,3)]) )$proj$A

# Step 1 -- make and compile template file
compile( "spatial_gamma.cpp" )

# Step 2 -- build inputs and Object
dyn.load( dynlib("spatial_gamma") )
# SPDE-based
# Build object
Data = list("c_i"=as.vector(dat$density),
            #"j_i"=mesh$idx$loc-1,
            "area_i"=dat$area,
            "M0"=spde$param.inla$M0,
            "M1"=spde$param.inla$M1,
            "M2"=spde$param.inla$M2,
            "A_is"=A_is)

Params = list( "beta0"=0,
               "ln_tau"=0,
               "ln_kappa"=0,
               "ln_alpha"=0,
               "omega_s"=rnorm(nrow(spde$param.inla$M0)) )

Obj = MakeADFun( data=Data, parameters=Params, random="omega_s", DLL="spatial_gamma" )

# Step 3 -- test and Optimize
# Optimize
Opt_spde = fit_tmb( obj=Obj, newtonsteps=1, bias.correct=TRUE )
h_spde = Obj$env$spHess(random=TRUE)
report_spde = Obj$report()

# Results
SD = sdreport( Obj ) # standard errors
SD$pdHess # check the validity of the Hessian matrix
plot(report_spde[["omega_i"]] ~ dat$density, xlab="Densité observée (kg/km²)", 
     ylab="Effet spatial projeté (omega_i)")

# likelihood profile diagnostic TMB model
prof_beta0  <- tmbprofile(Obj,name="beta0")
confint(prof_beta0)
plot(prof_beta0)

prof_ln_kappa  <- tmbprofile(Obj,name="ln_kappa")
confint(prof_ln_kappa)
plot(prof_ln_kappa)

prof_ln_tau  <- tmbprofile(Obj,name="ln_tau")
confint(prof_ln_tau)
plot(prof_ln_tau)

prof_ln_alpha  <- tmbprofile(Obj,name="ln_alpha")
confint(prof_ln_alpha)
plot(prof_ln_alpha)

# Prediction ####
# create projection matrix from vertices to grid
coordinates <- cbind(prediction_grid$long, prediction_grid$lat)
#A_gs = fm_evaluator(mesh,
#                    loc=st_coordinates(st_centroid(prediction_grid)))$proj$A
A_gs <- fm_evaluator(mesh,
                     loc=coordinates)$proj$A


# extract random spatial effect
omega_s <- report_spde$omega_s
# projection on grid
omega_grid <- A_gs %*% omega_s
# Association to grid
prediction_grid$prediction <- as.numeric(omega_grid)
# Visualization
plot(prediction_grid["prediction"], main = "Prédiction du champ gaussien")
ggplot(prediction_grid, aes(long, lat, fill = prediction)) +
  geom_raster() +
  scale_fill_gradient2() +
  coord_fixed()+
  theme(aspect.ratio = 3)

# create the prediction map 
prediction_grid$density <- exp(
  report_spde$beta0+prediction_grid$prediction
)

ggplot(prediction_grid, aes(long, lat, fill = density)) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed()+
  theme(aspect.ratio = 3)

Biomass_tot <- sum(prediction_grid$density*(0.5**2))/1000
