# this model project the point pattern from HOLOTVSPM2025 survey in the 
#projection system EPSG:4467  RGSPM06  UTM zone 21N of Saint-Pierre and Miquelon

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(fields)
library(gridExtra)

#############
#load data
#############

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

HOLOTVSPM2025 <- readRDS(paste(here(),
                               "/comparing_biomass_method/Data/03_HOLOTVSPM2025_summary.rds",
                               sep=""))

data_position_substrat <- readRDS(paste(here(),
                                        "/via3_data_exploration/Data/substrat/data_position_substrat.rds",
                                        sep=""))

data_substrate <- readRDS(paste(here(),
                                "/via3_data_exploration/Data/substrat/data_substrat_2025.rds",
                                sep=""))

########################
# Start with one PPP
########################

stn <- "96" # # stn="149" stn="179" stn="143"
PPP <- list_PPP[[stn]] # PPP<-list_PPP[["179"]] PPP<-list_PPP[["149"]]
# PPP<-list_PPP[["143"]]
PPP[["window"]][["units"]] <- list("metre","metres")
#PPP[["window"]][["yrange"]] <- PPP[["window"]][["yrange"]] - PPP[["window"]][["yrange"]][[1]]


##############
# math recall
##############
# Reflection matrix:
#    Reflection(θ)=[cos2θ   sin2θ
#                   sin2θ  −cos2θ]

# Rotation matrix:
#    Rotation(θ)=[cosθ  -sinθ
#                 sinθ   cosθ]
##############


# spatialization of the PPP ####
start_coord <- HOLOTVSPM2025[grep(stn,HOLOTVSPM2025$STN),
                             c("X.haul_start","Y.haul_start")]
end_coord <- HOLOTVSPM2025[grep(stn,HOLOTVSPM2025$STN),
                           c("X.shoot_end","Y.shoot_end")]

affine <-(end_coord$Y.shoot_end - start_coord$Y.haul_start)/
  (end_coord$X.shoot_end - start_coord$X.haul_start)

ordonnee <- mean(c(end_coord$Y.shoot_end - affine*end_coord$X.shoot_end,
                   start_coord$Y.haul_start - affine*start_coord$X.haul_start))

Cartesian_equation <- function(x,affine,ordonnee){
  y <- affine*x+ordonnee
  return(y)
}

# rotation of angle theta
theta <- atan(start_coord$X.haul_start/start_coord$Y.haul_start)
new_coord <- as.matrix(data.frame(X=PPP[["x"]],Y=PPP[["y"]]))
new_coord <- spdep::Rotation(new_coord,-theta) #sens horaire donc signe négatif

# translation at the start point of the haul
#matrix.start <- t(matrix(start_coord,2,PPP$n))
#or
matrix.start <- as.matrix(data.frame(
  X=rep(start_coord$X.haul_start,dim(new_coord)[1]),
  Y=rep(start_coord$Y.haul_start,dim(new_coord)[1])
))
new_coord <- new_coord + matrix.start

# create the polygonal window for point pattern
#The vertices must be listed anticlockwise. No vertex should be repeated 
#(i.e. do not repeat the first vertex).
poly_coord <- data.frame(
  x=c(0,1.5,1.5,0),
  y=c(PPP$window$yrange[1],PPP$window$yrange[1],
      PPP$window$yrange[2],PPP$window$yrange[2]))

poly_coord <- spdep::Rotation(poly_coord,-theta)
matrix.start <- as.matrix(data.frame(
  x=rep(start_coord$X.haul_start,nrow(poly_coord)),
  y=rep(start_coord$Y.haul_start,nrow(poly_coord))
))
poly_coord <- poly_coord + matrix.start

win <- owin(poly = poly_coord)
X <- ppp(as.data.frame(new_coord)$X,
         as.data.frame(new_coord)$Y, window = win)


############################
# Create a list with all PPP
############################

list_PPP_epsg4461 <- list()
data_position_epsg4461 <- data.frame()
stations <- names(list_PPP)

for (stn in stations){
  # call data ####
  PPP <- list_PPP[[stn]] 
  PPP[["window"]][["units"]] <- list("metre","metres")
  
  start_coord <- HOLOTVSPM2025[grep(stn,HOLOTVSPM2025$STN),
                               c("X.haul_start","Y.haul_start")]
  end_coord <- HOLOTVSPM2025[grep(stn,HOLOTVSPM2025$STN),
                             c("X.shoot_end","Y.shoot_end")]
  
  poly_coord <- data.frame(
    x=c(0,1.5,1.5,0),
    y=c(PPP$window$yrange[1],PPP$window$yrange[1],
        PPP$window$yrange[2],PPP$window$yrange[2]))
  
  # coord transformation ####
  new_coord <- as.matrix(data.frame(X=PPP[["x"]],
                                    Y=PPP[["y"]]))
  # rotation of angle theta.rotation and sign(slope)
  slope <- (end_coord$Y.shoot_end - start_coord$Y.haul_start)/
    (end_coord$X.shoot_end - start_coord$X.haul_start)
  theta.x <- atan(slope) 
  if (slope<0){
    theta.rotation <- (theta.x + 90*pi/180)
  } else {
    theta.rotation <- -(90*pi/180 - theta.x)#clockwise = negative sign
  }
  new_coord <- spdep::Rotation(new_coord, theta.rotation)
  poly_coord <- spdep::Rotation(poly_coord, theta.rotation)
  
  # translation at the start point of the haul
  matrix.start <- as.matrix(data.frame(
    X=rep(start_coord$X.haul_start,dim(new_coord)[1]),
    Y=rep(start_coord$Y.haul_start,dim(new_coord)[1])
  ))
  new_coord <- new_coord + matrix.start
  matrix.start <- as.matrix(data.frame(
    x=rep(start_coord$X.haul_start,nrow(poly_coord)),
    y=rep(start_coord$Y.haul_start,nrow(poly_coord))
  ))
  poly_coord <- poly_coord + matrix.start
  
  # create the point pattern ####
  #The vertices must be listed anticlockwise. No vertex should be repeated 
  #(i.e. do not repeat the first vertex).
  win <- owin(poly = poly_coord)
  X <- ppp(as.data.frame(new_coord)$X,
           as.data.frame(new_coord)$Y, window = win)
  X[["window"]][["units"]] <- list("metre","metres")
  
  # save the data ####
  list_PPP_epsg4461[[stn]] <- X
  data_position_epsg4461 <- rbind(data_position_epsg4461,
                                  data.frame(
                                    station = stn,
                                    X = X[["x"]],
                                    Y = X[["y"]]
                                  ))
  
}

saveRDS(list_PPP_epsg4461,
        paste(here(),"/point_process/Data/list_PPP_2025_epsg4461.rds",sep=""))
saveRDS(data_position_epsg4461,
        paste(here(),"/point_process/Data/data_position_2025_epsg4461",sep=""))
