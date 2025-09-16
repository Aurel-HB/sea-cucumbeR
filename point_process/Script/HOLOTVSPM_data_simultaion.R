# This code try to simulate a dataset that could be from the survey HOLOTVSPM

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(tweedie)
library(ggspatial)


#############################################################
# Part 1 - Simulation by dredge haul
#############################################################

# First the number of Sea cucumber/haul is simulate  probability distribution
# Inside each haul, a ppp simulate the distribution of the species

# creation of object of class owin
# the haul length is 500m and the gopro view 1m50
win <- owin(xrange = c(0, 1.5), yrange = c(0, 500))
plot(win)


#import the ZEE grid
shp_grid <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""),
                    layer = "grille_zee_spm")
# filter the square that are in the sampling zone
shp_grid <- shp_grid %>% dplyr::filter(!is.na(Echantillo))

#creation color vector based on Echantillo displaying the sampling plan
shp_grid <- shp_grid %>%
  mutate(Echantillo = ifelse(Id == 193, "VIDEO", Echantillo)) %>%
  mutate(Echantillo = ifelse(Id == 104, "PAS ECHANTILLONNAGE", Echantillo))

shp_grid$color_sampling <- c(NA)

shp_grid$color_sampling[grep(pattern = "VIDEO", shp_grid$Echantillo)] <- 
  "lightsalmon3"
shp_grid$color_sampling[grep(pattern = "PECHE", shp_grid$Echantillo)] <- 
  "seagreen3"
shp_grid$color_sampling[grep(pattern = "CROCHE", shp_grid$Echantillo)] <- 
  "#990000"
shp_grid$color_sampling[grep(pattern = "ECHANTILLONNAGE",
                             shp_grid$Echantillo)] <- "#fffbfd00"
shp_grid$color_sampling[is.na(shp_grid$Echantillo)] <- "#fffbfd00"

# save the sf object for mapping the area and doing beautiful map
saveRDS(shp_grid,
        file = paste(here(),"/SIG/SIG_Data/sf_sampling_grid.rds", sep=""))
shp_grid <- readRDS(paste(here(),"/SIG/SIG_Data/sf_sampling_grid.rds", sep=""))

# extract the centroid of each square to facilitate the plot of station
sampling_centroid <- st_centroid(shp_grid) %>%
  mutate(X = st_coordinates(.)[,1]) %>%
  mutate(Y = st_coordinates(.)[,2])


# load a example of data to calibrate the simulation
holo23 <- read.csv(paste(here(),"/point_process/Data/holo23.csv", sep=""),
                   sep = ";")
source(paste(here(),"/point_process/Script/fct_probability_distribution.R",
             sep=""))

# choice of a distribution for the biomass  ####
loi <- "Lognormale"
ggplot(data = holo23, aes(x = Biomass, y = after_stat(density))) +
  geom_histogram(colour = "black", fill = "grey", bins = 100) +
  geom_density(alpha = .3, fill = "blue") +
  geom_area(aes( y= probability_distribution(Biomass, loi)),
            color="darkslategray", fill = "darkseagreen",
            alpha = 0.4, linewidth = 1)+
  geom_area(aes( y= dgamma(Biomass, shape = 1.1, scale = 170)),
            color="darkslategray", fill = "darkorange",
            alpha = 0.2, linewidth = 1)+
  xlab(holo23$Biomass) +
  theme_minimal()

# choice of a distribution for the abundance ####
loi <- "Lognormale"
ggplot(data = holo23, aes(x = Abun, y = after_stat(density))) +
  geom_histogram(colour = "black", fill = "grey", bins = 100) +
  geom_density(alpha = .3, fill = "blue") +
  #geom_area(aes( y= probability_distribution(Abun, loi)),
  #          color="darkslategray", fill = "darkseagreen",
  #          alpha = 0.4, linewidth = 1)+
  geom_area(aes( y= dlnorm(Abun, 
                           meanlog = 5,
                           sdlog = 1)),
            color="darkslategray", fill = "darkseagreen",
            alpha = 0.4, linewidth = 1)+
  geom_area(aes( y= dgamma(Abun, shape = 1.1, scale = 400)),
            color="darkslategray", fill = "darkorange",
            alpha = 0.2, linewidth = 1)+
  xlab(holo23$Abun) +
  theme_minimal()

#### Prepare abundance table ####
# keep the station sampled
shp_grid_sampling <- merge(holo23["Station"], sampling_centroid,
                           by.x = "Station", by.y = "Id")

#set.seed(26) # for reproducibility
set.seed(13)

#add the abun from a probability distribution
#shp_grid_sampling$abun <- rgamma(n = length(shp_grid_sampling$Station),
#                                 shape = 1.1, scale = 400)
shp_grid_sampling$abun <- rlnorm(n = length(shp_grid_sampling$Station),
                                 meanlog = 5, sdlog = 1)
# transform abundance in individual per meter square that will be used for PPP
# indeed the poisson process used an intensity by unit area for simulation
win.area <- win$xrange[2]*win$yrange[2]
shp_grid_sampling$abun <- shp_grid_sampling$abun/win.area

# keep interested column and spatialised shp_grid_sampling
shp_grid_sampling <- data.frame(shp_grid_sampling[,1:2],
                                shp_grid_sampling[,8:9],
                                "Intensity" = shp_grid_sampling[,10],
                                shp_grid_sampling[,6])
shp_grid_sampling <- st_as_sf(shp_grid_sampling)



# save and load shp_grid_sampling
saveRDS(shp_grid_sampling,
        file = paste(here(),"/point_process/Output/holo_simu_intensity.rds",
                     sep=""))
shp_grid_sampling <- readRDS(
  paste(here(), "/point_process/Output/holo_simu_intensity.rds", sep=""))

#### Prepare PPP ####
# create a ppp per station using the density as intensity

#phom <- rpoispp(lambda = 1.366813, win = win)
#phom$n

set.seed(20) # for reproducibility
list_PPP <- list()
for (i in 1:length(shp_grid_sampling$Station)){
  #Homogeneous point process
  list_PPP[[i]] <- rpoispp(lambda = shp_grid_sampling$Intensity[i], win = win)
  #name by station
  names(list_PPP)[i] <-   paste("STN",shp_grid_sampling$Station[i], sep = "")
}

#### Prepare the table with abundance ####
# the abundance is the number of point per PPP in liste_PPP
data_abun <- st_as_sf(data.frame("station" = names(list_PPP),
                        "abun" = 0,
                        shp_grid %>%
                          dplyr::filter(Id %in% shp_grid_sampling$Station) %>%
                          select('geometry')
                        ))
for (i in 1:length(data_abun$abun)){
  data_abun$abun[i] <- list_PPP[[i]]$n
}
data_abun$station <- gsub(pattern = "STN",replacement = "",data_abun$station)

saveRDS(data_abun,
        file = paste(here(),"/point_process/Output/holo_simu_abun.rds", sep=""))


#### Prepare the table with individual position ####
###### test for one haul ###
##From ppp to sf
# ppp object
X <- list_PPP[[1]]
# create data frame with coordinates
ddf <- data.frame(presence = 1, x = X$x, y = X$y)
# spatialised in sampling area with centroid of a sampled square
ddf$x <- ddf$x + shp_grid_sampling$X[1]
ddf$y <- ddf$y + shp_grid_sampling$Y[1]
# create sf object with data frame and name of coordinates
d <- st_as_sf(ddf, coords = c("x", "y"))

## plot the PPP on the map ###
ggplot()+
  #geom_sf(data = shp_grid, color = "black",
  #        fill = shp_grid$color_sampling)+
  geom_point(data = ddf, aes(x = x, y= y, colour = "green"))+
  theme(title = element_text(color = "black",face = "bold"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(.3, "cm"))+
  annotation_scale(location = "bl", line_width = .3) +
  annotation_north_arrow(location = "tr", height = unit(0.5, "cm"),
                         width = unit(0.5, "cm")) +
  theme()

###### creation of data_position ###
##From ppp to sf

data_position <- data.frame()
 
for (i in 1:length(list_PPP)){
  # ppp object
  X <- list_PPP[[i]]
  # create data frame with coordinates
  ddf <- data.frame(presence = 1, station = names(list_PPP)[i],
                    x = X$x, y = X$y)
  # spatialised in sampling area with centroid of a sampled square
  ddf$x <- ddf$x + shp_grid_sampling$X[i]
  ddf$y <- ddf$y + shp_grid_sampling$Y[i]
  data_position <- rbind(data_position, ddf)
}
# create sf object with data frame and name of coordinates
data_position <- data_position %>%
  mutate(long = x) %>%
  mutate(lat = y)
data_position <- st_as_sf(data_position, coords = c("x", "y")) %>% 
  st_set_crs(., 4467)

## plot the PPP on the map ###
ggplot()+
  geom_sf(data = shp_grid, color = "black",
          fill = shp_grid$color_sampling)+
  geom_point(data = data_position, aes(x = long, y= lat))+
  theme(title = element_text(color = "black",face = "bold"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.key.size = unit(.3, "cm"))+
  annotation_scale(location = "bl", line_width = .3) +
  annotation_north_arrow(location = "tr", height = unit(0.5, "cm"),
                         width = unit(0.5, "cm")) +
  theme()

data_position$station <- gsub(pattern = "STN",
                              replacement = "",data_position$station)

saveRDS(data_position,
        file = paste(here(),"/point_process/Output/holo_simu_position.rds",
                     sep=""))

#############################################################
# Part 2 - Simulation of the complete zone by a point process
#############################################################


