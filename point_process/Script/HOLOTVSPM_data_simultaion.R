# This code try to simulate a dataset that could be from the survey HOLOTVSPM
# First the number of Sea cucumber/haul is simulate  probability distribution
# Inside each haul, a ppp simulate the distribution of the species
library(here)
library(sf)
library(spatstat)
library(dplyr)

# creation of object of class owin
# the haul length is 500m and the gopro view 1m50
win <- owin(xrange = c(0, 1.5), yrange = c(0, 500))
plot(win)


#import the ZEE grid
shp_grid <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""),
                    layer = "grille_zee_spm")
# filter the square that are in the sampling zone
shp_grid <- shp_grid %>% dplyr::filter(!is.na(Echantillo))

# extract the centroid of each square to facilitate the plot of station
sampling_centroid <- st_centroid(shp_grid) %>%
  mutate(X = st_coordinates(.)[,1]) %>%
  mutate(Y = st_coordinates(.)[,2])