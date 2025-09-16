## this script process the annotation data from via to a datatable and a sf ####
# the datatable will summarise the information from csv
# the format of the table is :
#### shape <- the shape of the annotation, 1 for point and 2 for square
#### x_pixel and y_pixel <- coordinates in pixel in the frame
#### time <- time in second
#### station <- represent the id of the sampled station
# the sf object will associate the gps data and the datatable
# it will give the a proxy of the exact position of each sea cucumber seen

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(XML)

####################################################
#PART 1 - Work on the datatable from the annotation 
####################################################

#############
#load data
#############
ici <- paste(here(),"/via3_data_exploration/Data/raw/", sep="")

#load the data
raw.data <- data.frame()
myxls <- list.files(ici)[grepl('.csv',list.files(ici))]


for (i in 1:length(myxls)){
  data <- read.csv(paste(ici,myxls[i], sep = ""), 
                   sep = ",", header = TRUE, skip = 9)
  raw.data <- rbind(raw.data, data)
}

data <- raw.data

#############
#prepare data
#############

# extract type of annotation ####
shape <- c()
for (variable in 1:length(data$metadata)) {
  test <- data$spatial_coordinates[variable]
  test <- substr(test,start = 2, stop = 2)
  shape <- c(shape,as.numeric(test))
}
data <- data.frame(data, shape)

# extraction coordinates ####
extract <- data$spatial_coordinates

x_pixel <- c()
y_pixel <- c()
for (variable in 1:length(extract)) {
  test <- extract[variable]
  test <- strsplit(test,split = ",")
  x_pixel <- c(x_pixel,as.numeric(test[[1]][2]))
  y <- test[[1]][3]
  y_pixel <- c(y_pixel,as.numeric(gsub(pattern = "]",replacement = "",y)))
}

data <- data.frame(data,x_pixel,y_pixel)

time <- c()
for (variable in 1:length(data$temporal_coordinates)) {
  test <- data$temporal_coordinates[variable]
  test <- substr(test,start = 2, stop = nchar(test))
  test <- gsub(pattern = "]",replacement = "",test)
  time <- c(time,as.numeric(test))
}

data <- data.frame(data, time)

# extract station ####
station <- c()
for (variable in 1:length(data$temporal_coordinates)) {
  test <- data$file_list[variable]
  test <- substr(test,start = 13, stop = 15)
  test <- gsub(pattern = "_",replacement = "",test)
  #station <- c(station,as.numeric(test))
  station <- c(station,test)
}

data <- data.frame(data, station)

#############
data <- data[,7:11]

saveRDS(data,
        paste(here(),
              "/via3_data_exploration/Data/processed/data_pixel_2025.rds",
                    sep=""))

data_pixel <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_pixel_2025.rds",
                            sep=""))

###############################################################################
#PART 2 - Work on the sf object to connect the datatable with gps position 
###############################################################################

# import and extract data #####
#import the tracks ###
gpx_parsed <- htmlTreeParse(file = paste(here(),
                                         "/SIG/SIG_Data/tracks_2025.gpx",sep=""),
                            useInternalNodes = TRUE)

tracks <- data.frame(
  latitude = as.numeric(
    xpathSApply(doc = gpx_parsed, path = "//trkpt", fun = xmlAttrs)["lat", ]),
  longitude = as.numeric(
    xpathSApply(doc = gpx_parsed, path = "//trkpt", fun = xmlAttrs)["lon", ]),
  time = xpathSApply(doc = gpx_parsed, path = "//trkpt/time", fun = xmlValue)
)

# import the annotation summary ###
count_tot <- read.csv(paste(here(),
                            "/via3_data_exploration/Data/HOLOTVSPM2025.csv",
                            sep=""), 
                      sep = ",", header = TRUE, skip = 13)
count_tot <- count_tot %>%
  filter(!is.na(x = Date))

# import the feuille_route to have stop and start coordinates ###
feuille_route <- read.csv(paste(here(),
                            "/via3_data_exploration/Data/feuille_route_2025.csv",
                            sep=""), sep = ";", header = TRUE,
                          fileEncoding='latin1', check.names=F)

feuille_route <- feuille_route %>% filter(`Validité du trait`== "oui")

#create 2 dataframe with the good coordinates and the important information ####
feuille_route_essential_start <- feuille_route[, c("Date",
                                                   "N° Station",
                                                   "Heure début Virage",
                                                   "Latitude début Virage (N)",
                                                   "Longitude début Virage (W)")]

feuille_route_essential_stop <- feuille_route[, c("Date",
                                                  "N° Station",
                                                  "Heure fin Filage",
                                                  "Latitude fin Filage (N)",
                                                  "Longitude fin Filage (W)")]

names(feuille_route_essential_start) <- c("Date", "Station", "heure",
                                          "latitude", "longitude") 
names(feuille_route_essential_stop) <- c("Date", "Station", "heure",
                                         "latitude", "longitude")
#transform the coordinates to use in sf
source(paste(here(),
             "/via3_data_exploration/fct_degree_decimal.R",
             sep=""))

for (i in 1:length(feuille_route_essential_start$Date)){
  feuille_route_essential_start$latitude[i] <-
    (deg_dec(feuille_route_essential_start$latitude[i]))
  feuille_route_essential_start$longitude[i] <-
    (-deg_dec(feuille_route_essential_start$longitude[i]))
  feuille_route_essential_stop$latitude[i] <-
    (deg_dec(feuille_route_essential_stop$latitude[i]))
  feuille_route_essential_stop$longitude[i] <-
    (-deg_dec(feuille_route_essential_stop$longitude[i]))
}

feuille_route_essential_start$latitude <- as.numeric(
  feuille_route_essential_start$latitude)
feuille_route_essential_start$longitude <- as.numeric(
  feuille_route_essential_start$longitude)
feuille_route_essential_stop$latitude <- as.numeric(
  feuille_route_essential_stop$latitude)
feuille_route_essential_stop$longitude <- as.numeric(
  feuille_route_essential_stop$longitude)

# transform coordinates to use the SPM projection
source(paste(here(),
             "/via3_data_exploration/fct_WG84_WGSPM.R",
             sep=""))

coord <- transfo_WG84_WGSPM(feuille_route_essential_start[,c(2,5,4)])
feuille_route_essential_start$latitude <- coord$coords.x2
feuille_route_essential_start$longitude <- coord$coords.x1
coord <- transfo_WG84_WGSPM(feuille_route_essential_stop[,c(2,5,4)])
feuille_route_essential_stop$latitude <- coord$coords.x2
feuille_route_essential_stop$longitude <- coord$coords.x1
rm(coord)

# transform the time of a haul in a number of seconds
temps <- as.POSIXct(feuille_route_essential_start$heure,
                    format = ("%Hh%M")) -
  as.POSIXct(feuille_route_essential_stop$heure, format = ("%Hh%M"))
temps <- as.numeric(temps)
temps <- temps*-60
feuille_route_essential_start <- feuille_route_essential_start %>%
  mutate(temps = temps)

# create polygons that represent the hauls
#### test for one haul : use of the station 191 ####
## Make a set of coordinates that represent vertices
## with longitude and latitude 
#x_start <- round(feuille_route_essential_start$longitude[1], digits = 3)
#y_start <- round(feuille_route_essential_start$latitude[1] , digits = 3)
#x_stop  <- round(feuille_route_essential_stop$longitude[1] , digits = 3)
#y_stop  <- round(feuille_route_essential_stop$latitude[1]  , digits = 3)
#
## we want a width of 1.5 meter so we used the properties that diagonal of a 
## square is sqrt(2)*side
#side <- round(1.5/sqrt(2), digits = 1)
#
#x_coords <- c(x_start, x_stop,x_stop+side,x_start+side,x_start)
#y_coords <- c(y_start, y_stop,y_stop+side,y_start+side,y_start)
#
## convert in polygon then into a Polygon class
#poly1 <- sp::Polygon(cbind(x_coords,y_coords))
#firstPoly <- sp::Polygons(list(poly1), ID = "A")
## convert into SpatialPolygons
#firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
#firstSpatialPoly
#plot(firstSpatialPoly)
#
#start <- st_as_sf(data.frame(x=x_start, y=y_start),
#                  coords = c("x","y"),
#                  crs = CRS("+init=epsg:4467"))
#stop <- st_as_sf(data.frame(x=x_stop, y=y_stop),
#                  coords = c("x","y"),
#                  crs = CRS("+init=epsg:4467"))
#distance <- as.numeric(st_distance(start,stop))
#length_pixel <- 1920 # the videos are the format 1920*540
#height <- 1/380 #for the first part of the video 380 pixel = 1m
#
#pixel_191 <- data_pixel %>% filter(station == "191")
#
#win <- owin(xrange = c(0, 1.5), yrange = c(0,distance))
#plot(win)
#
#pixel_191 <- pixel_191 %>%
#  mutate(Y = time*distance/temps[1]) %>%
#  mutate(X = x_pixel*1.5/length_pixel) %>%
#  mutate(test = time*distance/temps[1]+y_pixel*1.5/length_pixel)
## approximation that a pixel represent the same ratio in meter for 
## length and height
#
#ggplot()+
#  geom_point(data = pixel_191, aes(x = X, y= test), colour = "green")+
#  geom_point(data = pixel_191, aes(x = X, y= Y), colour = "black")
#
# spatialization for all the data ####

data_position <- data.frame()
for (stn in unique(data_pixel$station)){
  pixel <- data_pixel %>% filter(station == stn)
  data_start <- feuille_route_essential_start %>% filter(Station == stn)
  data_stop <- feuille_route_essential_stop %>% filter(Station == stn)
  
  # Make a set of coordinates that represent the haul of the station "stn" 
  x_start <- round(data_start$longitude[1], digits = 3) # digit 3 for millimeter
  y_start <- round(data_start$latitude[1] , digits = 3)
  x_stop  <- round(data_stop$longitude[1] , digits = 3)
  y_stop  <- round(data_stop$latitude[1]  , digits = 3)
  
  # calculate the distance of the haul in the SPM projection
  start <- st_as_sf(data.frame(x=x_start, y=y_start),
                    coords = c("x","y"),
                    crs = CRS("+init=epsg:4467"))
  stop <- st_as_sf(data.frame(x=x_stop, y=y_stop),
                   coords = c("x","y"),
                   crs = CRS("+init=epsg:4467"))
  distance <- as.numeric(st_distance(start,stop))
  
  length <- 1.5/1920 # the videos are the format 1920*540 and the the GoPro 
  # enable a vision of field of 1.5m so the ratio is 1.5/1920
  height <- 1/380 #ratio for the first part of the video 380 pixel = 1m
  
  pixel <- pixel %>%
    mutate(X = x_pixel*length) %>%
    mutate(Y = time*distance/data_start$temps[1] + y_pixel*height)
  
  data_position <- rbind(data_position, pixel[,c(5,1,6,7,4)])
}

saveRDS(data_position,
        paste(here(),
              "/via3_data_exploration/Data/processed/data_position_2025.rds",
              sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

#####
# create a data with the total of sea cucumber per station and
# the surface sampled
#####
tuyau_grid <- readRDS(paste(
  here(),"/SIG/SIG_Data/sf_sampling_grid.rds",
  sep=""))

sampling_grid <- tuyau_grid[
  grep(pattern = "VIDEO", tuyau_grid$Echantillo),] %>%
  mutate(abun = 0) %>%
  mutate(tot_time = 0) %>%
  mutate(distance = 0)

# extract the centroid of each square to facilitate the plot of station
sampling_grid <- st_centroid(sampling_grid) %>%
  mutate(X = st_coordinates(.)[,1]) %>%
  mutate(Y = st_coordinates(.)[,2])

for (stn in as.character(sampling_grid$Id)){
  #calculate the abundance
  data <- data_position %>% filter(station == stn)
  abun <- nrow(data)
  
  #calculate the distance needed for the surface
  data_start <- feuille_route_essential_start %>% filter(Station == stn)
  data_stop <- feuille_route_essential_stop %>% filter(Station == stn)
  #### Make a set of coordinates that represent the haul of the station "stn" 
  x_start <- round(data_start$longitude[1], digits = 3) # digit 3 for millimeter
  y_start <- round(data_start$latitude[1] , digits = 3)
  x_stop  <- round(data_stop$longitude[1] , digits = 3)
  y_stop  <- round(data_stop$latitude[1]  , digits = 3)
  
  #### calculate the distance of the haul in the SPM projection
  start <- st_as_sf(data.frame(x=x_start, y=y_start),
                    coords = c("x","y"),
                    crs = CRS("+init=epsg:4467"))
  stop <- st_as_sf(data.frame(x=x_stop, y=y_stop),
                   coords = c("x","y"),
                   crs = CRS("+init=epsg:4467"))
  distance <- as.numeric(st_distance(start,stop))
  
  
  
}
rm(data, abun, data_start, data_stop,x_start,y_start,x_stop,y_stop,
   start,stop,distance)





