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
#library(spatstat)
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

saveRDS(data, paste(here(),
                    "/via3_data_exploration/Data/processed/data_position_2025",
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

#transform the coordinates to use in sf
deg_dec <- function(coord){
  # coords is the shape xx°xx°xxx 
  dec <- as.numeric(substr(coord, start = 1, stop = 2))
  min <- as.numeric(substr(coord, start = 4, stop = 5))/60
  sec <- as.numeric(substr(coord, start = 7, stop = nchar(coord)))/3600
  return(dec+min+sec)
}

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


transfo_WG84_WGSPM = function(data,
                              src.proj = CRS("+init=epsg:4326"),
                              dst.proj = CRS("+init=epsg:4467")) {
  require(sp)
  #data is the shape (id,lon,lat) in this order
  names(data) <- c("id","lon","lat")
  as.data.frame(
    spTransform(
      SpatialPointsDataFrame(
        coords = data.frame(Xbng = data$lon,
                            Ybng = data$lat),
        data = data.frame(id = data$id,
                          Xlon = data$lon,
                          Ylat = data$lat),
        proj4string = src.proj), dst.proj))
  
}

coord <- transfo_WG84_WGSPM(feuille_route_essential_start[,c(2,5,4)])
feuille_route_essential_start$latitude <- coord$coords.x2
feuille_route_essential_start$longitude <- coord$coords.x1
coord <- transfo_WG84_WGSPM(feuille_route_essential_stop[,c(2,5,4)])
feuille_route_essential_stop$latitude <- coord$coords.x2
feuille_route_essential_stop$longitude <- coord$coords.x1
rm(coord)



