# This script prepare sf object with the sampling haul for 2021, 2022, 2023,
# 2025 surveys with density of sea cucumber, bottom temp, bathy and location

#############
#load packages
#############
library(here)
library(sf)
library(sp)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(XML)
library(viridis)

#############################################################
#PART 1 - Work on the  annotation summary of 2021, 2022, 2023
#############################################################

# import the annotation summary ###
count_tot_2021 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/HOLOSPMTV21.csv", sep=""), 
                      sep = ";", header = TRUE, skip = 12)
count_tot_2022 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/HOLOSPMTV22.csv", sep=""), 
                           sep = ";", header = TRUE, skip = 13)
count_tot_2023 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/HOLOSPMTV23.csv", sep=""), 
                           sep = ";", header = TRUE, skip = 13)


# import the feuille_route to have stop and start coordinates ###
feuille_route_2021 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/coordonnees_holospm2021.csv",
                                sep=""), sep = ";", header = TRUE,
                          fileEncoding='UTF-8', check.names=F)
feuille_route_2022 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/coordonnees_holospm2022.csv",
        sep=""), sep = ";", header = TRUE,
  fileEncoding='UTF-8', check.names=F)
feuille_route_2023 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/coordonnees_holospm2023.csv",
        sep=""), sep = ";", header = TRUE,
  fileEncoding='UTF-8', check.names=F)


#create dataframes with the good coordinates and the important information ####

list_count <- list("count_tot_2021"=count_tot_2021,
              "count_tot_2022"=count_tot_2022,
              "count_tot_2023"=count_tot_2023)
list_feuille_route <- list("feuille_route_2021"=feuille_route_2021,
                      "feuille_route_2022"=feuille_route_2022,
                      "feuille_route_2023"=feuille_route_2023)

for (indice in 1:3){
  feuille_route <- list_feuille_route[[names(list_feuille_route)[indice]]]  
  count <- list_count[[names(list_count)[indice]]]
  
  feuille_route_essential_start <- feuille_route[, c("Date",
                                                  "N° Station",
                                                  "Heure début Virage",
                                                  "Latitude début Virage (N)",
                                                  "Longitude début Virage (W)",
                                                  "Vitesse de traîne (nds)")]
  
  feuille_route_essential_stop <- feuille_route[, c("Date",
                                                    "N° Station",
                                                    "Heure fin Filage",
                                                    "Latitude fin Filage (N)",
                                                    "Longitude fin Filage (W)")]
  
  names(feuille_route_essential_start) <- c("Date", "Station", "heure",
                                            "latitude", "longitude","vitesse") 
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
  feuille_route_essential_start$latitude <- coord$y
  feuille_route_essential_start$longitude <- coord$x
  coord <- transfo_WG84_WGSPM(feuille_route_essential_stop[,c(2,5,4)])
  feuille_route_essential_stop$latitude <- coord$y
  feuille_route_essential_stop$longitude <- coord$x
  rm(coord)
  
  # transform the time of a haul in a number of seconds
  temps <- as.numeric(as.POSIXct(feuille_route_essential_start$heure,
                      format = ("%H:%M")) -
    as.POSIXct(feuille_route_essential_stop$heure, format = ("%H:%M")))
  temps <- temps*-60
  feuille_route_essential_start$temps <- temps
  
  #####
  # create a data with the total of sea cucumber per station and
  # the surface sampled
  #####
  data_position <- data.frame()
  for (stn in unique(count$STN)){
    Nombre <- count %>% filter(STN == stn)%>%
      select(STN,Nbre.adultes,taux.echant)
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
                      crs = 4467)
    stop <- st_as_sf(data.frame(x=x_stop, y=y_stop),
                     coords = c("x","y"),
                     crs = 4467)
    distance <- as.numeric(st_distance(start,stop))
    surface <- distance*1.5 #the field view is about 1.5m 
    
    # extract the centroid of the tracks
    X.track <- mean(c(x_start,x_stop))
    Y.track <- mean(c(y_start,y_stop))
    
    #extract the time of annotation start
    time_start <- as.POSIXlt(data_start[,3], format = "%H:%M:%S")
    time_start <- time_start$min*60 +time_start$sec
    
    # place the value in the table
    data_ligne <- data.frame(
      "station" = Nombre$STN[1],
      "abun" = Nombre$Nbre.adultes[1],
      "long" = X.track,
      "lat" = Y.track,
      "time" = data_start$temps[1],
      "distance" = distance,
      "speed_ms" = data_start$vitesse[1]*(1852/3600),
      "surface" = surface,
      "sampling.rate" = as.numeric(Nombre$taux.echant[1])
    )
    
    data_position <- rbind(data_position, data_ligne)
  }
  
  # there is a problem with the calculation of the distance using coordinates
  # we will approximate the distance with the speed and the time of sampling
  # when speed is NA -> mean(speed) and time is rate*600
  
  average_speed <- mean(data_position$speed_ms[!is.na(data_position$speed_ms)])
  data_position <- data_position %>%
    mutate(distance = ifelse(is.na(speed_ms),
                             average_speed*sampling.rate*600,
                             speed_ms*sampling.rate*600)) %>%
    mutate(surface = distance*1.5)
  assign(paste("data_position",2020+indice,sep = "_"), data_position)
}

saveRDS(data_position,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_2021.rds",
              sep=""))
saveRDS(data_position,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_2022.rds",
              sep=""))
saveRDS(data_position,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_2023.rds",
              sep=""))

#############################################################
#PART 2 - Create a data table with all the year
#############################################################

# load data
data_abun_2025 <- readRDS(
  paste(here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
        sep=""))
centroid_2025 <- read.csv(paste(
  here(),"/SIG/SIG_data/centroid_tracks.csv", sep=""))

data_position_2021 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_2021.rds",
        sep="")
)

data_position_2022 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_2022.rds",
        sep="")
)

data_position_2023 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_2023.rds",
        sep="")
)

# create a data frame with station,X,Y,abun,intensity,area for each year ####
# 2025
for (i in 1:nrow(data_abun_2025)){
  coords <-centroid_2025 %>% filter(Station == data_abun_2025$station[i])
  data_abun_2025$X[i] <- coords$X.track[1]
  data_abun_2025$Y[i] <- coords$Y.track[1]
}

data_abun_2025 <- data_abun_2025 %>%
  as.data.frame()%>%
  select(station,X,Y,abun,intensity,area)

# 2021
data_abun_2021 <- data_position_2021 %>%
  mutate(intensity = abun/surface) %>%
  select(station,long,lat,abun,intensity,surface)
names(data_abun_2021) <- c("station","X","Y","abun","intensity","area")

# 2022
data_abun_2022 <- data_position_2022 %>%
  mutate(intensity = abun/surface) %>%
  select(station,long,lat,abun,intensity,surface)
names(data_abun_2022) <- c("station","X","Y","abun","intensity","area")

# 2023
data_abun_2023 <- data_position_2023 %>%
  mutate(intensity = abun/surface) %>%
  select(station,long,lat,abun,intensity,surface)
names(data_abun_2023) <- c("station","X","Y","abun","intensity","area")


# create a sf that have all the abun data over the year ####
data_abun_2021$year <- 2021
data_abun_2022$year <- 2022
data_abun_2023$year <- 2023
data_abun_2025$year <- 2025

data_abun_tot <- rbind(data_abun_2021,
                       data_abun_2022,
                       data_abun_2023,
                       data_abun_2025)

# convert in sf
data_abun_tot <- data_abun_tot %>%
  mutate(long = X) %>%
  mutate(lat = Y) %>%
  st_as_sf(.,coords = c("long","lat"),crs = 4467)

saveRDS(data_abun_tot,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_abun_tot.rds",
              sep=""))

#############################################################
#PART 3 - Add covariable to the sf table with all the abundance
#############################################################

#load data
bathy_spm <- readRDS(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
                           "Bathy_3PS.rds",sep = ""))
sea_floor_temp_spm <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/",
        "BottomT_500.rds",sep = ""))
sea_floor_temp_spm <- sea_floor_temp_spm[,grep(pattern = "May",
                                               names(sea_floor_temp_spm))]
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

data_abun <- list(
  "2021" = data_abun_2021,
  "2022" = data_abun_2022,
  "2023" = data_abun_2023,
  "2025" = data_abun_2025
)

for (annee in c(2021,2022,2023,2025)){
  data <- sea_floor_temp_spm[,grep(pattern = as.character(annee),
                                   names(sea_floor_temp_spm))]
  names(data) <- c("temp","geometry")
  data_abun[[as.character(annee)]] <- data_abun[[as.character(annee)]] %>%
    mutate(long = X) %>%
    mutate(lat = Y) %>%
    st_as_sf(.,coords = c("long","lat"),crs = 4467) 
  
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              bathy_spm)
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              data)
}


data_abun_tot_cov <- rbind(data_abun[["2021"]],
                           data_abun[["2022"]],
                           data_abun[["2023"]],
                           data_abun[["2025"]])

# keep only point in the study area
data_abun_tot_cov <- st_join(calcul_area,data_abun_tot_cov)
data_abun_tot_cov <- as.data.frame(data_abun_tot_cov)
data_abun_tot_cov <- data_abun_tot_cov[,1:(length(names(data_abun_tot_cov))-1)]
data_abun_tot_cov <- data_abun_tot_cov %>%
  mutate(long = X) %>%
  mutate(lat = Y) %>%
  st_as_sf(.,coords = c("long","lat"),crs = 4467)


saveRDS(data_abun_tot_cov,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_abun_tot_cov.rds",
              sep=""))

ggplot(data_abun_tot_cov)+
  geom_sf(aes(color=temp))+
  scale_color_viridis()+
  geom_sf(data=calcul_area, fill = "#11111111")+
  facet_wrap(~year, ncol=2)+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0d1e6"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))

ggplot(data_abun_tot_cov)+
  geom_sf(aes(color=log(intensity)), size = 3)+
  scale_color_viridis()+
  geom_sf(data=calcul_area, fill = "#11111111")+
  facet_wrap(~year, ncol=2)+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0d1e6"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))
