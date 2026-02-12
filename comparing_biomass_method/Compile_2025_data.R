# This code creates a compilation of all data from the 2025 survey.

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


## Prepare the dataframe to compil the information ####
dat <-  merge(data_annot,data_coordinates,by="STN")
#calculate surface sampled for density
dat$annotation_distance <- dat$annotation_rate*dat$haul_distance
dat$surface <- dat$annotation_distance*field_view

#convert the time of annotation start in second
dat$Time.annotation_start <- as.POSIXlt(dat$Time.annotation_start,
                                        format = "%H:%M:%S")
dat$Time.annotation_start <- dat$Time.annotation_start$min*60 +
  dat$Time.annotation_start$sec

dat$Distance.annotation_start <- (dat$Time.annotation_start/dat$haul_duration)*
  dat$haul_distance

dat$biomass_density <- dat$Total_Number*mean_weight/(
  dat$surface)#/1e+06) #density in kg/km2

dat <- dat%>% select(
  c(1,2,3,4,16,26,5,6,7,17,27,29,10,30,28,11,9,8,12,13,14,15,18,19,20,21,22,23,24,25)
)

saveRDS(dat, 
        paste(here(),
              "/comparing_biomass_method/Data/03_HOLOTVSPM2025_summary.rds",
        sep=""))
