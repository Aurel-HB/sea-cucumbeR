# this script focus on the extraction of the substrate annotation
# the objective is to create a dataframe that compile substrate and merge with 
# sea cucumber observation

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

#############
#load data
#############
ici <- paste(here(),"/via3_data_exploration/Data/substrat/", sep="")

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

# filter only the substrate annotations 
data <- data[grep("substrat",data$metadata),]

time_start <- c()
time_stop <- c()
for (variable in 1:length(data$temporal_coordinates)) {
  test <- data$temporal_coordinates[variable]
  test <- strsplit(test,split = ",")
  test[[1]][1] <- substr(test[[1]][1],start = 2,
                                    stop = nchar(test[[1]][1]))
  test[[1]][2] <- substr(test[[1]][2],start = 1,
                         stop = nchar(test[[1]][2])-1)
  time_start <- c(time_start,as.numeric(test[[1]][1]))
  time_stop <- c(time_stop,as.numeric(test[[1]][2]))
}

data <- data.frame(data, time_start, time_stop)

# extract station ####
station <- c()
for (variable in 1:length(data$temporal_coordinates)) {
  test <- data$file_list[variable]
  test <- substr(test,start = 20, stop = 22)
  if (!is.na(as.numeric(test))){
    test <- as.numeric(test)
  }
  station <- c(station,test)
}

data <- data.frame(data, station)

# extract substrate ####
substrat <- c()
for (variable in 1:length(data$metadata)) {
  test <- data$metadata[variable]
  test <- substr(test,start = 7, stop = 16)
  substrat <- c(substrat,test)
}

data <- data.frame(data, substrat)


#############
data <- data[,7:10]

saveRDS(data,
        paste(ici,
              "data_substrat_2025.rds",
              sep=""))

data_substrate <- readRDS(paste(
  ici,
  "data_substrat_2025.rds",
  sep=""))

###########################
# convert time in distance
###########################
HOLOTVSPM2025 <- readRDS(paste(here(),
                               "/comparing_biomass_method/Data/03_HOLOTVSPM2025_summary.rds",
                               sep=""))
data_substrate$y_start <- NA
data_substrate$y_stop <- NA
for (stn in unique(data_substrate$station)){
  stn_info <- HOLOTVSPM2025 %>% filter(STN==as.character(stn))
  substrat_station <- data_substrate %>% filter(station==stn)
  for(segment in 1:nrow(substrat_station)){
    substrat_station$y_start[segment] <- substrat_station$time_start[
      segment]*stn_info$haul_distance[1]/stn_info$haul_duration[1]
    substrat_station$y_stop[segment] <- substrat_station$time_stop[
      segment]*stn_info$haul_distance[1]/stn_info$haul_duration[1]
  }
  data_substrate[grep(stn,data_substrate$station),] <- substrat_station
  rm(stn_info,substrat_station)
}

saveRDS(data_substrate,
        paste(ici,
              "data_substrat_2025.rds",
              sep=""))

#############################################################
#associate the substrate with the sea cucumbers annotations
#############################################################
data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

data_position$substrate <- NA
for (stn in unique(data_substrate$station)){
  data_temporary <- data_position %>% filter(station==as.character(stn))
  substrat_station <- data_substrate %>% filter(station==stn)
  substrat_station$time_start[grep(min(substrat_station$time_start),
                        substrat_station$time_start)] <- trunc(
                          min(substrat_station$time_start))
  
  for (annot in 1:nrow(data_temporary)){
    annot_time <- data_temporary$time[annot]
    for (segment in 1:nrow(substrat_station)){
      if (substrat_station$time_start[segment]<annot_time &&
          annot_time < substrat_station$time_stop[segment]){
        data_temporary$substrate[annot] <- substrat_station$substrat[segment]
      }
    }
  }
  data_position$substrate[
    grep(stn,data_position$station)] <- data_temporary$substrate
}


saveRDS(data_position,
        paste(ici,
              "data_position_substrat.rds",
              sep=""))
