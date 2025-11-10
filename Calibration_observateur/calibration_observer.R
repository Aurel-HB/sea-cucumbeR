# comparison observers count ####

# package
library(dplyr)
library(here)

path <- paste(here(),"/Calibration_observateur/Data/", sep="")

#list of 10 stations used for the calibration between observers
list_station <- c(104,106,125,157,161,147,96,131,163,191) 

#enter the initial of the observers
list_observer <- c("AHB", "LL", "RM")
#list_observer <- c("AHB", "RM")
#list_observer <- c("AHB", "LL")

#### Proto 1.0 ####
folder <- "annot_proto_1.0/"
# summarize the counts in a dataframe ####
data_count <- data.frame(list_station,0,0,0)
names(data_count) <- c("station",list_observer)

for (i_stn in 1:length(list_station)){
  for (observer in list_observer){
    data <- read.csv(paste(path,folder,"HOLO2025_STN",
                           list_station[i_stn],"_",observer,".csv", sep = ""), 
                  sep = ",", header = TRUE, skip = 9)
    number <- length(data[,1])
    data_count[i_stn,observer] <- number
  }
}
rm(data)

# start with begining stats ####

data_count$mean <- rowMeans(data_count[,2:(length(list_observer)+1)])
data_count$sd <- 0
data_count$lecture_indic <- 0
for (i_stn in 1:length(data_count$station)){
  data_count$sd[i_stn] <- sd(as.numeric((data_count[i_stn, 
                                                   1:length(list_observer)+1])))
  data_count$lecture_indic[i_stn] <- 
    max(as.numeric((data_count[i_stn, 
                               1:length(list_observer)+1]))) -
    min(as.numeric((data_count[i_stn, 
                               1:length(list_observer)+1])))
}

#lecture indicator <- difference between the highest/smallest numbers
#index <- lecture_indicator/mean * 100
data_count$index <- round((data_count$lecture_indic/data_count$mean)*100,1)
data_count1.0 <- data_count

#qunatile9 <- quantile(data_count$index, probs = 0.9)
error <- mean(data_count$index)

#### Proto 1.3 ####
folder <- "annot_proto_1.3/"
list_station <- c(104,106,125,157,161,147,96,131,163,191) 
# summarize the counts in a dataframe ####
data_count <- data.frame(list_station,0,0,0,0,0,0,0,0,0)
names(data_count) <- c("station",
                       paste(c("Point","Rectangle","Tot"),
                             list_observer[1] , sep = "_"),
                       paste(c("Point","Rectangle","Tot"),
                             list_observer[2] , sep = "_"),
                       paste(c("Point","Rectangle","Tot"),
                             list_observer[3] , sep = "_")
                       #paste(c("Rectangle"),list_observer , sep = "_"),
                       #paste(c("Tot"),list_observer , sep = "_")
                       )

for (i_stn in 1:length(list_station)){
  for (observer in list_observer){
    data <- read.csv(paste(path,folder,"HOLO2025_STN",
                           list_station[i_stn],"_",observer,".csv", sep = ""), 
                     sep = ",", header = TRUE, skip = 9)
    number <- length(data[,1])
    data_count[i_stn,paste("Tot",observer, sep = "_")] <- number
    shape <- substr(data$spatial_coordinates,start = 2, stop = 2)
    data_count[i_stn,
               paste("Point",observer, sep = "_")] <- length(
                 grep(pattern = "1", shape))
    data_count[i_stn,
               paste("Rectangle",observer, sep = "_")] <- length(
                 grep(pattern = "2", shape))
  }
}

# start with begining stats ####
data_count_rectangle <- data_count[,c(1,grep(pattern = "Rectangle",
                                             x = names(data_count)))]
data_count_point <- data_count[,c(1,grep(pattern = "Point",
                                             x = names(data_count)))]
data_count_tot <- data_count[,c(1,grep(pattern = "Tot",
                                         x = names(data_count)))]

table_stat <- function(data_count, list_observer){
  data_count$mean <- rowMeans(data_count[,2:(length(list_observer)+1)])
  data_count$sd <- 0
  data_count$lecture_indic <- 0
  for (i_stn in 1:length(data_count$station)){
    data_count$sd[i_stn] <- sd(as.numeric((data_count[i_stn, 
                                                      1:length(list_observer)+1])))
    data_count$lecture_indic[i_stn] <- 
      max(as.numeric((data_count[i_stn, 
                                 1:length(list_observer)+1]))) -
      min(as.numeric((data_count[i_stn, 
                                 1:length(list_observer)+1])))
  }
  
  #lecture indicator <- difference between the highest/smallest numbers
  #index <- lecture_indicator/mean * 100
  data_count$index <- round((data_count$lecture_indic/data_count$mean)*100,1)
  return(data_count)
}

data_count_rectangle <- table_stat(data_count_rectangle, list_observer)
data_count_point <- table_stat(data_count_point, list_observer)
data_count_tot <- table_stat(data_count_tot, list_observer)


qunatile9 <- quantile(data_count$index, probs = 0.9)
error <- mean(data_count_tot$index)

#### Proto 1.4 ####
# the station 106, 157, 163 do not have been count with the proto 1.4
list_station <- c(104,125,161,147,96,131,191)

data_count_1.4 <- data_count_tot %>% 
  filter(station %in% list_station)
error_1.4 <- mean(data_count_1.4$index)
