#library
library(here)
library(stringr)
library(dplyr)

#path
path_input <- paste(here(),
                    "/survey_data_exploration/anotate_Data/row_data/", sep="")
path_output <- paste(here(),
                     "/survey_data_exploration/anotate_Data/process_data/",
                     sep="")



info <- read.csv2(paste(path_input,"Info_station_2021.csv", sep = ""))


name_video <- "HOLO21_STN079"
#start process on 1 data ####
data.row <- read.csv(paste(path_input,name_video,".csv", sep = ""), sep = "," )

data <- data.row[2:length(data.row[,1]),]

data$temps <- substr(data$X2..Video.or.Image.Identifier, start = 1, stop = 8)
data$temps <- as.POSIXct(data$temps, format = "%H:%M:%S")

station <- as.numeric(str_sub(string = name_video,start = -3,end = -1))
for (i in (1:length(info[,1]))){
  if (info$STN[i] == station){
    start_tps <- info$Start[i]
    end_tps <- info$Stop[i]
  }
}

start_tps <- as.POSIXct(start_tps, format = "%M:%S")
end_tps <- as.POSIXct(end_tps, format = "%M:%S")

data.filter <- data %>% dplyr::filter(temps>= start_tps) %>%
  dplyr::filter(temps< end_tps) %>%
  mutate(minute = as.numeric(substr(as.character(temps),start = 15,stop = 16)))

for (indice in (1:length(data.filter[,1]))) {
  start_tps <- data.filter$minute[1]
  minute <- data.filter$minute[indice]
  time <- data.filter$X2..Video.or.Image.Identifier[indice]
  new_minute <- minute - start_tps
  if (minute<10){
    new_time <- paste(substr(time,start = 1,stop = 4),
                      new_minute, str_sub(string = time,
                                          start = 6,end = nchar(time)), sep="")
    data.filter$X2..Video.or.Image.Identifier[indice] <- new_time
  }
  if (minute>=10){
    new_time <- paste(substr(time,start = 1,stop = 3),"0",
                      new_minute, str_sub(string = time,
                                          start = 6,end = nchar(time)), sep="")
    data.filter$X2..Video.or.Image.Identifier[indice] <- new_time
  }
}

start_frame <- as.numeric(
  data.filter$X3..Unique.Frame.Identifier[1]
  ) - as.numeric(
  str_sub(string = data.filter$X3..Unique.Frame.Identifier[1],
          start=-2, end=-1))
                      

data.filter <- data.filter %>% 
  mutate(X3..Unique.Frame.Identifier = as.numeric(
    X3..Unique.Frame.Identifier) - start_frame)



data.process <- rbind(data.row[1,], data.filter[,1:dim(data.row)[2]])
# export csv ###
write.csv(data.process, paste(path_output,name_video,".csv", sep =""), row.names=FALSE)


#process for all the data ####
# list the name of the files
names_video <- c("HOLO21_STN081","HOLO21_STN120","HOLO21_STN137","HOLO21_STN161")

for (name_video in names_video) {
  #start process on 1 data ####
  data.row <- read.csv(paste(path_input,name_video,".csv", sep = ""), sep = "," )
  
  data <- data.row[2:length(data.row[,1]),]
  
  data$temps <- substr(data$X2..Video.or.Image.Identifier, start = 1, stop = 8)
  data$temps <- as.POSIXct(data$temps, format = "%H:%M:%S")
  
  station <- as.numeric(str_sub(string = name_video,start = -3,end = -1))
  for (i in (1:length(info[,1]))){
    if (info$STN[i] == station){
      start_tps <- info$Start[i]
      end_tps <- info$Stop[i]
    }
  }
  
  start_tps <- as.POSIXct(start_tps, format = "%M:%S")
  end_tps <- as.POSIXct(end_tps, format = "%M:%S")
  
  data.filter <- data %>% dplyr::filter(temps>= start_tps) %>%
    dplyr::filter(temps< end_tps)%>%
    mutate(minute = as.numeric(substr(as.character(temps),start = 15,stop = 16)))
  
  for (indice in (1:length(data.filter[,1]))) {
    start_tps <- data.filter$minute[1]
    minute <- data.filter$minute[indice]
    time <- data.filter$X2..Video.or.Image.Identifier[indice]
    new_minute <- minute - start_tps
    if (minute<10){
      new_time <- paste(substr(time,start = 1,stop = 4),
                        new_minute, str_sub(string = time,
                                            start = 6,end = nchar(time)), sep="")
      data.filter$X2..Video.or.Image.Identifier[indice] <- new_time
    }
    if (minute>=10){
      new_time <- paste(substr(time,start = 1,stop = 3),"0",
                        new_minute, str_sub(string = time,
                                            start = 6,end = nchar(time)), sep="")
      data.filter$X2..Video.or.Image.Identifier[indice] <- new_time
    }
  }
  
  start_frame <- as.numeric(
    data.filter$X3..Unique.Frame.Identifier[1]
  ) - as.numeric(
    str_sub(string = data.filter$X3..Unique.Frame.Identifier[1],
            start=-2, end=-1))
  
  
  data.filter <- data.filter %>% 
    mutate(X3..Unique.Frame.Identifier = as.numeric(
      X3..Unique.Frame.Identifier) - start_frame)
  
  data.process <- rbind(data.row[1,], data.filter[,1:dim(data.row)[2]])
  # export csv ###
  write.csv(data.process, paste(path_output,name_video,".csv", sep =""), row.names=FALSE)
  
}
