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
  dplyr::filter(temps<= end_tps)

data.process <- rbind(data.row[1,], data.filter[,1:dim(data.row)[2]])
# export csv ###
write.csv(data.process, paste(path_output,name_video,".csv", sep =""), row.names=FALSE)
