# This code transform the Confidence Threshold of the correction .csv to 1

#library
library(here)
library(stringr)
library(dplyr)

#path
path_input <- paste(here(),
                    "/viame_data_exploration/anotate_Data/row_data/", sep="")
path_output <- paste(here(),
                     "/viame_data_exploration/anotate_Data/process_data/",
                     sep="")



#process for all the video ####
# list the name of the files
names_video <- c("HOLO22_STN112","HOLO22_STN147","HOLO22_STN176","HOLO22_STN205")


for (name_video in names_video) {
  row_correct <- read.csv(paste(path_input,name_video,"_detections_corrected.csv", 
                                sep = ""))
  
  correct <- row_correct[2:length(row_correct[,1]),]
  correct$X3..Unique.Frame.Identifier <- as.numeric(
    correct$X3..Unique.Frame.Identifier
  )
  correct$X..1..Detection.or.Track.id <- as.numeric(
    correct$X..1..Detection.or.Track.id
  )
  
  correct <- correct %>% dplyr::mutate(X8..Detection.or.Length.Confidence = 1)
  # export csv ###
  write.csv(data.process, paste(path_output,name_video,
                                "_detections_corrected.csv", sep =""),
            row.names=FALSE)
}