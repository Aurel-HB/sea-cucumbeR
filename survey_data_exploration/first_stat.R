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

name_video <- "HOLO22_STN147"
detect <- read.csv(paste(path_input,name_video,"_detections.csv", sep = ""), 
                   sep = "," )
correct <- read.csv(paste(path_input,name_video,"_detections_corrected.csv", 
                          sep = ""))
