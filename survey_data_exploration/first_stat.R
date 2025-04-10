# This code compare the detection and the correction to determine the quality 
# of the AI model

#library
library(here)
library(stringr)
library(dplyr)
library(ggplot2)
library(viridis)

#path
path_input <- paste(here(),
                    "/survey_data_exploration/anotate_Data/row_data/", sep="")
path_output <- paste(here(),
                     "/survey_data_exploration/anotate_Data/process_data/",
                     sep="")

# test work for one video ####
name_video <- "HOLO22_STN112"
row_detect <- read.csv(paste(path_input,name_video,"_detections.csv", sep = ""), 
                   sep = "," )
row_correct <- read.csv(paste(path_input,name_video,"_detections_corrected.csv", 
                          sep = ""))

detect <- row_detect[2:length(row_detect[,1]),]
correct <- row_correct[2:length(row_correct[,1]),]
correct$X3..Unique.Frame.Identifier <- as.numeric(
  correct$X3..Unique.Frame.Identifier
)
correct$X..1..Detection.or.Track.id <- as.numeric(
  correct$X..1..Detection.or.Track.id
)
detect$X..1..Detection.or.Track.id <- as.numeric(
  detect$X..1..Detection.or.Track.id
)

cucumaria.detect <- detect[grep(
  "cucumaria",detect$X10.11...Repeated.Species),]
cucumaria.correct <- correct[grep(
  "cucumaria",correct$X10.11...Repeated.Species),]

# check if there tracking in the corrected file
duplica <- {correct} %>%
  dplyr::group_by(X..1..Detection.or.Track.id
  ) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L)

#test1 <- c()
false_negative <- 0
true_positive <- 0 
for (i in cucumaria.correct$X..1..Detection.or.Track.id){
  if(i %in% cucumaria.detect$X..1..Detection.or.Track.id){
    true_positive <- true_positive + 1
    #test1 <- c(i,test1)
  } else {
    false_negative <- false_negative + 1
  }
}

#test2 <- c()
false_positive <- 0
#true_positive <- 0 
for (i in cucumaria.detect$X..1..Detection.or.Track.id){
  if(i %in% cucumaria.correct$X..1..Detection.or.Track.id){
    #true_positive <- true_positive + 1
    #test2 <- c(i,test2)
  } else {
    false_positive <- false_positive + 1
  }
}


#process for all the video ####
# list the name of the files
names_video <- c("HOLO22_STN112","HOLO22_STN147","HOLO22_STN176","HOLO22_STN205")
indicator <- data.frame(video_ID = names_video,
                        tot_detect = 0,
                        tot_correct = 0,
                        false_positive = 0,
                        false_negative = 0,
                        true_positive = 0,
                        true_positive_prop = 0,
                        false_positive_prop = 0,
                        false_negative_prop = 0)

for (name_video in names_video) {
  row_detect <- read.csv(paste(path_input,name_video,"_detections.csv", sep = ""), 
                         sep = "," )
  row_correct <- read.csv(paste(path_input,name_video,"_detections_corrected.csv", 
                                sep = ""))
  
  detect <- row_detect[2:length(row_detect[,1]),]
  correct <- row_correct[2:length(row_correct[,1]),]
  correct$X3..Unique.Frame.Identifier <- as.numeric(
    correct$X3..Unique.Frame.Identifier
  )
  correct$X..1..Detection.or.Track.id <- as.numeric(
    correct$X..1..Detection.or.Track.id
  )
  detect$X..1..Detection.or.Track.id <- as.numeric(
    detect$X..1..Detection.or.Track.id
  )
  
  cucumaria.detect <- detect[grep(
    "cucumaria",detect$X10.11...Repeated.Species),]
  cucumaria.correct <- correct[grep(
    "cucumaria",correct$X10.11...Repeated.Species),]
  
  false_negative <- 0
  true_positive <- 0 
  for (i in cucumaria.correct$X..1..Detection.or.Track.id){
    if(i %in% cucumaria.detect$X..1..Detection.or.Track.id){
      true_positive <- true_positive + 1
    } else {
      false_negative <- false_negative + 1
    }
  }
  
  
  false_positive <- 0
  for (i in cucumaria.detect$X..1..Detection.or.Track.id){
    if(i %in% cucumaria.correct$X..1..Detection.or.Track.id){
    } else {
      false_positive <- false_positive + 1
    }
  }
  
  tot_detect <- length(cucumaria.detect$X..1..Detection.or.Track.id)
  tot_correct <- length(cucumaria.correct$X..1..Detection.or.Track.id)
  
  
  for (i in 1:length(indicator$video_ID)){
    if (indicator$video_ID[i] == name_video){
      indicator$tot_detect[i] <- tot_detect
      indicator$tot_correct[i] <- tot_correct
      indicator$false_positive[i] <- false_positive
      indicator$false_negative[i] <- false_negative
      indicator$true_positive[i] <- true_positive
      indicator$false_positive_prop[i] <- (false_positive/tot_detect)*100
      indicator$true_positive_prop[i] <- (true_positive/tot_detect)*100
      indicator$false_negative_prop[i] <- (false_negative/tot_correct)*100
    }
  }
  
}

#graph to see result ####
ggplot(data = indicator)+
  geom_point(aes(video_ID,true_positive_prop), colour = "tomato3", shape = 17)+
  geom_point(aes(video_ID,false_positive_prop), colour = "steelblue3", 
             shape = 18)+
  geom_point(aes(video_ID,false_negative), colour = "seagreen3", shape = 19)+
  geom_line(aes(1:4,true_positive_prop), colour = "tomato3")+
  geom_line(aes(1:4,false_positive_prop), colour = "steelblue3")+
  geom_line(aes(1:4,false_negative), colour = "seagreen3")+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="")+
  theme_classic()

ggplot(data = indicator)+
  geom_point(aes(video_ID,true_positive_prop), colour = "tomato3", shape = 17)+
  geom_point(aes(video_ID,false_positive_prop), colour = "steelblue3", 
             shape = 15)+
  geom_point(aes(video_ID,false_negative_prop), colour = "seagreen3", shape = 19)+
  geom_line(aes(1:4,true_positive_prop), colour = "tomato3")+
  geom_line(aes(1:4,false_positive_prop), colour = "steelblue3")+
  geom_line(aes(1:4,false_negative_prop), colour = "seagreen3")+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="")+
  theme_classic()

indic2 <- rbind(
  data.frame(video_ID = indicator$video_ID, 
             value = indicator$true_positive_prop,
             type = "true_positive_proportion"),
  data.frame(video_ID = indicator$video_ID, 
             value = indicator$false_positive_prop,
             type = "false_positive_proportion"),
  data.frame(video_ID = indicator$video_ID, 
             value = indicator$false_negative_prop,
             type = "false_negative_proportion")
  )

ggplot(data = indic2)+
  geom_point(aes(x = video_ID, y = value, colour = type, shape = type))+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="",
    fill="")+
  theme_classic()


ggplot(data = indic2)+
  geom_col(aes(x=video_ID, y = value, fill = type), position = "dodge")+
  scale_fill_discrete(type = c("seagreen3","steelblue3","tomato3") )+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="",
    fill="")+
  theme_bw()

ggplot(data = indic2 %>% dplyr::filter(type != "false_negative_proportion"))+
  geom_col(aes(x=video_ID, y = value, fill = type), position = "dodge")+
  scale_fill_discrete(type = c("steelblue3","tomato3") )+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="",
    fill="")+
  theme_bw()

ggplot(data = indic2 %>% dplyr::filter(type == "false_negative_proportion"))+
  geom_col(aes(x=video_ID, y = value, fill = type), position = "dodge")+
  scale_fill_discrete(type = c("seagreen3") )+
  labs(
    title = "Summary of the AI model's indicators",
    x="",
    y="",
    fill="")+
  theme_bw()
