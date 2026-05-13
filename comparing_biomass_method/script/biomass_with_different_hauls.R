# This code work on the decrease of information by shorter dredge haul

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(boot)


#############################################################
# work on simulated data
#############################################################

shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))
data_position <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_position.rds",
        sep=""))
data_abun <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_abun.rds", sep=""))
predata_intensity <- readRDS(
  paste(here(), "/point_process/Output/holo_simu_intensity.rds", sep=""))


data<- data.frame()
abun <- c()
reduction <- 10#0#2#4#10
for (i in unique(data_position$station)){
  data_used <- data_position[data_position$station == i,]
  distance <- max(data_used$lat) - min(data_used$lat)
  threshold <- max(data_used$lat) - distance/reduction
  data_used <- data_used %>% filter(lat > threshold)
  abun <- c(abun, nrow(data_used))
  data <- rbind(data, data_used)
}

data_abun$abun <- abun


# Calculate biomass with extrapolation of density of haul to square #####


# param for real data
vitesseNavire <- 1.4 # (in mile/h) can be changed by a vector if available
intervalTemps <- 10 # (in min)  can be changed by a vector if available
champVision <- 1.5/1000 # (in kilometer)
surfCompte <- as.numeric(intervalTemps*champVision*vitesseNavire*1852/60)
mean_weight <- 0.4 # (in kg)
surfsquare <- 1.852*1.852 # in kilometer square

# param for simulated data 
surfCompte <- 0.5*champVision/reduction # in kilometer square

## Prepare the dataframe
dat <- as.data.frame(
  cbind(data_abun,as.data.frame(predata_intensity[,c("X","Y")])[,1:2]))

## Prepare grid for prediction used in the sdmtmb approach
## Enable to have the same number of point in the total stock of the 2 method
# use the sampling grid to generate a extrapolation grid
grid <- shp_grid %>% 
  dplyr::filter(Id > 85) %>%
  dplyr::filter(Id != 188) %>%
  dplyr::filter(Id != 189) %>%
  st_centroid() %>%
  mutate(long = st_coordinates(.)[,1]/1000) %>%
  mutate(lat = st_coordinates(.)[,2]/1000) %>%
  as.data.frame() %>%
  select(long,lat) 

plot(grid)


## Analyse per zone
## number of cucumber relative to total area
## Density = x/y with x = number individual et y = surface counted
dat$density <- dat$abun/surfCompte

## RESULT ABUNDANCE --------------------------------------
## abundance in number of individuals relative to the total area of each square
dat$Abun_square <- dat$density*surfsquare

## RESULT BIOMASSE --------------------------------------
## biomass in tonnes relative to the total area of each zone
dat$bioAdult <- dat$Abun_square*mean_weight/1000
## total biomass of the stock
bio_stock <- mean(dat$bioAdult)*length(grid$long)
## --------------------------------------------------------

## Estimation of sampling variances ####
### Variances with bootstrap
# N = number of replica to define
N <- 10000

# function to calculate the total biomass
bio_tot <- function(donnees, indices) {
  abondance_echantillon <- donnees[indices]
  abondance_m2 <- abondance_echantillon / surfCompte
  biomasse_carre <- abondance_m2 * surfsquare * mean_weight/1000 # in ton
  biomasse_totale <- sum(biomasse_carre)*length(grid$long)/
    length(dat$abun) # sampling on the grid 
  return(biomasse_totale)
}

# execute the bootstrap
set.seed(13) # for reproducibility
bootobject <- boot(data = dat$abun, statistic = bio_tot, R = N)

## FIGURE : ABUNDANCE ADULTS
# used statistics
boot_bio <- as.numeric(bootobject[["t"]])
median_value <- median(boot_bio)
quantile_5 <- quantile(boot_bio, 0.05)
quantile_95 <- quantile(boot_bio, 0.95)
mean_value <- mean(boot_bio)
data_stat <- data.frame(
  reduction = reduction,
  median_value = median(boot_bio),
  quantile_5 = quantile(boot_bio, 0.05),
  quantile_95 = quantile(boot_bio, 0.95),
  mean_value = mean(boot_bio)
)

# Créer l'histogramme avec ggplot2
ggplot(data.frame(Biomass = boot_bio), aes(x = Biomass)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", color = "black") +
  geom_density(col='red') +
  geom_vline(aes(xintercept = median_value), color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile_5), color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = quantile_95), color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = mean_value), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Histogramme des Biomasses obtenues par Bootstrap",
       x = "Biomasse",
       y = "Densité") +
  scale_x_continuous(labels = scales::label_comma())+
  theme_minimal()

# keep result of simulation 
# after every simulation just add the stat in the table

#data_compare <- data.frame() # run only the first time
data_compare <- rbind(data_compare, data_stat)
saveRDS(data_compare,
        file = paste(here(),
                     "/point_process/Output/holo_simu_reduction_comparison.rds",
                     sep=""))

#############################################################
# work on real data
#############################################################

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
  station <- c(station,as.numeric(test))
}

data <- data.frame(data, station)

#############
data <- data[,7:11]

saveRDS(data, paste(here(),
                    "/via3_data_exploration/Data/processed/data_position_2025",
                    sep=""))

data <- readRDS(paste(here(),
                      "/via3_data_exploration/Data/processed/data_position_2025",
                      sep=""))
#############
#test with high, medium and low density (133, 98 and 106)
#############
data.low <- data %>% filter(station==106)
data.medium <- data %>% filter(station==98)

#### low density ####
data.compare.low <- data.frame(reduction = 0, abun = nrow(data.low))

data.compare.low <- rbind(data.compare.low,
                          data.frame(
                            reduction = c(2,2),
                            abun = c(
                              nrow(data.low %>% filter(time>300))*2,
                              nrow(data.low %>% filter(time<=300))*2
                            )
                          ))

#### medium density ####
data.compare.medium <- data.frame(reduction = 0, abun = nrow(data.medium))

# reduction per 2
reduction <- 2
start <- seq(0,600-300,30)
stop <- start + 300
abun <- c()
for (i in 1:length(start)){
  value <- nrow(data.medium %>% 
                  filter(time > start[i]) %>%
                  filter(time <= stop[i]))
  abun <- c(abun,value*2)
}

data.compare.medium <- rbind(data.compare.medium,
                             data.frame(
                               reduction = reduction,
                               abun = abun))

# used statistics
stat_medium_2 <- data.compare.medium$abun[data.compare.medium$reduction==2]
data_stat_2 <- data.frame(
  reduction = reduction,
  median_value = median(stat_medium_2),
  quantile_5 = quantile(stat_medium_2, 0.05),
  quantile_95 = quantile(stat_medium_2, 0.95),
  mean_value = mean(stat_medium_2),
  sd = sd(stat_medium_2)
)

# reduction per 4
reduction <- 4
start <- seq(0,600-150,1)
stop <- start + 150
abun <- c()
for (i in 1:length(start)){
  value <- nrow(data.medium %>% 
                  filter(time > start[i]) %>%
                  filter(time <= stop[i]))
  abun <- c(abun,value*4)
}

data.compare.medium <- rbind(data.compare.medium,
                             data.frame(
                               reduction = reduction,
                               abun = abun))

# used statistics
stat_medium_4 <- data.compare.medium$abun[data.compare.medium$reduction==4]
data_stat_4 <- data.frame(
  reduction = reduction,
  median_value = median(stat_medium_4),
  quantile_5 = quantile(stat_medium_4, 0.05),
  quantile_95 = quantile(stat_medium_4, 0.95),
  mean_value = mean(stat_medium_4),
  sd = sd(stat_medium_4)
)

data_stat <- rbind(data_stat_2,data_stat_4)

#### high density ####
high.station <- c(114,
                  133,
                  149,
                  167)
data_stat <- data.frame()
for (indice in high.station){
  data.high <- data %>% filter(station==indice)
  
  data.compare.high <- data.frame(reduction = 0, abun = nrow(data.high))
  
  # reduction per 2
  reduction <- 2
  start <- seq(round(min(data.high$time)/60, digits = 0)*60,
               600-300,30)
  stop <- start + 300
  abun <- c()
  for (i in 1:length(start)){
    value <- nrow(data.high %>% 
                    filter(time > start[i]) %>%
                    filter(time <= stop[i]))
    abun <- c(abun,value*2)
  }
  
  data.compare.high <- rbind(data.compare.high,
                             data.frame(
                               reduction = reduction,
                               abun = abun))
  
  # used statistics
  stat_high_2 <- data.compare.high$abun[data.compare.high$reduction==2]
  data_stat_2 <- data.frame(
    station = indice,
    reduction = reduction,
    median_value = median(stat_high_2),
    quantile_5 = quantile(stat_high_2, 0.05),
    quantile_95 = quantile(stat_high_2, 0.95),
    mean_value = mean(stat_high_2),
    sd = sd(stat_high_2)
  )
  
  # reduction per 4
  reduction <- 4
  start <- seq(round(min(data.high$time)/60, digits = 0)*60,
               600-150,10)
  stop <- start + 150
  abun <- c()
  for (i in 1:length(start)){
    value <- nrow(data.high %>% 
                    filter(time > start[i]) %>%
                    filter(time <= stop[i]))
    abun <- c(abun,value*4)
  }
  
  data.compare.high <- rbind(data.compare.high,
                             data.frame(
                               reduction = reduction,
                               abun = abun))
  
  # used statistics
  stat_high_4 <- data.compare.high$abun[data.compare.high$reduction==4]
  data_stat_4 <- data.frame(
    station = indice,
    reduction = reduction,
    median_value = median(stat_high_4),
    quantile_5 = quantile(stat_high_4, 0.05),
    quantile_95 = quantile(stat_high_4, 0.95),
    mean_value = mean(stat_high_4),
    sd = sd(stat_high_4)
  )
  
  
  # reduction per 10
  reduction <- 10
  start <- seq(round(min(data.high$time)/60, digits = 0)*60,
               600-60,10)
  stop <- start + 60
  abun <- c()
  for (i in 1:length(start)){
    value <- nrow(data.high %>% 
                    filter(time > start[i]) %>%
                    filter(time <= stop[i]))
    abun <- c(abun,value*10)
  }
  
  data.compare.high <- rbind(data.compare.high,
                             data.frame(
                               reduction = reduction,
                               abun = abun))
  
  # used statistics
  stat_high_10 <- data.compare.high$abun[data.compare.high$reduction==10]
  data_stat_10 <- data.frame(
    station = indice,
    reduction = reduction,
    median_value = median(stat_high_10),
    quantile_5 = quantile(stat_high_10, 0.05),
    quantile_95 = quantile(stat_high_10, 0.95),
    mean_value = mean(stat_high_10),
    sd = sd(stat_high_10)
  )
  
  data_stat <- rbind(data_stat,data_stat_2,data_stat_4,data_stat_10)
  
}
