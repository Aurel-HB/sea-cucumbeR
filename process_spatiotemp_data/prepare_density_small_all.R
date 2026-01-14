# This script prepare sf object with the sampling haul for 2021, 2022, 2023,
# 2025 surveys with density of small sea cucumber, bottom temp, bathy and location

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
  paste(here(),"/process_spatiotemp_data/Data/raw/JUV_HOLOSPMTV21.csv", sep=""), 
                      sep = ";", header = TRUE, skip = 12)
count_tot_2022 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/JUV_HOLOSPMTV22.csv", sep=""), 
                           sep = ";", header = TRUE, skip = 13)
count_tot_2023 <- read.csv(
  paste(here(),"/process_spatiotemp_data/Data/raw/JUV_HOLOSPMTV23.csv", sep=""), 
                           sep = ";", header = TRUE, skip = 13)


# import the data_position of the big sea cucumber to only change the abundance ###
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


#create dataframes by replacing the adult abun by the small ones abun ####
for (i_row in 1:nrow(data_position_2021)){
  STN <- data_position_2021$station[i_row]
  data_position_2021$abun[i_row] <- count_tot_2021[
    count_tot_2021$STN==STN,c("Nbre.petits")]
}

for (i_row in 1:nrow(data_position_2022)){
  STN <- data_position_2022$station[i_row]
  data_position_2022$abun[i_row] <- count_tot_2022[
    count_tot_2022$STN==STN,c("Nbre.petits")]
}

for (i_row in 1:nrow(data_position_2023)){
  STN <- data_position_2023$station[i_row]
  data_position_2023$abun[i_row] <- count_tot_2023[
    count_tot_2023$STN==STN,c("Nbre.petits")]
}


saveRDS(data_position_2021,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_small_2021.rds",
              sep=""))
saveRDS(data_position_2022,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_small_2022.rds",
              sep=""))
saveRDS(data_position_2023,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_position_small_2023.rds",
              sep=""))

#############################################################
#PART 2 - Create a data table with all the year
#############################################################

# load data
data_position_2021 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_small_2021.rds",
        sep="")
)

data_position_2022 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_small_2022.rds",
        sep="")
)

data_position_2023 <- readRDS(
  paste(here(),
        "/process_spatiotemp_data/Data/processed/data_position_small_2023.rds",
        sep="")
)

# create a data frame with station,X,Y,abun,intensity,area for each year ####

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

data_abun_tot <- rbind(data_abun_2021,
                       data_abun_2022,
                       data_abun_2023)

# convert in sf
data_abun_tot <- data_abun_tot %>%
  mutate(long = X) %>%
  mutate(lat = Y) %>%
  st_as_sf(.,coords = c("long","lat"),crs = 4467)

saveRDS(data_abun_tot,
        paste(here(),
              "/process_spatiotemp_data/Data/processed/data_abun_tot_small.rds",
              sep=""))

data_abun_tot <- readRDS(paste(here(),
              "/process_spatiotemp_data/Data/processed/data_abun_tot_small.rds",
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
chla_spm <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/",
        "Chla_500.rds",sep = ""))
chla_spm <- chla_spm[,grep(pattern = "May",names(chla_spm))]

vo_spm <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/",
        "Current_vo_500.rds",sep = ""))
vo_spm <- vo_spm[,grep(pattern = "May",names(vo_spm))]

uo_spm <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/",
        "Current_uo_500.rds",sep = ""))
uo_spm <- uo_spm[,grep(pattern = "May",names(uo_spm))]


calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

data_abun <- list(
  "2021" = data_abun_2021,
  "2022" = data_abun_2022,
  "2023" = data_abun_2023
)

for (annee in c(2021,2022,2023)){
  #select year for each covariates
  data_temp <- sea_floor_temp_spm[,grep(pattern = as.character(annee),
                                   names(sea_floor_temp_spm))]
  names(data_temp) <- c("temp","geometry")
  data_chla <- chla_spm[,grep(pattern = as.character(annee),
                                        names(chla_spm))]
  names(data_chla) <- c("chla","geometry")
  data_uo <- uo_spm[,grep(pattern = as.character(annee),
                                        names(uo_spm))]
  names(data_uo) <- c("uo","geometry")
  data_vo <- vo_spm[,grep(pattern = as.character(annee),
                                        names(vo_spm))]
  names(data_vo) <- c("vo","geometry")
  
  
  data_abun[[as.character(annee)]] <- data_abun[[as.character(annee)]] %>%
    mutate(long = X) %>%
    mutate(lat = Y) %>%
    st_as_sf(.,coords = c("long","lat"),crs = 4467) 
  
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              bathy_spm)
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              data_temp)
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              data_chla)
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              data_uo)
  data_abun[[as.character(annee)]] <- st_join(data_abun[[as.character(annee)]],
                                              data_vo)
}


data_abun_tot_cov <- rbind(data_abun[["2021"]],
                           data_abun[["2022"]],
                           data_abun[["2023"]])

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
              "/process_spatiotemp_data/Data/processed/data_abun_tot_small_cov.rds",
              sep=""))

ggplot(data_abun_tot_cov)+
  geom_sf(aes(color=temp))+
  scale_color_viridis()+
  geom_sf(data=calcul_area, fill = "#11111111",size = 2)+
  facet_wrap(~year, nrow=1)+
  theme(aspect.ratio = 3,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#dadaeeaa"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))


