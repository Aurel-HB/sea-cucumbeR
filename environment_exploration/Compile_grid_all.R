# This script compile all the covariate on the prediction grid for 2021 to 2025

# load package
library(here)
library(sf)
library(dplyr)
library(ggplot2)

# load grid
grid_bathy <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/grid_bathy.rds",sep = "")
)
grid_tmp <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/grid_bottom_tmp.rds",sep = "")
)
grid_chl <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/grid_chl.rds",sep = "")
)
grid_uo <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/grid_uo.rds",sep = "")
)
grid_vo <- readRDS(
  paste(here(),"/environment_exploration/Environment_Data/processed/grid_vo.rds",sep = "")
)

# The survey is only in may so we keep the may data
grid_tmp <- grid_tmp[,grep(pattern = "May",names(grid_tmp))]
grid_chl <- grid_chl[,grep(pattern = "May",names(grid_chl))]
grid_uo <- grid_uo[,grep(pattern = "May",names(grid_uo))]
grid_vo <- grid_vo[,grep(pattern = "May",names(grid_vo))]

# The grid have to be replicated for each year with the environmental covariate values
grid_tot <- data.frame()
for (i in 1:5){
  data <- st_join(grid_tmp[,i],grid_bathy)
  data <- st_join(data,grid_chl[,i])
  data <- st_join(data,grid_uo[,i])
  data <- st_join(data,grid_vo[,i])
  data <- data %>% mutate(year = 2020+i)
  names(data) <- c("bottomT","long","lat","bathy","chla","uo","vo","geometry","year")
  grid_tot <- rbind(grid_tot,data[,c(1,4,5,6,7,2,3,9,8)])
} # the lat and long are in km for modelling

saveRDS(grid_tot,
        paste(here(),
              "/environment_exploration/Environment_Data/processed/grid_tot.rds",sep = "")
        )
