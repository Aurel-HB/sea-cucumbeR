# this code use the chapter 2 of Spatial Point Process from Baddeley
# the aim is to explore the different haul to see the PPP that define the 
# sea cucumber population.

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)

#############
#load data
#############
data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

#############
#Intensity
#############

#############
#Correlation
#############

#############
#Spacing
#############