# This code test an apporach of modelisation on point process to characterize it

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

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

HOLOTVSPM2025 <- readRDS(paste(here(),
              "/comparing_biomass_method/Data/03_HOLOTVSPM2025_summary.rds",
              sep=""))

############
# Start with one PPP
############

PPP <- list_PPP[[1]]
stn <- names(list_PPP)[1]
d <- 

# spatialization of the PPP




