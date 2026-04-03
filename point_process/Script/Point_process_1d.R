# This code test an approach of modelling on 1d point process and 
# line transect analysis

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
#library(spatstat.model)
#library(spastat.linnet)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(fields)
library(gridExtra)

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

data_position_substrat <- readRDS(paste(here(),
                                        "/via3_data_exploration/Data/substrat/data_position_substrat.rds",
                                        sep=""))

data_substrate <- readRDS(paste(here(),
                                "/via3_data_exploration/Data/substrat/data_substrat_2025.rds",
                                sep=""))

########################
# Start with one PPP
########################