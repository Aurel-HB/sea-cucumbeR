# This code create a first approach to variography with the HOLOTVSPM2025 data

library(gstlearn)
library(ggplot2)
library(here)
library(dplyr)

## Load data ####
data_abun <- readRDS(
  paste(here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
        sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))