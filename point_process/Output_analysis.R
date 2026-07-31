# This code analyse the output from the script Point_Process_Modelling

### ### ### ### ###
#Load packages ####
### ### ### ### ###
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(FactoMineR)

### ### ### ###
#Load data ####
### ### ### ###
list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

scoring_tot <- readRDS(paste(
  here(),"/point_process/Output/scoring_tot.rds",sep=""))

residual_tot <- readRDS(paste(
  here(),"/point_process/Output/residual_tot.rds",sep=""))

smoothed_residuals_raw <- readRDS(
  paste(here(),"/point_process/Output/smoothed_residuals.rds",sep=""))

process_list <- readRDS(paste(
  here(),"/point_process/Output/process_list.rds",sep=""))

### ### ### ###
#Analysis  ####
### ### ### ###

## Prepare data ####
metrics <- scoring_tot %>% 
  mutate(cor.coef = as.numeric(residual_tot$cor.coef)) %>%
  mutate(lm.coef = as.numeric(residual_tot$lm.coef))

## Select model ####
data_select <- data.frame(
  STN = unique(metrics$STN),
  lambda_score = NA,
  K_score = NA,
  cor.coef = NA,
  lm.coef = NA
)
for (stn in unique(metrics$STN)){
  data <- metrics %>% filter(STN == stn) %>%
    filter(!is.na(lambda_score))
  
  data_select$lambda_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$lambda_score), data$lambda_score)]
  
  data_select$K_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$K_score), data$K_score)]
  
  data_select$cor.coef[grep(stn, data_select$STN)] <- 
    data$model[grep(min(abs(1-data$cor.coef)), abs(1-data$cor.coef))]
  
  data_select$lm.coef[grep(stn, data_select$STN)] <- 
    data$model[grep(min(abs(1-data$lm.coef)), abs(1-data$lm.coef))]
  
}






