# This code create a biomass index by the extrapolation by square then bootstrap

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(boot)


shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))
predata_intensity <- readRDS(
  paste(here(), "/point_process/Output/holo_simu_intensity.rds", sep=""))
data_abun <- readRDS(
  paste(here(),"/point_process/Output/holo_simu_abun.rds", sep=""))

# param for real data
vitesseNavire <- 1.4 # (in mile/h) can be changed by a vector if available
intervalTemps <- 10 # (in min)  can be changed by a vector if available
champVision <- 1.5/1000 # (in kilometer)
surfCompte <- as.numeric(intervalTemps*champVision*vitesseNavire*1852/60)
mean_weight <- 0.4 # (in kg)
surfsquare <- 1.852*1.852 # in kilometer square

# param for simulated data 
surfCompte <- 0.5*champVision # in kilometer square

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
