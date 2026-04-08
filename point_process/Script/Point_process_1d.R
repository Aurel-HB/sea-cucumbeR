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

###################################
# Start with one PPP with spatstat
###################################
stn=96# stn=149 stn=179 stn=143
PPP <- PPP<-list_PPP[[as.character(stn)]]
PPP[["window"]][["units"]] <- list("metre","metres")
line_data <- data.frame(x = 0,
                        y = data_position[grep(
                          as.character(stn),data_position$station),]$Y)
lPPP <- ppp(line_data$x,line_data$y, window = owin(xrange = c(0,0),
                                                   yrange = PPP$window$yrange))
lPPP[["window"]][["units"]] <- list("metre","metres")

# explanatory stat test ####
A <- hopskel(lPPP)# random:A=1, cluster:A<1, regularity:A>1

y_sorted <- sort(line_data$y)
# Compute Sn(y) = cumulative proportion
n <- length(y_sorted)
Sn <- cumsum(rep(1/n, n)) # Sn(y)** is the cumulative proportion equivalent to 
#the empirical distribution function in the KS test.
# KS test against a uniform distribution (conditional on n)
ks_result <- ks.test(y_sorted, "punif", min(y_sorted), max(y_sorted))
# empirical CDF
plot(ecdf(y_sorted), main = "Empirical vs Uniform CDF")
# theoretical CDF
curve(punif(x, min = min(y_sorted), 
            max = max(y_sorted)), add = TRUE, col = "red")

#####

endpoints <- ppp(x=c(0,0), y=PPP$window$yrange, 
                 window = owin(c(-0.1,0.1), PPP$window$yrange))
L <- linnet(endpoints, edges = matrix(c(1,2),ncol = 2))
#X <- rpoislpp(lambda = 0.1, L = L)
#plot(X, main = "")
#axis(2)
X <- lpp(lPPP,L)
plot(X, main = "")
axis(2)

#Poisson model ####
fit0 <- lppm(X~1)
fit1 <- lppm(X~y)
#Likelihood Ratio Test
anova(fit0, fit1, test="LR")
summary(fit1)

