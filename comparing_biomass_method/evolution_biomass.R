# This code compare the biomass between sdmTMB and bootstrap 
# It show also the evolution between the different survey

library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)

# load data
evolution_biomass <- readRDS(
  paste(here(),"/comparing_biomass_method/Data/evolution_biomass.rds", sep=""))
  
# show evolution
ggplot(evolution_biomass[c(1,2,3,5),],aes(x = date, y = biomass))+
  geom_point(colour = "grey30")+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2)+
  geom_point(aes(x = date, y = TAC))+
  scale_y_log10()+
  labs(x = "", y = "Biomass")

test <- evolution_biomass[c(1,2,3,5),]
test <- rbind(test[,1:4],
              data.frame(date=test$date,biomass=test$TAC,lwr=NA,upr=NA))
test$indic <- c(rep("Biomass from survey",4),rep("TAC",4))
test <- rbind(test, data.frame(date=2024,biomass=2260,lwr=NA,upr=NA,indic="TAC"))

ggplot(test,aes(x = date, y = biomass))+
  geom_point(colour = "grey30")+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2)+
  facet_wrap(~indic, nrow = 2,scales="free_y")+
  labs(x = "", y = "Biomass (tonnes)")+
  theme_bw()

test$origine <- "black"
test <- rbind(test, data.frame(date=2025.1,biomass=evolution_biomass$biomass[6],
                               lwr=evolution_biomass$lwr[6],
                               upr=evolution_biomass$upr[6],
                               indic="Biomass from survey",
                               origine="lightblue"))

ggplot(test,aes(x = as.numeric(date), y = biomass,colour = origine))+
  geom_point(show.legend = F)+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2,show.legend = F)+
  facet_wrap(~indic, nrow = 2,scales="free_y")+
  labs(x = "", y = "Biomass (tonnes)")+
  theme_bw()
