# This code convert and extract the position of haul's start and shoot's end
# with the corresponding time from the roadmap of the survey

# load packages ####
library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(XML)

# load function ####
source(paste(here(), "/via3_data_exploration/fct_degree_decimal.R", sep=""))
source(paste(here(),
             "/via3_data_exploration/fct_WG84_WGSPM.R",
             sep=""))

# load data ####
feuille_route <- read.csv(paste(here(),
                          "/via3_data_exploration/Data/feuille_route_2025.csv",
                                sep=""), sep = ";", header = TRUE,
                          fileEncoding='latin1', check.names=F)

feuille_route <- feuille_route %>% filter(`Validité du trait`== "oui")

#create 2 dataframe with the good coordinates and the important information ####
feuille_route_essential_start <- feuille_route[, c("Date",
                                                   "N° Station",
                                                   "Heure début Virage",
                                                   "Latitude début Virage (N)",
                                                   "Longitude début Virage (W)")]

feuille_route_essential_stop <- feuille_route[, c("Date",
                                                  "N° Station",
                                                  "Heure fin Filage",
                                                  "Latitude fin Filage (N)",
                                                  "Longitude fin Filage (W)")]

names(feuille_route_essential_start) <- c("Date", "Station", "heure",
                                          "latitude", "longitude") 
names(feuille_route_essential_stop) <- c("Date", "Station", "heure",
                                         "latitude", "longitude")
#transform the coordinates to use in sf
for (i in 1:length(feuille_route_essential_start$Date)){
  feuille_route_essential_start$latitude[i] <-
    (deg_dec(feuille_route_essential_start$latitude[i]))
  feuille_route_essential_start$longitude[i] <-
    (-deg_dec(feuille_route_essential_start$longitude[i]))
  feuille_route_essential_stop$latitude[i] <-
    (deg_dec(feuille_route_essential_stop$latitude[i]))
  feuille_route_essential_stop$longitude[i] <-
    (-deg_dec(feuille_route_essential_stop$longitude[i]))
}

feuille_route_essential_start$latitude <- as.numeric(
  feuille_route_essential_start$latitude)
feuille_route_essential_start$longitude <- as.numeric(
  feuille_route_essential_start$longitude)
feuille_route_essential_stop$latitude <- as.numeric(
  feuille_route_essential_stop$latitude)
feuille_route_essential_stop$longitude <- as.numeric(
  feuille_route_essential_stop$longitude)

# transform coordinates to use the SPM projection

coord <- transfo_WG84_WGSPM(feuille_route_essential_start[,c(2,5,4)])
feuille_route_essential_start$longitude <- coord$x
feuille_route_essential_start$latitude <- coord$y
coord <- transfo_WG84_WGSPM(feuille_route_essential_stop[,c(2,5,4)])
feuille_route_essential_stop$longitude <- coord$x
feuille_route_essential_stop$latitude <- coord$y
rm(coord)

# compile in one dataframe ####
names(feuille_route_essential_start)[3:5] <- c("time.haul_start",
                                               "Y.haul_start","X.haul_start")
names(feuille_route_essential_stop)[3:5] <- c("time.shoot_end",
                                               "Y.shoot_end","X.shoot_end")

roadmap <- merge(feuille_route_essential_start, feuille_route_essential_stop,
                 by = intersect(names(feuille_route_essential_start),
                                      names(feuille_route_essential_stop))
                 )
roadmap <- roadmap[order(as.numeric(roadmap$Station)),]

st_write(roadmap,
         paste(here(),"/SIG/SIG_data/roadmap_2025.csv", sep=""),
         append=FALSE)

# check on a map ####
maps <- ggplot(roadmap)+
  geom_point(aes(x = X.haul_start, y= Y.haul_start),color = "yellow")+
  geom_point(aes(x = X.shoot_end, y= Y.shoot_end),
             colour = "red")+
  scale_x_continuous(breaks = seq(-56.40, -56.15, by = 0.125)) +
  annotate("text", x = 544500, y = 5074897, colour = "red",
           size = 4, label = "ZEE")+
  annotate("text", x = 568341, y = 5021189, colour = "darkblue",
           size = 2.2, label = "IFREMER")+
  labs(x = "", y = "", colour = "", breaks = c("1", "3"))+#, caption = "IFREMER")+
  theme(title = element_text(color = "black",face = "bold"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        legend.background = element_rect(fill = "white"),
        legend.key = element_rect(fill = "white", color = NA),
        legend.position.inside = c(0.3, 0.15),
        legend.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        #plot.caption = element_text(hjust = 1, vjust = 1,
        #                            size = 8, color = "darkblue"),
        plot.margin = margin(t = 10,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0), # Left margin
        legend.key.size = unit(.3, "cm"))+
  annotation_scale(location = "bl", line_width = .3) +
  annotation_north_arrow(location = "tr", height = unit(0.5, "cm"),
                         width = unit(0.5, "cm")) +
  theme()

maps
