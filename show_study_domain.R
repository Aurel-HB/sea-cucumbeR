library(sf)
library(here)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(marmap)
library(tidyverse)
library(RColorBrewer)
#library(mapsf)

world <- ne_download(scale = 10, type = 'countries', returnclass = "sf")

#import the grid and the area
shp_grid <- st_read(dsn=paste(here(),"/SIG_Data", sep=""), layer = "grille_zee_spm")
#zee <- st_transform(shp, crs = "ESRI:102002")
grid_zee <- st_transform(shp_grid, crs = "EPSG:4326")
grid_zee <- grid_zee$geometry

# from www.marineregions.org
shp <- st_read(dsn=paste(here(),"/SIG_Data", sep=""), layer = "eez")
zee <- st_transform(shp, crs = "EPSG:4326")
zee <- zee$geometry

#bathy ####
#file to heavy
#shp_bathy <- st_read(dsn=paste(here(),"/SIG_Data", sep=""), layer = "bathy_10m")
#bathy <- st_transform(shp_bathy, crs = "EPSG:4326")
#bathy <- bathy$geometry
#ggplot()+geom_sf(data=bathy, fill = "navyblue")

spm <- getNOAA.bathy(lon1 = -58, lon2 = -54,
                        lat1 = 43, lat2 = 48, resolution =1)
summary(spm)
plot(spm)
#test autoplot
autoplot.bathy(spm, geom=c("tile","contour")) +
  scale_fill_gradient2(low="dodgerblue4", mid="gainsboro", high="darkgreen") +
  #geom_point(data = ctd, aes(x = Longitude, y = Latitude),
  #           colour = 'black', size = 3, alpha = 1, shape = 15) +
  labs(y = "Latitude", x = "Longitude", fill = "Elevation") +
  coord_cartesian(expand = 0)+
  ggtitle("A marmap map with ggplot2") 

Bathy <- as.matrix(spm)
class(Bathy) <- "matrix"

Bathy_offshore <- Bathy %>%
  as.data.frame() %>%
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric)) %>%
  dplyr::filter(value <= -500) %>%
  mutate(value = value *-1)

Bathy_coast <- Bathy %>%
  as.data.frame() %>%
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric)) %>% 
  dplyr::filter(between(value,-500,0)) %>%
  mutate(value = value *-1)


#lat <- round(unique(Bathy$lat), digits = 3)
ggplot()+
  geom_contour(data = Bathy_coast, aes(x = lon, y = lat, z = value), bins = 4, colour = "black") +
  geom_contour(data = Bathy_offshore, aes(x = lon, y = lat, z = value), bins = 6, colour = "black") +
  coord_map()

#map ####

ggplot(zee)+
  geom_sf(color = "royalblue", fill = "#78c679BB" )+
  geom_sf(data=world, fill = "tan")+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0d1e6"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  coord_sf(xlim = c(-58,-54), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tr", height = unit(0.7, "cm"), width = unit(0.7, "cm")) + 
  theme()

ggplot(zee)+
  geom_contour_filled(data = Bathy_coast, aes(x = lon, y = lat, z = value),
               bins = 4, colour = "black", show.legend = FALSE) +
  geom_contour_filled(data = Bathy_offshore, aes(x = lon, y = lat, z = value),
               bins = 6, colour = "black", show.legend = FALSE) +
  discrete_scale('fill', 'myscale', colorRampPalette(c("#deebf7", "#63a4d8")))+
  geom_sf(color = "royalblue", fill = "#78c679BB" )+
  geom_sf(data=world, fill = "tan")+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "#d0d1e6"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  coord_sf(xlim = c(-58,-54), ylim = c(43,48), expand = FALSE)+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tr", height = unit(0.7, "cm"), width = unit(0.7, "cm")) + 
  theme() 
