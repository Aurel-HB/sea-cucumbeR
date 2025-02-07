library(sf)
library(here)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(marmap)
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
                        lat1 = 43, lat2 = 48, resolution = 10)
summary(spm)
plot(spm)

#map ####
ggplot(zee)+
  geom_sf(color = "royalblue", fill = "#cb181dBB" )+
  geom_sf(data=world, fill = "tan")+
  theme(aspect.ratio = 1,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "gray90"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  coord_sf(xlim = c(-58,-54), ylim = c(43,48), expand = FALSE)+
  scale_color_gradient(low = "white", high = "red")+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")+
  annotation_scale(location = "bl", line_width = .5) +
  annotation_north_arrow(location = "tr", height = unit(0.7, "cm"), width = unit(0.7, "cm")) + 
  theme() 
