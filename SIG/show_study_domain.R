library(sf)
library(here)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(marmap)
library(tidyverse)
library(RColorBrewer)
library(viridis)
#library(mapsf)

world <- ne_download(scale = 10, type = 'countries', returnclass = "sf")

#import the grid and the area
shp_grid <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""), layer = "grille_zee_spm")
#zee <- st_transform(shp, crs = "ESRI:102002")
grid_zee <- st_transform(shp_grid, crs = "EPSG:4326")
grid_zee <- grid_zee$geometry

# from www.marineregions.org
shp <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""), layer = "eez")
zee <- st_transform(shp, crs = "EPSG:4326")
zee <- zee$geometry

#bathy ####
#file to heavy
#shp_bathy <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""), layer = "bathy_10m")
#bathy <- st_transform(shp_bathy, crs = "EPSG:4326")
#bathy <- bathy$geometry
#ggplot()+geom_sf(data=bathy, fill = "navyblue")

spm <- getNOAA.bathy(lon1 = -58, lon2 = -54,
                        lat1 = 43, lat2 = 48, resolution =1)
summary(spm)
plot(spm)
#test autoplot
autoplot.bathy(spm, geom=c("tile","contour")) +
  scale_fill_gradient2(low="navy", mid="lightblue", high="darkgreen") +
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

Bathy_zee <- Bathy %>%
  as.data.frame() %>%
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric))
Bathy_zee <- st_as_sf(Bathy_zee, crs = "EPSG:4326", coords = c("lon","lat"))

Bathy_tot <- Bathy %>%
  as.data.frame() %>%
  rownames_to_column(var = "lon") %>%
  gather(lat, value, -1) %>%
  mutate_all(funs(as.numeric))

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
  geom_sf(data = Bathy_zee, aes(color=value))+
  scale_color_viridis(option = "mako")+
  geom_sf(color = "royalblue", fill = "#fcfbfdAA" )+
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
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) + 
  theme() 

ggplot(zee)+
  geom_contour_filled(data = Bathy_coast, aes(x = lon, y = lat, z = value),
               bins = 4, colour = "black", show.legend = FALSE) +
  geom_contour_filled(data = Bathy_offshore, aes(x = lon, y = lat, z = value),
               bins = 6, colour = "black", show.legend = FALSE) +
  discrete_scale('fill', 'myscale', colorRampPalette(c("#deebf7", "#306094")))+#63a4d8
  geom_sf(color = "royalblue", fill = "#e41a1cAA" )+
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
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) + 
  theme() 

#sea floor tmp
sea_floor_tmp <- readRDS(paste(here(), "/SIG/SIG_Data/sea_floor_tmp.rds",
                               sep=""))

ggplot(zee)+
  geom_sf(data = sea_floor_tmp['tmpJan'], aes(color=tmpJan), size = 6)+
  scale_color_viridis()+
  geom_contour(data = Bathy_tot, aes(x = lon, y = lat, z = value),
                      bins = 25, colour = "black") +
  geom_sf(color = "royalblue", fill = "#fcfbfdAA" )+
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
  annotation_north_arrow(location = "tr", height = unit(1, "cm"), width = unit(1, "cm")) + 
  theme() 

# We want to keep the data of the domain of study
location <- st_intersects(zee,sea_floor_tmp)
obs <- sea_floor_tmp[location[[1]],]

ggplot()+
  geom_sf(data = zee, color = "royalblue", fill = "#fcfbfdAA" )+
  geom_sf(data = obs['tmpJan'], aes(color=tmpJan))+
  scale_color_viridis()+
  theme_bw()

location <- st_intersects(zee,Bathy_zee)
obs <- Bathy_zee[location[[1]],]

ggplot()+
  geom_sf(data = obs, aes(color=value))+
  geom_sf(data = zee, color = "royalblue", fill = "#fcfbfdAA" )+
  theme_bw()


# create the polygon of the tuyau study area
poly_tuyau <- data.frame(WGSPMX = c(557407,566667), WGSPMY = c(5077675,5023967))
poly_tuyau <- poly_tuyau %>% 
  st_as_sf(coords = c("WGSPMX", "WGSPMY"), 
           crs = 4467) %>% 
  st_bbox() %>% 
  st_as_sfc() %>%
  st_sf()
#st_area(poly_tuyau) # calculate the area of the polygon
#saveRDS(poly_tuyau,
# "C:/Users/ahebertb/Documents/sea-cucumbeR/SIG/SIG_Data/study_calcul_area.rds")

# line that separate the zone 1 and the zone 2 of sea cucumber fishery 
# legislation at Saint Pierre and Miquelon

Point_A <- c("46°31ʹ19","56°47ʹ59")
Point_B <- c("46°31ʹ15","56°28ʹ54")
Point_C <- c("46°42ʹ03","56°28ʹ48")
Point_D <- c("46°41ʹ19","55°55ʹ25")

deg_dec <- function(coord){
  # coords is the shape xx°xx°xxx 
  dec <- as.numeric(substr(coord, start = 1, stop = 2))
  min <- as.numeric(substr(coord, start = 4, stop = 5))/60
  sec <- as.numeric(substr(coord, start = 7, stop = nchar(coord)))/3600
  return(dec+min+sec)
}

delimitation <- data.frame(
  Id = c("Point_A","Point_B","Point_C","Point_D"),
  long = c(Point_A[2],Point_D[2],Point_C[2],Point_D[2]),
  lat = c(Point_A[1],Point_D[1],Point_C[1],Point_D[1])
)

for (i in 1:nrow(delimitation)){
  for (y in 2:ncol(delimitation)){
    delimitation[i,y] <- deg_dec(delimitation[i,y])
  }
}
delimitation$long <- (-1)*as.numeric(delimitation$long)
delimitation$lat <- as.numeric(delimitation$lat)

sf_delim <- delimitation %>% 
  st_as_sf(coords = c("long", "lat"), 
           crs = 4326)



ggplot(zee)+
  geom_sf(color = "royalblue", fill = "#78c679BB" )+
  geom_sf(data=sf_delim)+
  geom_line(data = delimitation, aes(x = long, y = lat))+
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
  theme()
