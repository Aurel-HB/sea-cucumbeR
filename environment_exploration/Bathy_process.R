#This code extract the bathy and process the data to keep the study area

#load the packages
library(here)
library(marmap)
library(lattice)
library(RColorBrewer)
library(sf)
library(viridis)
library(ggplot2)
library(raster)
library(stars)
library(dplyr)

#Study area ####
#load the study area
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

# use the study area polygon to generate an extrapolation grid
coordinate <- as.data.frame(calcul_area[[1]][[1]][[1]])/1000
names(coordinate) <- c("X","Y")
res <- 0.5
grid <- expand.grid(
  long = seq(min(coordinate$X), max(coordinate$X), by = res),
  lat = seq(min(coordinate$Y), max(coordinate$Y), by = res)
)
plot(grid)

grid <- grid %>%
  mutate(X = long*1000)%>%
  mutate(Y= lat*1000)%>%
  st_as_sf(crs = "EPSG:4467",
           coords = c("X","Y"))

#Bathy ####
#extract bathy in the area
spm <- getNOAA.bathy(lon1 = -58, lon2 = -54,
                     lat1 = 43, lat2 = 48, resolution =1)
summary(spm)
plot(spm)

#Convert bathymetric data to a raster layer
ncrast <- marmap::as.raster(spm)
plot(ncrast)

# convert the coordinate of the raster to the epsg 4467 of SPM
ncrast_SPM <- projectRaster(ncrast, crs = "epsg:4467")

# save the raster to open it as a stars object and convert in sf
writeRaster(ncrast_SPM,
            paste(here(),"/environment_exploration/Environment_Data/processed/",
                  "Bathy_3PS",sep = ""),
            format = "GTiff", overwrite=TRUE)
ncstars=read_stars(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
  "Bathy_3PS.tif",sep = ""))

# transform stars object in sf
bathy_spm <- st_as_sf(ncstars, crs = "EPSG:4467",
                              coords = c("lon","lat"))
names(bathy_spm) <- c("bathy","geometry")
saveRDS(bathy_spm,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "Bathy_3PS.rds",sep = ""))

# finally join the prediction grid with the sea_floor_tmp_SPM sf
grid_bathy <- 
  st_join(grid, bathy_spm) # to get intersection of points and poly

ggplot()+
  geom_sf(data = grid_bathy, aes(color = bathy))+
  scale_color_viridis()+
  geom_sf(data=calcul_area, fill = "#11111100")+
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
                                        colour = "white"))

saveRDS(grid_bathy,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "grid_bathy.rds",sep = ""))
