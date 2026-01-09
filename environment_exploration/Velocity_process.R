#This code extract the Eastward & Northward sea water velocity in sea water from 
# a NetCDF file and process the data to keep the period between 2021 and 2025 
# at shallower bathy in the study area: (45.849N, 45.3656S, -56.1402E, -56.26W)
#    (-37.92902m)


# load the package
#library(ncdf4)
library(here)
#library(CFtime)
library(lattice)
library(RColorBrewer)
library(sf)
library(viridis)
library(ggplot2)
library(raster)
library(stars)
library(dplyr)

#load the study area
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

#open a netCDF file 
ncrast_east <- raster::stack(paste(here(),
                              "/environment_exploration/Environment_Data/row/",
                        "cmems_mod_glo_phy_my_0.083deg_velocity_34.43m.nc",
                              sep=""), varname = "uo")

ncrast_north <- raster::stack(paste(here(),
                              "/environment_exploration/Environment_Data/row/",
                         "cmems_mod_glo_phy_my_0.083deg_velocity_34.43m.nc",
                                   sep=""), varname = "vo")

print(ncrast_east)
print(ncrast_north)


# Get a data.frame with raster cell values ####
#Eastward sea water velocity
uo_df <- as.data.frame(ncrast_east, xy = TRUE) 

name <- c()
for (i in 1:dim(ncrast_east)[3]){
  months <- c("Jan","Feb","Mar","Apr","May",
              "Jun","Jul","Aug","Sep","Oct","Nov",
              "Dec")
  date <- ncrast_east@z[["Date/time"]][i]
  month <- as.numeric(substr(date,start = 6,stop = 7))
  year <- substr(date,start = 1,stop = 4)
  name <- c(name,paste(months[month], year,sep = "_"))
}
names(uo_df) <- c("lon","lat",name)

#Northward sea water velocity
vo_df <- as.data.frame(ncrast_north, xy = TRUE) 

name <- c()
for (i in 1:dim(ncrast_north)[3]){
  months <- c("Jan","Feb","Mar","Apr","May",
              "Jun","Jul","Aug","Sep","Oct","Nov",
              "Dec")
  date <- ncrast_north@z[["Date/time"]][i]
  month <- as.numeric(substr(date,start = 6,stop = 7))
  year <- substr(date,start = 1,stop = 4)
  name <- c(name,paste(months[month], year,sep = "_"))
}
names(vo_df) <- c("lon","lat",name)


# transform dataframe in sf
velocity_uo <- st_as_sf(uo_df, crs = "EPSG:4326", coords = c("lon","lat"))
velocity_vo <- st_as_sf(vo_df, crs = "EPSG:4326", coords = c("lon","lat"))

#fuse the raster with the grid of prediction used ####

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

# convert the coordinate of the raster to the epsg 4467 of SPM
ncrast_east_SPM <- projectRaster(ncrast_east, crs = "epsg:4467", res = 500)
uo_df_SPM <- as.data.frame(ncrast_east_SPM, xy = TRUE) 
names(uo_df_SPM) <- c("lon","lat",name)

plot(ncrast_east_SPM)

ncrast_north_SPM <- projectRaster(ncrast_north, crs = "epsg:4467", res = 500)
vo_df_SPM <- as.data.frame(ncrast_north_SPM, xy = TRUE) 
names(vo_df_SPM) <- c("lon","lat",name)

plot(ncrast_north_SPM)

# save the rasters to open it as a stars object and convert in sf
writeRaster(ncrast_east_SPM,
            paste(here(),"/environment_exploration/Environment_Data/processed/",
                  "Current_uo_500",sep = ""),
            format = "GTiff", overwrite=TRUE)
ncstars_east=read_stars(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
  "Current_uo_500.tif",sep = ""))

writeRaster(ncrast_north_SPM,
            paste(here(),"/environment_exploration/Environment_Data/processed/",
                  "Current_vo_500",sep = ""),
            format = "GTiff", overwrite=TRUE)
ncstars_north=read_stars(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
  "Current_vo_500.tif",sep = ""))

# transform stars object in sf
uo_SPM <- st_as_sf(ncstars_east, crs = "EPSG:4467",
                         coords = c("lon","lat"))
names(uo_SPM) <- c(name,"geometry")
saveRDS(uo_SPM,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "Current_uo_500.rds",sep = ""))

vo_SPM <- st_as_sf(ncstars_north, crs = "EPSG:4467",
                         coords = c("lon","lat"))
names(vo_SPM) <- c(name,"geometry")
saveRDS(vo_SPM,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "Current_vo_500.rds",sep = ""))

# finally join the prediction grid with the mass_chl_SPM sf
grid_uo <- 
  st_join(grid, uo_SPM) # to get intersection of points and poly

grid_vo <- 
  st_join(grid, vo_SPM) # to get intersection of points and poly

ggplot(grid_uo)+
  geom_sf(aes(color=May_2025))+
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

saveRDS(grid_uo,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "grid_uo.rds",sep = ""))

ggplot(grid_vo)+
  geom_sf(aes(color=May_2025))+
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

saveRDS(grid_vo,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "grid_vo.rds",sep = ""))
