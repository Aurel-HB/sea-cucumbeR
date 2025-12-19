#This code extract the bottom temperature from a NetCDF file and process the
# data to keep the period between 2021 and 2025 in the study area : 
# (45.849N, 45.3656S, -56.1402E, -56.26W)

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
ncrast <- raster::stack(paste(here(),
                        "/environment_exploration/Environment_Data/row/",
                        "cmems_mod_glo_phy_my_0.083deg_P1M-m_BottomT.nc",
                        sep=""), varname = "bottomT")
  
print(ncrast)

# Get a data.frame with raster cell values ####
tmp_df <- as.data.frame(ncrast, xy = TRUE) 

name <- c()
for (i in 1:dim(ncrast)[3]){
  months <- c("Jan","Feb","Mar","Apr","May",
              "Jun","Jul","Aug","Sep","Oct","Nov",
              "Dec")
  date <- ncrast@z[["Date/time"]][i]
  month <- as.numeric(substr(date,start = 6,stop = 7))
  year <- substr(date,start = 1,stop = 4)
  name <- c(name,paste(months[month], year,sep = "_"))
}
names(tmp_df) <- c("lon","lat",name)


# transform dataframe in sf
sea_floor_tmp <- st_as_sf(tmp_df, crs = "EPSG:4326", coords = c("lon","lat"))


#map the tmp ####
ggplot(sea_floor_tmp['May_2025'])+
  geom_sf(aes(color=May_2025))+
  scale_color_viridis()+
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
  coord_sf(xlim = c(-56.5,-56), ylim = c(45,46), expand = FALSE)+
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")

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
ncrast_SPM <- projectRaster(ncrast, crs = "epsg:4467", res = 500)
tmp_df_SPM <- as.data.frame(ncrast_SPM, xy = TRUE) 
names(tmp_df_SPM) <- c("lon","lat",name)

plot(ncrast_SPM)

# transform dataframe in sf
sea_floor_tmp_SPM <- st_as_sf(tmp_df_SPM, crs = "EPSG:4467",
                              coords = c("lon","lat"))

#map the tmp #
ggplot(sea_floor_tmp_SPM['May_2025'])+
  geom_sf(aes(color=May_2025))+
  scale_color_viridis()+
  geom_sf(data=calcul_area, fill = "#11111100")+
  geom_point(data=grid, aes(x=long*1000,y=lat*1000), color = "red")+
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
  xlab("")+ylab("")+
  labs(title = "Banc de Saint Pierre et Zone d'étude du projet")

# save the raster to open it as a stars object and convert in sf
writeRaster(ncrast_SPM,
            paste(here(),"/environment_exploration/Environment_Data/processed/",
                  "BottomT_500",sep = ""),
                             format = "GTiff", overwrite=TRUE)
ncstars=read_stars(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
                         "BottomT_500.tif",sep = ""))

# transform stars object in sf
sea_floor_tmp_SPM <- st_as_sf(ncstars, crs = "EPSG:4467",
                              coords = c("lon","lat"))
names(sea_floor_tmp_SPM) <- c(name,"geometry")
saveRDS(sea_floor_tmp_SPM,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
                                 "BottomT_500.rds",sep = ""))

# finally join the prediction grid with the sea_floor_tmp_SPM sf
grid_bottom_tmp <- 
  st_join(grid, sea_floor_tmp_SPM) # to get intersection of points and poly

ggplot(grid_bottom_tmp)+
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

saveRDS(grid_bottom_tmp,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "grid_bottom_tmp.rds",sep = ""))
