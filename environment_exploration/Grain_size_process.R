#This code extract the grain size and process the data to keep the study area

#load the packages
library(here)
library(marmap)
library(lattice)
library(RColorBrewer)
library(sf)
library(viridis)
library(ggplot2)
library(terra)
library(gstat)
library(stars)
library(dplyr)
library(readxl)

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


#import data
Atlantic_grain_size <- 
  read_excel(paste(here(),
    "/environment_exploration/Environment_Data/",
    "row/M.Li_atl_grainsize_gridded.xlsx", sep=""))

Canadian_grain_size <- 
  read_excel(paste(here(),
                   "/environment_exploration/Environment_Data/",
                   "row/M.Li_Canada_080603_Atlantic_grain size.xlsx", sep=""))

ggplot(Canadian_grain_size, aes(x=STN_LON,y=STN_LAT))+
  geom_point(aes(color=SAND))

Canadian_grain_size <- Canadian_grain_size %>%
  mutate(lat=STN_LAT) %>%
  mutate(long=STN_LON) %>%
  st_as_sf(., coords = c("STN_LON","STN_LAT"), crs = st_crs(4326))

ggplot()+
  geom_point(data = Canadian_grain_size %>%
               filter(45<lat) %>%
               filter(lat<46) %>%
               filter(-57<long) %>%
               filter(long<(-56)),
             aes(x=long,y=lat,color=SAND))+
  scale_color_viridis()


long_lim <- c(lon1 = -58, lon2 = -54)
lat_lim <- c(lat1 = 43, lat2 = 48)

area_grain_size <- Atlantic_grain_size %>%
  filter(lat_lim[1]<lat) %>%
  filter(lat<lat_lim[2]) %>%
  filter(long_lim[1]<long) %>%
  filter(long<long_lim[2])

ggplot()+
  geom_point(data= area_grain_size, aes(x=long, y=lat, 
                                        color=GrainSize_mm))#color=SAND))+
  #geom_sf(data = calcul_area)

area_grain_size <- area_grain_size %>%
  mutate(x=long)%>%
  mutate(y=lat)%>%
  st_as_sf(., coords = c("x","y"), crs = st_crs(4326))%>%
  st_transform(., crs = st_crs(4467))

area_grain_size$X <- st_coordinates(area_grain_size)[,1]
area_grain_size$Y <- st_coordinates(area_grain_size)[,2]

ggplot()+
  geom_sf(data = calcul_area)+
  geom_point(data = area_grain_size %>%
               filter(45<lat) %>%
               filter(lat<46) %>%
               filter(-57<long) %>%
               filter(long<(-56)),
             aes(x=X,y=Y,color=log(GrainSize_mm)))+
  scale_color_viridis()

ras_dom<-raster(xmn=long_lim[1], xmx=long_lim[2], ymn=lat_lim[1], ymx=lat_lim[2],
                crs="+proj=longlat +datum=WGS84 +no_defs ",
                resolution=c(0.1,0.1), vals=NA)
layer <- as.data.frame(area_grain_size)[,c("GrainSize_mm","long","lat")]
coordinates(layer) <- ~long+lat
result <- rasterize(layer,ras_dom,"GrainSize_mm", update=TRUE)
plot(result)

# Kriging the grain size data to add the infromation into the prediciton grid ####
# --- Step 1: Create Sample Data (Replacing your actual dataframe) ---
df <- as.data.frame(area_grain_size)[,c("GrainSize_mm","long","lat")]

# --- Step 2: Convert Dataframe to a Spatial Object ---
# We convert to 'sf' (simple features) which is the modern standard for vector data
df_sf <- st_as_sf(df, coords = c("long", "lat"), crs = 4326)

# --- Step 3: Create a Grid (The template for your raster) ---
# Define the extent and resolution of your output raster
# Resolution is in degrees here because we are using WGS84 (lat/long)
grid_sf <- rast(ext(long_lim[1], long_lim[2], lat_lim[1], lat_lim[2]),
             res = 0.01, crs = "EPSG:4326")

# Convert the empty grid to a spatial object for interpolation
grid_sf <- st_as_sf(as.points(grid))

# --- Step 4: Interpolation (Method A: IDW - Fast & Simple) ---
# IDW assumes closer points have more influence than distant points.
idw_model <- gstat(formula = GrainSize_mm ~ 1, locations = df_sf)
idw_result <- predict(idw_model, grid_sf)

# Convert the result back to a terra SpatRaster
# Create the empty raster template
raster_idw <- grid 
# Assign the calculated values to it
values(raster_idw) <- idw_result$var1.pred

# --- Step 5: Interpolation (Method B: Kriging - More Accurate/Statistical) ---
# Kriging is better if your data has a spatial structure (autocorrelation)
# Note: For true Kriging, you should technically model a variogram first.
vgm <- variogram(GrainSize_mm ~ 1, df_sf)
vgm_fit <- fit.variogram(vgm, model = vgm("Sph")) # Fitting a Spherical model
krig_model <- gstat(formula = GrainSize_mm ~ 1, locations = df_sf, model = vgm_fit)
krig_result <- predict(krig_model, grid_sf)

# Create the empty raster template
raster_krig <- grid 
# Assign the calculated values to it
values(raster_krig) <- krig_result$var1.pred

# --- Step 6: Visualization ---
par(mfrow = c(1, 2))
plot(raster_idw, main = "IDW Interpolation")
points(df_sf, col = "red", pch = 20)

plot(raster_krig, main = "Kriging Interpolation")
points(df_sf, col = "red", pch = 20)
par(mfrow = c(1, 1))

# convert the coordinate of the raster to the epsg 4467 of SPM
ncrast_SPM <- project(raster_krig, "EPSG:4467")
ncrast_SPM <- project(raster_idw, "EPSG:4467")

# save the raster to open it as a stars object and convert in sf
writeRaster(ncrast_SPM,
            paste(here(),"/environment_exploration/Environment_Data/processed/",
                  "GrainSize_3PS.tif",sep = ""), overwrite=TRUE)
ncstars=read_stars(paste(
  here(),"/environment_exploration/Environment_Data/processed/",
  "GrainSize_3PS.tif",sep = ""))

# transform stars object in sf
GrainSize_spm <- st_as_sf(ncstars, crs = "EPSG:4467",
                      coords = c("lon","lat"))
names(GrainSize_spm) <- c("GrainSize","geometry")
saveRDS(GrainSize_spm,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "GrainSize_3PS.rds",sep = ""))

# finally join the prediction grid with the sea_floor_tmp_SPM sf
grid_GrainSize <- 
  st_join(grid, GrainSize_spm) # to get intersection of points and poly

ggplot()+
  geom_sf(data = grid_GrainSize, aes(color = GrainSize))+
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

saveRDS(grid_GrainSize,
        paste(here(),"/environment_exploration/Environment_Data/processed/",
              "grid_GrainSize.rds",sep = ""))





# test ####

area_grain_size <- st_join(
  calcul_area, area_grain_size)

ggplot()+
  geom_sf(data = area_grain_size)+
  geom_point(data = area_grain_size, aes(x=X,y=Y))
