# This code create a data sf with the centroid of the tracks from 2025 
# with the polygon's coordinates of each station.

# load packages ####
library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(XML)

# load function ####
source(paste(here(), "/via3_data_exploration/fct_WG84_WGSPM.R", sep=""))

# load data ####
gpx_parsed <- htmlTreeParse(file = paste(here(),
                                         "/SIG/SIG_Data/tracks_2025",
                                         ".gpx",sep=""),
                            useInternalNodes = TRUE)

coords <- xpathSApply(doc = gpx_parsed, path = "//trkpt", fun = xmlAttrs)
time <- xpathSApply(doc = gpx_parsed, path = "//trkpt/time", fun = xmlValue)

sf_sampling_grid <- readRDS(paste(here(),"/SIG/SIG_Data/sf_sampling_grid.rds",
                            sep=""))

# prepare data ####
#tracks
df <- data.frame(
  latitude = as.numeric(coords["lat", ]),
  longitude = as.numeric(coords["lon", ]),
  time = time
)

# conversion of coordinates
coord_tracks <- data.frame(id = 1:length(df$lat), df[,c(2,1)])
coord_tracks <- transfo_WG84_WGSPM(coord_tracks)
names(coord_tracks)[2:3] <- c("WGSPMX", "WGSPMY")
coord_tracks <- data.frame(coord_tracks, time = df$time)

# preparation of the sf table
df_sf <- data.frame(coord_tracks,
                    X.track = coord_tracks$WGSPMX,
                    Y.track = coord_tracks$WGSPMY)
df_sf <- st_as_sf(df_sf, coords = c("WGSPMX","WGSPMY"))
df_sf <- st_set_crs(df_sf, 4467)

df_sf$time <- as.factor(substr(df_sf$time, start = 1, stop = 10))


#grid ##
# filter the square that are in the sampling zone
poly <- readRDS(
  paste(here(),"/SIG/SIG_Data/study_calcul_area.rds",sep="")
)
sf_sampling_grid <- st_filter(sf_sampling_grid, poly)
sf_sampling_grid <- sf_sampling_grid %>% 
  select(Echantillo,Id,color_sampling,geometry)
names(sf_sampling_grid) = c("Echantillo","Station","color_sampling","geometry")

# join and keep only the key data
sf_tot <- st_join(sf_sampling_grid,df_sf)

# process data ####
# for having the station, the type of sampling, the centroïd of each video tracks
# and the geometry of the statistic square associated

sf_temp <- sf_tot %>% select(Station,X.track,Y.track)%>%
  filter(!is.na(X.track))

# select centroïd of the tracks
centroid.tracks <- data.frame(station = unique(sf_temp$Station),
                              X.track = 0,
                              Y.track = 0) %>%
  filter(!is.na(station))

for (i in 1:nrow(centroid.tracks)){
  stn <- centroid.tracks$station[i]
  data <- sf_temp %>% filter(Station == stn)
  centroid.tracks$X.track[i] <- mean(data$X.track)
  centroid.tracks$Y.track[i] <- mean(data$Y.track)
}
rm(i,stn, data, sf_temp)

# prepare the final sf dataframe
sf_final <- sf_sampling_grid %>%
  mutate(X.track = NA) %>%
  mutate(Y.track = NA)

for (i in 1:nrow(sf_final)){
  stn <- sf_final$Station[i]
  if (stn %in% centroid.tracks$station){
    data <- centroid.tracks %>% filter(station == stn)
    sf_final$X.track[i] <- data$X.track[1]
    sf_final$Y.track[i] <- data$Y.track[1]
  }
}
rm(i,stn, data)

# correct the information ####
# the station is a artifact due to some track point of the station 163
sf_final$X.track[92] <- NA
sf_final$Y.track[92] <- NA
# due to a wrong manipulation, we don't have the tracks of the station 144
# but we have the start (564521.7,5059426) and stop (564620.4,5060177) of the sample
sf_final$X.track[74] <- mean(c(564521.7,564620.4))
sf_final$Y.track[74] <- mean(c(5059426,5060177))

# check on a map ####
maps <- ggplot()+
  scale_fill_manual(
    name = paste("Plan d'échantillonnage",2025,sep = " "),
    values = c("#990000", "#fffbfd00", "lightsalmon3", "seagreen3"))+
  guides(fill = guide_legend(override.aes = list(size=4)))+
  geom_sf(data = sf_final, color = "black",
          fill = sf_final$color_sampling)+
  geom_point(data = sf_final, aes(x = X.track, y= Y.track),
                                  colour = "black")+
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

# lay out the final shape of the sf tuyau area ####
sf_final <- sf_final %>%
  select(Station,X.track, Y.track, geometry)%>%
  filter(!is.na(X.track)) %>%
  st_centroid() %>%
  mutate(X.centroid_station = st_coordinates(.)[,1]) %>%
  mutate(Y.centroid_station = st_coordinates(.)[,2])


# and the station of the boite à pétoncle have to add ####
# if wanted but all the analysis will most only on the tuyau area
bap_grid <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""),
                    layer = "grille_boite_petoncle_rectangle")
bap_grid <- bap_grid %>% st_centroid() %>%
  mutate(X.centroid_station = st_coordinates(.)[,1]) %>%
  mutate(Y.centroid_station = st_coordinates(.)[,2]) %>%
  filter(idbis %in% c("43B","50B","57B"))
# we just have the stop and start of the each sample
bap <- data.frame(bap_grid$idbis,
                  c(NA,NA,NA),
                  c(NA,NA,NA),
                  bap_grid$geometry)  %>%
  st_as_sf()
                  
names(bap)[1:3] <- names(sf_final)[1:3]

bap <- bap %>% 
  mutate(X.centroid_station = bap_grid$X.centroid_station) %>%
  mutate(Y.centroid_station = bap_grid$Y.centroid_station)

# 43B : start (516757.7, 5165033) and stop (517131.6,5164992)
bap$X.track[1] <- mean(c(516757.7,517131.6))
bap$Y.track[1] <- mean(c(5165033,5164992))
# 50B : start (518338.1,5165069) and stop (518842.0,5165037)
bap$X.track[2] <- mean(c(518338.1,518842.0))
bap$Y.track[2] <- mean(c(5165069,5165037))
# 57B : start (519837.6,5163713) and stop (520073.7,5163310)
bap$X.track[3] <- mean(c(519837.6,520073.7))
bap$Y.track[3] <- mean(c(5163713,5163310))

sf_final <- rbind(sf_final,bap)

# export ####
#sf_final <- data.frame(sf_final)
#sf_final$geometry <- as.character(sf_final$geometry)
#
#st_write(sf_final %>%
#           select(Station,Echantillo,X.track, Y.track, geometry),
#         paste(here(),"/SIG/SIG_data/centroid_tracks.csv", sep=""),
#         append=FALSE)#,
#         #layer_options = "GEOMETRY=AS_XY")

st_write(sf_final,
         paste(here(),"/SIG/SIG_data/centroid_tracks.csv", sep=""),
         append=FALSE)#,layer_options = "GEOMETRY=AS_XY")
