# this script compile all the work test for estimate the biomass of 
# sea cucumber to produce the analysis and result for the paper 
# model-based versus design-based

# load packages ####
library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(boot)
library(sdmTMB)
library(INLA)
library(RColorBrewer)
library(XML)

# load data ####
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))
data_holotv <- readRDS(
  paste(here(),"/process_spatiotemp_data/Data/processed/data_abun_tot.rds",
        sep = "")
)
grid_proj <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_bathy.rds",
        sep="")
)

# prepare data ####
# synchronize variable name
names(grid_proj)[1:2] <- c("X","Y")

# keep only station in the study area
data_holotv <- 
  st_join(calcul_area,data_holotv) # to get intersection of points and poly

# convert abundance in biomass
data_holotv$biomass <- data_holotv$abun*0.4
data_holotv$density.t_km2 <- (data_holotv$biomass/1000)/(data_holotv$area/1e6)

# convert distance from meter to kilometer
data_holotv$X <- data_holotv$X/1000
data_holotv$Y <- data_holotv$Y/1000
data_holotv$area <- data_holotv$area/1e6

# keep only the column needed
data_estimate <- as.data.frame(data_holotv) %>% 
  select(station,X,Y,density.t_km2,year)
grid_proj <- as.data.frame(grid_proj)[,1:2]

# extract the surface of the study area
stock_surface <- as.numeric(st_area(calcul_area)/1e6)

###-###-###-###-###-###-###-###-###
# Show survey sampling design ####
###-###-###-###-###-###-###-###-###
#import the ZEE limitation
zee <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""), layer = "eez")
zee <- st_transform(zee, crs = 4467)
zee <- zee$geometry
#import the ZEE grid
shp_grid <- st_read(dsn=paste(here(),"/SIG/SIG_Data", sep=""),
                    layer = "grille_zee_spm")
shp_grid <- st_transform(shp_grid, crs = 4467)
shp_grid <- shp_grid %>% mutate(number = c(1:nrow(shp_grid)))

shp_grid_study <- st_centroid(shp_grid) %>%
  dplyr::filter(!is.na(x = Echantillo)) %>%
  mutate(X = st_coordinates(.)[,1]) %>%
  mutate(Y = st_coordinates(.)[,2]) 
shp_grid_study <- 
  st_join(calcul_area,shp_grid_study)
shp_grid_study$geometry <- shp_grid$geometry[
  shp_grid$number %in% shp_grid_study$number]

source(paste(here(),"/via3_data_exploration/fct_WG84_WGSPM.R", sep=""))

## import individual position
data_position <- readRDS(
  paste(here(),"/point_process/Data/data_position_2025_epsg4461.rds",sep=""))

### get meters ###
# Using a sample point in London
df <- data.frame(lon = -56, lat = 45.8)
df_sf <- st_as_sf(df, coords = c("lon", "lat"), crs = 4467)
df_projected <- st_transform(df_sf, crs = 3857) # 3857 is Web Mercator (Meters)
rm(df, df_sf)

## creation of the map ####
maps_tracks <- ggplot()+
  geom_sf(data = shp_grid_study, color = "black",
          fill = "#00000000")+
  geom_sf(data = zee, color = "red", fill = NA, linewidth = 1,
          linetype = "longdash", show.legend = FALSE)+
  geom_point(data=data_position, aes(X,Y), color="black")+
  coord_sf(xlim = c(max(shp_grid_study$X)+1852*2,min(shp_grid_study$X)-1852*2),
           ylim = c(max(shp_grid_study$Y)+1852*2,min(shp_grid_study$Y)-1852*2), 
           expand = FALSE, crs = 4467)+
  scale_x_continuous(breaks = seq(-56.40, -56.15, by = 0.125)) +
  annotate("text", x = 544500, y = 5074897, colour = "red",
           size = 4, label = "ZEE")+
  annotate("text", x = 568000, y = 5021800, colour = "darkblue",
           size = 2.2, label = "IFREMER")+
  labs(x = "", y = "", colour = "", breaks = c("1", "3"))+#, caption = "IFREMER")+
  theme(title = element_text(color = "black",face = "bold"),
        panel.border = element_blank(),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.major = element_line(colour = "lightblue"),
        panel.grid.minor = element_line(colour = "lightblue"),
        plot.margin = margin(t = 10,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0), # Left margin
        )+
  annotation_scale(location = "bl", line_width = .3) +
  annotation_north_arrow(location = "tr", height = unit(0.5, "cm"),
                         width = unit(0.5, "cm"))


###-###-###-###-###-###-###-###
# Design-based estimation ####
###-###-###-###-###-###-###-###

## the arithmetic mean ####
bio_tot <- function(donnees, indices) {
  biomasse_carre <- donnees[indices]  # in ton
  biomasse_totale <- mean(biomasse_carre)*stock_surface
  return(biomasse_totale)
}
biomass_stock_2021 <- mean(
  data_estimate$density.t_km2[grep(2021,data_estimate$year)]
  )*stock_surface
biomass_stock_2022 <- mean(
  data_estimate$density.t_km2[grep(2022,data_estimate$year)]
)*stock_surface
biomass_stock_2023 <- mean(
  data_estimate$density.t_km2[grep(2023,data_estimate$year)]
)*stock_surface
biomass_stock_2025 <- mean(
  data_estimate$density.t_km2[grep(2025,data_estimate$year)]
)*stock_surface

## Variances with bootstrap ####
N <- 10000 # number of replica to define
# execute the bootstrap
set.seed(13) # for reproducibility
bootobject_2021 <- boot(
  data = data_estimate$density.t_km2[grep(2021,data_estimate$year)],
  statistic = bio_tot, R = N)
bootobject_2022 <- boot(
  data = data_estimate$density.t_km2[grep(2022,data_estimate$year)],
  statistic = bio_tot, R = N)
bootobject_2023 <- boot(
  data = data_estimate$density.t_km2[grep(2023,data_estimate$year)],
  statistic = bio_tot, R = N)
bootobject_2025 <- boot(
  data = data_estimate$density.t_km2[grep(2025,data_estimate$year)],
  statistic = bio_tot, R = N)

## extract result ####
boot_result <- data.frame()

extract_boot_stat <- function(bootobject){
  boot_bio <- as.numeric(bootobject[["t"]])
  median_value <- median(boot_bio)
  quantile_5 <- quantile(boot_bio, 0.05)
  quantile_95 <- quantile(boot_bio, 0.95)
  mean_value <- mean(boot_bio)
  return(data.frame(mean_value,median_value,quantile_5,quantile_95))
}

boot_result <- rbind(boot_result,extract_boot_stat(bootobject_2021))
boot_result <- rbind(boot_result,extract_boot_stat(bootobject_2022))
boot_result <- rbind(boot_result,extract_boot_stat(bootobject_2023))
boot_result <- rbind(boot_result,extract_boot_stat(bootobject_2025))
row.names(boot_result) <- c(2021,2022,2023,2025)

## show result ####
# Créer l'histogramme avec ggplot2
ggplot(data.frame(Biomass = as.numeric(bootobject_2025[["t"]])), aes(x = Biomass)) +
  geom_histogram(aes(y = after_stat(density)), bins = 30, fill = "lightblue", color = "black") +
  geom_density(col='red') +
  geom_vline(aes(xintercept = boot_result["2025",]$median_value), 
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = boot_result["2025",]$quantile_5), 
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = boot_result["2025",]$quantile_95), 
             color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(aes(xintercept = boot_result["2025",]$mean_value), 
             color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Histogramme des Biomasses obtenues par Bootstrap",
       x = "Biomasse",
       y = "Densité") +
  scale_x_continuous(labels = scales::label_comma())+
  theme_minimal()


###-###-###-###-###-###-###-###
# Model-based estimation ####
###-###-###-###-###-###-###-###

## creation of the mesh with R-INLA ####
#using the grid of the area
bnd <- INLA::inla.nonconvex.hull(cbind(grid_proj$long, grid_proj$lat),
                                 convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(grid_proj$long, grid_proj$lat),
                                 convex = -0.15)

mesh_inla <- INLA::inla.mesh.2d(
  loc = as.matrix(as.data.frame(grid_proj)[,c(1,2)]),
  boundary = list(bnd,bnd2),
  cutoff = 1.852, # minimum triangle edge length
  max.edge = c(3*1.852, 20*1.852), # inner and outer max triangle lengths
) # 1.852 is the distance in meter of a mile nautic

mesh_2021 <- make_mesh(data_estimate %>% filter(year==2021),
                       c("X", "Y"), mesh = mesh_inla)
mesh_2022 <- make_mesh(data_estimate %>% filter(year==2022),
                       c("X", "Y"), mesh = mesh_inla)
mesh_2023 <- make_mesh(data_estimate %>% filter(year==2023),
                       c("X", "Y"), mesh = mesh_inla)
mesh_2025 <- make_mesh(data_estimate %>% filter(year==2025),
                       c("X", "Y"), mesh = mesh_inla)

plot(mesh_2021$mesh, main = NA, edge.color = "grey60", asp = 1)
points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)

## Write and fit the model ####
fit_2021 <- sdmTMB(
  density.t_km2 ~ 1,
  data = data_estimate %>% filter(year==2021),
  mesh = mesh_2021,
  family = Gamma(link = "log"),
  spatial = "on")

fit_2022 <- sdmTMB(
  density.t_km2 ~ 1,
  data = data_estimate %>% filter(year==2022),
  mesh = mesh_2022,
  family = Gamma(link = "log"),
  spatial = "on")

fit_2023 <- sdmTMB(
  density.t_km2 ~ 1,
  data = data_estimate %>% filter(year==2023),
  mesh = mesh_2023,
  family = Gamma(link = "log"),
  spatial = "on")

fit_2025 <- sdmTMB(
  density.t_km2 ~ 1,
  data = data_estimate %>% filter(year==2025),
  mesh = mesh_2025,
  family = Gamma(link = "log"),
  spatial = "on")

## check and validate the model ####
sanity(fit_2021)
sanity(fit_2022)
sanity(fit_2023)
sanity(fit_2025)
# check quantile residuals
data_estimate$resids <- NA
for (i in 1:4){
  fit <- list(fit_2021,fit_2022,fit_2023,fit_2025)[[i]]
  year <- c(2021,2022,2023,2025)[i]
  data_estimate$resids[grep(year,data_estimate$year)] <- residuals(fit) 
  # randomized quantile residuals
  qqnorm(data_estimate$resids)
  qqline(data_estimate$resids)
}
# check residual pattern
ggplot(data_estimate, aes(X, Y, col = resids)) +
  scale_colour_gradient2() +
  geom_point() +
  facet_wrap(~year, nrow=1)+
  coord_fixed()

## extract the biomass estimation ####
mod_result <- data.frame()
map_predict <- data.frame()
for (i in 1:4){
  fit <- list(fit_2021,fit_2022,fit_2023,fit_2025)[[i]]
  year <- c(2021,2022,2023,2025)[i]
  p_sdm <- predict(fit, newdata = grid_proj, return_tmb_object = TRUE,
                   se_fit = TRUE)
  index <- get_index(p_sdm, area = 0.25, bias_correct = TRUE)
  mod_result <- rbind(mod_result,index)
  map_predict <- rbind(map_predict, data.frame(date=year,p_sdm$data))
}
mod_result$`_sdmTMB_time`<-c(2021,2022,2023,2025)


###-###-###-###-###-###-###-###
# Show result of biomass ####
###-###-###-###-###-###-###-###
compare_biomass <- data.frame(
  date = c(2021,2021.2,2022,2022.2,2023,2023.2,2024,2025,2025.2),
  methods = c(rep(c("Design","Model"),3),NA,c("Design","Model")),
  biomass = NA,
  lwr = NA,
  upr = NA
)

for (i in 1:nrow(compare_biomass)){
  if (!is.na(compare_biomass$methods[i])){
    if (compare_biomass$methods[i] == "Design"){
      compare_biomass$biomass[i]<-boot_result[
        as.character(compare_biomass$date[i]),"mean_value"]
      compare_biomass$lwr[i]<-boot_result[
        as.character(compare_biomass$date[i]),"quantile_5"]
      compare_biomass$upr[i]<-boot_result[
        as.character(compare_biomass$date[i]),"quantile_95"]
    }
    
    if (compare_biomass$methods[i] == "Model"){
      compare_biomass$biomass[i]<-mod_result$est[
        grep(trunc(compare_biomass$date[i]),mod_result$`_sdmTMB_time`)]
      compare_biomass$lwr[i]<-mod_result$lwr[
        grep(trunc(compare_biomass$date[i]),mod_result$`_sdmTMB_time`)]
      compare_biomass$upr[i]<-mod_result$upr[
        grep(trunc(compare_biomass$date[i]),mod_result$`_sdmTMB_time`)]
    } 
  }
}

ggplot(compare_biomass,aes(x = date, y = biomass,colour = methods))+
  geom_point(show.legend = T)+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2,show.legend = F)+
  labs(x = "", y = "Biomass (tonnes)", title = "Biomass from survey")+
  theme_bw()

saveRDS(compare_biomass,paste(
  here(),"/comparing_biomass_method/Data/compare_biomass_estimation.rds",sep=""
))


ggplot(map_predict, aes(X, Y, fill = est)) +
  geom_raster() +
  scale_fill_viridis_c() +
  facet_wrap(~date,nrow=1)+
  coord_fixed(expand = FALSE)+
  theme(aspect.ratio = 3)

ggplot(map_predict, aes(X, Y, fill = est_se)) +
  geom_raster() +
  scale_fill_gradient2() +
  facet_wrap(~date,nrow=1)+
  coord_fixed(expand = FALSE)+
  theme(aspect.ratio = 3)
