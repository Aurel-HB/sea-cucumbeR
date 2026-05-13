# This code create a first approach to geostatistic with gstlearn 
# for biomass with abunx0.4


library(gstlearn)
library(ggplot2)
library(here)
library(dplyr)

## Load data ####
data_abun <- readRDS(
  paste(here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
        sep=""))

shp_grid <-readRDS(
  paste(here(),"/point_process/Output/shp_grid.rds", sep=""))

mean_weight <- 0.4 # weight of an adult sea cucumber in kg
surfsquare <- 1.852*1.852 # in kilometer square


## Prepare data ####
datcsv <- as.data.frame(data_abun)%>% select(X,Y,intensity)
datcsv$intensity <- datcsv$intensity*mean_weight #intensity in kg/m²
filecsv <- "Density_2025.csv"
write.csv(datcsv,
          paste(here(),"/comparing_biomass_method/Data/",filecsv,
                sep=""), row.names = FALSE)
# creating Db oject for gstlearn analysis
csv <- CSVformat_create(flagHeader=TRUE, naString = "MISS")
dat <- Db_createFromCSV(paste(here(),"/comparing_biomass_method/Data/",filecsv,
                              sep=""), csv=csv)
dat$setLocators(c("X","Y"), ELoc_X())
dat$setLocator("intensity", ELoc_Z(), cleanSameLocator=TRUE)
# plot the Db
ggplot(as.data.frame(dat[]))+
  geom_point(aes(x=X,y=Y, size = intensity))

# We define a projection, based on the center of gravity of the data.
# Modify the coordinates of the Data Base and of the Polygon consequently.
projec = Projection(TRUE, dat)
err = projec$operateOnDb(dat)
err = projec$operateOnPolygons(dat)

# Display in the projected system
ggplot(as.data.frame(dat[]))+
  geom_point(aes(x=X,y=Y, size = intensity))

## calculating the Spatial Indices####
# Instantiate the class for calculating the Spatial Indices
cg = SpatialIndices(dat)

# Calculate the center of gravity, the Inertia and the Isotropy of the densities
err = cg$computeCGI("intensity")
cgAxesI0 = cg$getAxes()
# Get the results of the calculation of the Center of gravity 
# (the calculations are performed in the projected space). 
# Results are back transformed into the initial system prior to displaying them.
vec = projec$operateInvert(cg$getCenter())
cat("Center of gravity = ",vec,"\n")
cat("Inertia =           ",cg$getInertia(),"\n")
cat("Isotropy =          ",cg$getIso(),"\n")
# Compute the center of gravity of the samples 
# (without considering any specific variable).
err = cg$computeCGI("")
cgAxesP0 = cg$getAxes()
vec = projec$operateInvert(cg$getCenter())
cat("Center of gravity = ",vec,"\n")
cat("Inertia =           ",cg$getInertia(),"\n")
cat("Isotropy =          ",cg$getIso(),"\n")
# Display the sample within the Polygon, together with the two inertia 
# calculated beforehand: the one using the sample locations only (in blue) and 
# the one using the variable intensity into account (in red).
ggplot(as.data.frame(dat[]))+
  geom_point(aes(x=X,y=Y, size = intensity))+
  geom_line(data = data.frame(X=cgAxesI0[[1]], Y=cgAxesI0[[2]]),
            aes(x = X, y = Y),linewidth=1, color="red")+
  geom_line(data = data.frame(X=cgAxesI0[[3]], Y=cgAxesI0[[4]]),
            aes(x = X, y = Y),linewidth=1, color="red")+
  geom_line(data = data.frame(X=cgAxesP0[[1]], Y=cgAxesP0[[2]]),
            aes(x = X, y = Y),linewidth=1, color="blue")+
  geom_line(data = data.frame(X=cgAxesP0[[3]], Y=cgAxesP0[[4]]),
            aes(x = X, y = Y),linewidth=1, color="blue")+
  labs(title = "Center of gravity and Inertia of densities and Samples")

# Calculate: - the Abundance index - the Positive area - the Spreading area and
# the Equivalent area
err = cg$spatial("intensity")
VDD = cg$getQT("intensity")
ggplot() +
  geom_point(data = data.frame(X=VDD[[1]], Y=VDD[[2]]),
            aes(x = X, y = Y), color="black") +
  geom_line(data = data.frame(X=VDD[[1]], Y=VDD[[2]]),
            aes(x = X, y = Y),linewidth=1, color="black")+
  labs(title = "Curve: (Q - Q(T)) / Q",
       x = "T", y = "(Q-Q(T))/Q")

# Calculate the patches and represent them by color.
centers = cg$getPatches("intensity", 500000, 10)
ggplot(as.data.frame(dat[])) + 
  geom_point(aes(x=X,y=Y, size = intensity, color = as.factor(Patch)))

## calculate influence area ####
# Calculer le diagramme de Voronoï
library(deldir)
library(sf)
stations <- data.frame(
     x = data_abun$X,
     y = data_abun$Y
)

voronoi <- deldir(stations$x, stations$y,
                  rw=c(min(stations$x)-1,
                       max(stations$x)+1,
                       min(stations$y)-1,
                       max(stations$y))+1)

# Extraire les polygones de Voronoï
tiles <- tile.list(voronoi)

# Fonction pour convertir un polygone de Voronoï en un polygone sf
tile_to_sf <- function(tile) {
  # Extraire les coordonnées du polygone
  coords <- cbind(tile$x, tile$y)
  # Ajouter la première coordonnée à la fin pour fermer le polygone
  coords <- rbind(coords, coords[1, ])
  # Créer un polygone sf
  st_polygon(list(coords))
}

# Appliquer la fonction à chaque polygone de Voronoï
voronoi_sf <- st_as_sfc(lapply(tiles, tile_to_sf))

# Convertir en un objet sf avec une colonne d'identifiant
voronoi_sf <- st_as_sf(voronoi_sf)

# Calculer l'aire de chaque polygone (en unités carrées)
voronoi_sf$area <- st_area(voronoi_sf)

# Visualisation
plot(st_geometry(voronoi_sf))
points(stations$x, stations$y, col = "red", pch = 19)


## Prepare data ####
datcsv$area <- as.numeric(voronoi_sf$area)
datcsv$intensity <-datcsv$intensity /1000 #t/m²
filecsv <- "Density_2025.csv"
write.csv(datcsv,
          paste(here(),"/comparing_biomass_method/Data/",filecsv,
                sep=""), row.names = FALSE)
# creating Db oject for gstlearn analysis
csv <- CSVformat_create(flagHeader=TRUE, naString = "MISS")
dat <- Db_createFromCSV(paste(here(),"/comparing_biomass_method/Data/",filecsv,
                              sep=""), csv=csv)
dat$setLocators(c("X","Y"), ELoc_X())
dat$setLocator("intensity", ELoc_Z(), cleanSameLocator=TRUE)
dat$setLocator("area", ELoc_W())
# plot the Db
ggplot(as.data.frame(dat[]))+
  geom_point(aes(x=X,y=Y, size = intensity))

# We define a projection, based on the center of gravity of the data.
# Modify the coordinates of the Data Base and of the Polygon consequently.
projec = Projection(TRUE, dat)
err = projec$operateOnDb(dat)
err = projec$operateOnPolygons(dat)

# Display in the projected system
ggplot(as.data.frame(dat[]))+
  geom_point(aes(x=X,y=Y, size = intensity))

## calculating the Spatial Indices####
# Instantiate the class for calculating the Spatial Indices
cg = SpatialIndices(dat)

# Calculate: - the Abundance index - the Positive area - the Spreading area and
# the Equivalent area
err = cg$spatial("intensity")
VDD = cg$getQT("intensity")
ggplot() +
  geom_point(data = data.frame(X=VDD[[1]], Y=VDD[[2]]),
             aes(x = X, y = Y), color="black") +
  geom_line(data = data.frame(X=VDD[[1]], Y=VDD[[2]]),
            aes(x = X, y = Y),linewidth=1, color="black")+
  labs(title = "Curve: (Q - Q(T)) / Q",
       x = "T", y = "(Q-Q(T))/Q")
