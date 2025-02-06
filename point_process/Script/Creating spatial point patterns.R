#https://www.paulamoraga.com/book-spatial/
library(spatstat)

# creation of object of class owin
win <- owin(xrange = c(0, 1), yrange = c(0, 2))
#plot(win)

##method 1
# simulation of 100 points in the window
x <- runif(100, 0, 1)
y <- runif(100, 0, 2)
# creation of the ppp object
X <- ppp(x = x, y = y, window = win)
X

##method 2
X <- runifpoint(n = 100, win = win)

#plot the result
plot(X)
axis(1)
axis(2)

#extract the obs window
Window(X)

#Marks denoting associated information of events
# set marks with marks()
marks(X) <- 1:npoints(X)
# alternatively, set marks with %mark%
X <- X %mark% 1:npoints(X)
plot(X)

#Test whether a set of points lie within a particular observation window
win <- owin() # unit square observation window
marks(X) <- inside.owin(X, w = win)
plot(X)
axis(1)
axis(2)

##From ppp to sf
library(sf)
# ppp object
X <- longleaf
# create data frame with coordinates and marks
ddf <- data.frame(x = X$x, y = X$y, m = marks(X))
# create sf object with data frame and name of coordinates
d <- st_as_sf(ddf, coords = c("x", "y"))

##From sf to ppp
X <- as.ppp(st_coordinates(d), st_bbox(d))
marks(X) <- d$m # alternatively we can use X <- X %mark% d$m
plot(X)

#use a sf polygon
library(rnaturalearth)
library(rnaturalearthdata)
map <- ne_countries(type = "countries", country = "Saint Pierre and Miquelon",
                    scale = "medium", returnclass = "sf")
map <- st_transform(map, crs = "EPSG:29172")
win <- as.owin(map)
X <- runifpoint(100, win = win)
plot(X)

library(here)
shp <- st_read(dsn=paste(here(),"/point_process/Data", sep=""), layer = "grille_zee_spm")
zee <- st_transform(shp, crs = "ESRI:102002")
zee <- zee$geometry
win <- as.owin(zee)
X <- runifpoint(100, win = win)
plot(X)
