#https://www.paulamoraga.com/book-spatial/
library(spatstat)

#Homogeneous point process
phom <- rpoispp(lambda = 100,
                win = owin(xrange = c(0, 1), yrange = c(0, 2)))
phom$n
plot(phom, main = "Homogeneous")

#Independent random points
punif <- runifpoint(n = 200,
                    win = owin(xrange = c(0, 1), yrange = c(0, 2)))
punif$n
plot(punif,  main = "Random point")

#Inhomogeneous point process
# intensity function
lambda <- function(x){return(10 + 100 * x[1] + 200 * x[2])}

# grid
xseq <- seq(0, 1, length.out = 50)
yseq <- seq(0, 2, length.out = 100)
grid <- expand.grid(xseq, yseq)

# evaluation of the function on a grid
z <- apply(grid, 1, lambda)

# plot
library(fields)
zmat <- matrix(z, 50, 100)
fields::image.plot(xseq, yseq, zmat, xlab = "x", ylab = "y",
                   main = "lambda(x, y)", asp = 1)

fnintensity <- function(x, y){return(10 + 100 * x + 200 * y)}
pinhom <- rpoispp(lambda = fnintensity,
                  win = owin(xrange = c(0, 1), yrange = c(0, 2)))
pinhom$n
plot(pinhom, main = "Inhomogeneous")



#test
library(sf)
library(here)
shp <- st_read(dsn=paste(here(),"/point_process/Data", sep=""), layer = "grille_zee_spm")
zee <- st_transform(shp, crs = "ESRI:102002")
zee <- zee$geometry
win <- as.owin(zee)

#problem in definition unit area beacause lambda (points per unit area)
plot(phom, main = "Homogeneous")
phom <- rpoispp(lambda = 0.00000001, win = win)
plot(phom, main = "Homogeneous")

fnintensity <- function(x, y){return(0.0000000001 + 0.000000000001 * log(x) + 0.000000002 * log(y))}
pinhom <- rpoispp(lambda = fnintensity,
                  win = win)
plot(pinhom, main = "Inhomogeneous")