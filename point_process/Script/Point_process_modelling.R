# This code test an approach of modelling on point process to characterize it

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(fields)

#############
#load data
#############

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

HOLOTVSPM2025 <- readRDS(paste(here(),
              "/comparing_biomass_method/Data/03_HOLOTVSPM2025_summary.rds",
              sep=""))

########################
# Start with one PPP
########################

PPP <- list_PPP[[1]]
PPP[["window"]][["units"]] <- list("metre","metres")
PPP[["window"]][["yrange"]] <- PPP[["window"]][["yrange"]] - PPP[["window"]][["yrange"]][[1]]
stn <- names(list_PPP)[1]


##############
# math recall
##############
# Reflection matrix:
#    Reflection(θ)=[cos2θ   sin2θ
#                   sin2θ  −cos2θ]

# Rotation matrix:
#    Rotation(θ)=[cosθ  -sinθ
#                 sinθ   cosθ]
##############


# spatialization of the PPP ####
start_coord <- HOLOTVSPM2025[grep(96,HOLOTVSPM2025$STN),
                             c("X.haul_start","Y.haul_start")]
end_coord <- HOLOTVSPM2025[grep(96,HOLOTVSPM2025$STN),
                           c("X.shoot_end","Y.shoot_end")]

affine <-(end_coord$Y.shoot_end - start_coord$Y.haul_start)/
  (end_coord$X.shoot_end - start_coord$X.haul_start)

ordonnee <- mean(c(end_coord$Y.shoot_end - affine*end_coord$X.shoot_end,
                   start_coord$Y.haul_start - affine*start_coord$X.haul_start))

Cartesian_equation <- function(x,affine,ordonnee){
  y <- affine*x+ordonnee
  return(y)
}

# rotation of angle theta
theta <- atan(start_coord$X.haul_start/start_coord$Y.haul_start)
new_coord <- as.matrix(data.frame(X=PPP[["x"]],Y=PPP[["y"]]))
new_coord <- spdep::Rotation(new_coord,-theta) #sens horaire donc signe négatif

# translation at the start point of the haul
matrix.start <- t(matrix(start_coord,2,316))
matrix.start <- as.matrix(data.frame(
  X=rep(start_coord$X.haul_start,dim(new_coord)[1]),
  Y=rep(start_coord$Y.haul_start,dim(new_coord)[1])
  ))
new_coord <- new_coord + matrix.start





# statistic summary ####
# spatial inhomogeneity : kernel smoothed estimate of intensity 
D1 <- density.ppp(PPP,bw.ppl(PPP))
D2 <- density.ppp(PPP,bw.diggle(PPP))
dX <- density(PPP, sigma=1, at="points")
# Quadrat counting test of homogeneity
tS <- quadrat.test(PPP)
tS <- quadrat.test(PPP,1,5)
tS <- quadrat.test(PPP,1,10)
tS
# Variance‑to‑Mean Ratio (VMR)
Q   <- quadratcount(PPP, nx = 1, ny = 10)
N   <- as.vector(Q)
VMR <- var(N) / mean(N)
VMR
#> VMR ≈ 1 → CSR; VMR ≫ 1 → strong clustering/inhomogeneity; VMR < 1 → regularity.


# inference statistic test ####
#Poisson model
fit0 <- ppm(PPP~1)
fit1 <- ppm(PPP~x + y)
#Likelihood Ratio Test
anova(fit0, fit1, test="LR")
#This example shows the likelihood ratio test of the null hypothesis of CSR 
#against the alternative of an inhomogeneous Poisson process with intensity that 
#is a loglinear function of the coordiantes. The p-value, shown under 
#the heading Pr(>Chi), is extremely small, indicating rejection of  CSR in favour 
#of the alternative. Note that 2e-16 or 2 × 10−16 is the smallest detectable 
#difference  between ‘real numbers’ on the 32-bit computer which produced this 
#output, so the output says that the p-value is effectively zero.
AIC(fit0);AIC(fit1)

X <- simulate(fit1); plot(X[[1]])
lambda <- function(x){return(
  fit1$coef[[1]]+fit1$coef[[2]]*x[1]+fit1$coef[[3]]*x[2]
)}
xseq <- seq(0, 1.5, length.out = 50)
yseq <- seq(60, 660, length.out = 6000)
grid <- expand.grid(xseq, yseq)
z <- apply(grid, 1, lambda)
zmat <- matrix(z, 50, 6000)
#library(fields)
image.plot(xseq, yseq, zmat, xlab = "x", ylab = "y",main = "lambda(x)")

#select best model for intensity function
bigfit <- ppm(PPP ~ polynom(x,y,3))
formula(bigfit)
goodfit <- step(bigfit, trace=0)
formula(goodfit)
AIC(goodfit)
# simulation of the process on another surface
# 1) extract the intensity function
lambda <- function(x){return(
  goodfit$coef[[1]]+
    goodfit$coef[[2]]*x[1]+
    goodfit$coef[[3]]*x[2]+
    goodfit$coef[[4]]*x[1]^2+
    goodfit$coef[[5]]*x[2]^2+
    goodfit$coef[[6]]*x[1]^3+
    goodfit$coef[[7]]*x[2]^3
)}
#2) create a grid of the intensity function on the surface of simulation
xseq <- seq(0, 1.5, length.out = 50)
yseq <- seq(60, 660, length.out = 6000)
grid <- expand.grid(xseq, yseq)
z <- apply(grid, 1, lambda)
zmat <- matrix(z, 50, 6000)
image.plot(xseq, yseq, zmat, xlab = "x", ylab = "y",main = "log(lambda(x))")
ggplot(data = data_position %>% filter(station==stn))+
  geom_point(aes(x=X,y=Y))+
  ylab("")
#3) create the new point process pattern
fnintensity <- function(x, y){return(
  exp(goodfit$coef[[1]]+
    goodfit$coef[[2]]*x+
    goodfit$coef[[3]]*y+
    goodfit$coef[[4]]*x^2+
    goodfit$coef[[5]]*y^2+
    goodfit$coef[[6]]*x^3+
    goodfit$coef[[7]]*y^3)
)}
pinhom <- rpoispp(lambda = fnintensity,
                  win = owin(xrange = c(0, 1.5), yrange = c(60,660)))
plot(x = pinhom$x,y = pinhom$y)
#plot(pinhom, main = "Inhomogeneous")




# Cox model
fitox0 <- kppm(PPP~1, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
fitox1 <- kppm(PPP~x + y, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
AIC(fitox0);AIC(fitox1)
#simulation
X <- simulate(fitox0); plot(X[[1]])
plgcp <- rLGCP("matern", fitox0$lambda, var = fitox0$clustpar[1],
               scale = fitox0$clustpar[2], nu=fitox0$covmodel$margs$nu,
               win = owin(xrange = c(0, 150), yrange = c(60,660)))
plot(attr(plgcp, "Lambda"))
plgcp <- rLGCP("matern", fitox0$lambda, var = fitox0$clustpar[1],
               scale = fitox0$clustpar[2], nu=fitox0$covmodel$margs$nu,
               win = owin(xrange = c(0, 1.5), yrange = c(60,660)))
plot(x = plgcp$x,y = plgcp$y)
m <- as.im(function(x,y){
  exp(
    fitox1$po$coef[1]+
      x*fitox1$po$coef[2]+
      y*fitox1$po$coef[3]
  )
}, W=owin(xrange = c(0, 1.5), yrange = c(60,660)))
plgcp <- rLGCP("matern", m,
               var = fitox1$clustpar[1],
               scale = fitox1$clustpar[2], nu=fitox1$covmodel$margs$nu,
               win = owin(xrange = c(0, 1.5), yrange = c(60,660)))
plot(x = plgcp$x,y = plgcp$y)

fitox2 <- kppm(PPP~1, clusters = "LGCP", method="clik2", model="gauss")
plgcp <- rLGCP("gauss", fitox2$lambda, var = fitox2$clustpar[1],
               scale = fitox2$clustpar[2],
               win = owin(xrange = c(0, 1.50), yrange = c(60,160)))
plot(attr(plgcp, "Lambda"))
plgcp <- rLGCP("gauss", fitox2$lambda, var = fitox2$clustpar[1],
               scale = fitox2$clustpar[2],
               win = owin(xrange = c(0, 1.50), yrange = c(160,260)))
plot(attr(plgcp, "Lambda"))
