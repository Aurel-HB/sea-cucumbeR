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
library(gridExtra)

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

data_position_substrat <- readRDS(paste(here(),
              "/via3_data_exploration/Data/substrat/data_position_substrat.rds",
              sep=""))

data_substrate <- readRDS(paste(here(),
              "/via3_data_exploration/Data/substrat/data_substrat_2025.rds",
              sep=""))

########################
# Start with one PPP
########################

PPP <- list_PPP[[1]] # PPP<-list_PPP[["179"]] PPP<-list_PPP[["149"]]
# PPP<-list_PPP[["143"]]
PPP[["window"]][["units"]] <- list("metre","metres")
#PPP[["window"]][["yrange"]] <- PPP[["window"]][["yrange"]] - PPP[["window"]][["yrange"]][[1]]
stn <- names(list_PPP)[1]
# stn=149 stn=179 stn=143

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
plot(D1);plot(D2);plot(dX)
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
#Poisson model ####
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
Show_result <- rbind(
  as.data.frame(cbind(X$`Simulation 1`$x,X$`Simulation 1`$y,
                      rep("X",X$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  #ylim(PPP$window$yrange)+
  #theme(axis.ticks = element_blank(),
  #      axis.text = element_blank())
  theme()

# extract point process to use it on a another surface
#lambda <- function(x){return(
#  fit1$coef[[1]]+fit1$coef[[2]]*x[1]+fit1$coef[[3]]*x[2]
#)}
#xseq <- seq(0, 1.5, length.out = 50)
#yseq <- seq(60, 660, length.out = 6000)
#grid <- expand.grid(xseq, yseq)
#z <- apply(grid, 1, lambda)
#zmat <- matrix(z, 50, 6000)
##library(fields)
#image.plot(xseq, yseq, zmat, xlab = "x", ylab = "y",main = "lambda(x)")

#select best model for intensity function
bigfit <- ppm(PPP ~ polynom(x,y,3))
formula(bigfit)
goodfit <- step(bigfit, trace=0)
formula(goodfit)
AIC(goodfit)
anova(fit1, goodfit, test="LR")

X <- simulate(goodfit); plot(X[[1]])
Show_result <- rbind(
  as.data.frame(cbind(X$`Simulation 1`$x,X$`Simulation 1`$y,
                      rep("X",X$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

# simulation of the process on another surface 
# 1) extract the intensity function
#lambda <- function(x){return(
#  goodfit$coef[[1]]+
#    goodfit$coef[[2]]*x[2]+
#    goodfit$coef[[3]]*x[2]^2+
#    goodfit$coef[[4]]*x[1]^3+
#    goodfit$coef[[5]]*x[2]^3
#)}
#2) create a grid of the intensity function on the surface of simulation
#xseq <- seq(0, 1.5, length.out = 50)
#yseq <- seq(60, 660, length.out = 6000)
#grid <- expand.grid(xseq, yseq)
#z <- apply(grid, 1, lambda)
#zmat <- matrix(z, 50, 6000)
#image.plot(xseq, yseq, zmat, xlab = "x", ylab = "y",main = "log(lambda(x))")
#ggplot(data = data_position %>% filter(station==stn))+
#  geom_point(aes(x=X,y=Y))+
#  ylab("")
#3) create the new point process pattern
#fnintensity <- function(x, y){return(
#  exp(goodfit$coef[[1]]+
#    goodfit$coef[[2]]*y+
#    goodfit$coef[[3]]*y^2+
#    goodfit$coef[[4]]*x^3+
#    goodfit$coef[[5]]*y^3)
#)}
#pinhom <- rpoispp(lambda = fnintensity,
#                  win = owin(xrange = c(0, 1.5), yrange = c(60,660)))
#plot(x = pinhom$x,y = pinhom$y)
#plot(pinhom, main = "Inhomogeneous")


# Check the model ####
#Envelope Tests (K-Function)
fit_checked <- goodfit
# Simulate 99 replicates from the fitted model
n_sim <- 99
simulated_patterns <- simulate(fit_checked,n_sim)
# Compute K-functions for observed and simulated data
K_obs <- Kest(PPP)
K_sim <- Kest(simulated_patterns[[1]])
for (i in 2:length(simulated_patterns)){
  K_sim <- rbind(K_sim,Kest(simulated_patterns[[i]]))
}
# Plot envelope test
plot(K_obs, main = "Envelope Test for Poisson Model")
plot(K_sim, add = TRUE, col = "gray", lty = 2)
plot(K_obs,add = TRUE, main = "Envelope Test for Poisson Model")
legend("topleft", legend = c("Observed", "Simulated"), lty = c(1, 2),
       col = c("black", "gray"))

#QQ-Plots for Nearest Neighbor Distances
# Compute distances to nearest neighbor for observed and simulated data
nn_obs <- nndist(PPP)
nn_sim <- matrix(unlist(lapply(1:n_sim, function(i) nndist(simulated_patterns[[i]]))), ncol = n_sim)
# Flatten simulated distances and sort
nn_sim_flat <- sort(as.vector(nn_sim))
# Sort observed distances
nn_obs_sorted <- sort(nn_obs)
# QQ-Plot
plot(qnorm(ppoints(length(nn_obs_sorted))), nn_obs_sorted,
     main = "QQ-Plot: Observed vs. Theoretical Quantiles",
     xlab = "Theoretical Quantiles (Normal)", ylab = "Observed Quantiles")
abline(0, 1, col = "red")
# Ensure lengths match by resampling the observed distances
nn_obs_sorted <- sort(sample(nn_obs, length(nn_sim_flat), replace = TRUE))
# QQ-Plot: Observed vs. Simulated (empirical quantiles)
plot(nn_sim_flat, nn_obs_sorted,
     main = "QQ-Plot: Observed vs. Simulated Nearest Neighbor Distances",
     xlab = "Simulated Quantiles", ylab = "Observed Quantiles")
abline(0, 1, col = "red")

#Residual Analysis (Pearson Residuals)
# Extract fitted intensity at data points
fitted_intensity <- predict(fit_checked)
# Residuals for unmarked patterns: log(fitted_intensity) - log(observed intensity)
# Since observed intensity is 1 at each point, use:
residuals <- log(fitted_intensity + 1e-10) - 0
#plot residual
plot(residuals, pch = ".", main = "Spatial Residuals")
data_resid <- data.frame()
for (i in 1:residuals$dim[1]){
  data_resid <- as.data.frame(rbind(data_resid,
                                    data.frame(
                                    "x"=rep(residuals$xcol[i],residuals$dim[1]),
                                    "y"=residuals$yrow,
                                    "r"=as.data.frame(residuals[["v"]])[,i]
                                    )))
}
ggplot(data = data_resid)+
  geom_point(aes(x=x,y=y,color=r))+
  scale_color_continuous(type = "viridis")+
  ylab("")+xlab("")+labs(title = "Spatial Residuals")+
  theme()


# Cox model ####
fitox0 <- kppm(PPP~1, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
formule <- formula(goodfit)
fitox1 <- kppm(PPP,formule, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
fitox1 <- kppm(PPP~y, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
AIC(fitox0);AIC(fitox1)
#simulation
X_LGCP <- simulate(fitox0); plot(X[[1]])
Show_result <- rbind(
  as.data.frame(cbind(X$`Simulation 1`$x,X$`Simulation 1`$y,
                      rep("X",X$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

# extract point process to use it on a another surface
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


### Locally_weighted Poisson point process model ####
pp <- PPP
# 1) find an estimate of lambda Poisson
fit <- ppm(pp ~ polynom(x,y,3))
fit <- step(fit, trace=0)
# 2) calculate the local version of the K-function
localKfunction <- localK(pp, verbose = TRUE)
# 3) find an estimate of the interaction term
theo <- localKfunction$theo
nX <- npoints(pp)
r <- localKfunction$r
dr <- c(diff(r), diff(r)[nX - 1])
l1 <- vector(l = nX)
for(i in 1:nX){
  progressreport(i,nX)
  l1[i] <- exp(sum((sign(localKfunction[[i]] - theo) * (
    (localKfunction[[i]] - theo) ^ 2)) / theo * dr, 
                   na.rm = TRUE))
}
pp1 <- pp
pp1$marks <- l1
# 4) smooth the surface with the Nadaraya-Watson smoother
im0_kernel <- Smooth(pp1, sigma = bw.CvL(pp1))
im0_kernel <- Smooth(pp1, sigma = bw.CvL(pp1), positive = TRUE)
im0_kernel <- Smooth(pp1, sigma = bw.diggle(pp1), positive = TRUE)
#bw.CvL=Cronie and van Lieshout's 
#Criterion for Bandwidth Selection for Kernel Density

#Log-transform the kernel (if all values are positive)
im0_kernel_log <- log(im0_kernel + 1e-6)
# Standardizing to avoid extreme values 
im0_kernel_log_scaled <- im0_kernel_log
im0_kernel_log_scaled[["v"]] <- scale(im0_kernel_log[["v"]])
summary(im0_kernel_log_scaled)

# 5) Fit the final model, maintaining the same covariate dependence specification
formula <- as.formula(paste("pp",as.character(fit$trend)[1],
                            as.character(fit$trend)[2],
                            " + offset(im0_kernel_log_scaled)",sep="")) 
lwppp <- ppm(formula)
#lwppp <- ppm(pp~y+offset(im0_kernel_log_scaled))
AIC(lwppp)

X_inhom <- simulate(fit)
X_lwppp <- simulate(lwppp)
Show_result <- rbind(
  as.data.frame(cbind(X_lwppp$`Simulation 1`$x,X_lwppp$`Simulation 1`$y,
                      rep("X_lwppp",X_lwppp$`Simulation 1`$n))),
  as.data.frame(cbind(X_inhom$`Simulation 1`$x,X_inhom$`Simulation 1`$y,
                      rep("X_inhom",X_inhom$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())


# point process with substrate as covariate ####
ppp_sub <- PPP
#ppp_sub$substrate <- data_position_substrat$substrate[grep(
#  pattern = stn,data_position_substrat$station)]
# (c) covariate values in data frame
substrate_function<-function(x,y){
  temporary <- data_substrate %>% filter(station==stn)
  substrate <- c()
  for (i in 1:length(y)){
    n=1
    if (y[i]>temporary$y_stop[n]){
      n <- n+1
    }
    substrate <- c(substrate,temporary$substrat[n])
  }
  return(as.factor(substrate))
}
Q <- quadscheme(ppp_sub)
xQ <- x.quad(Q)
yQ <- y.quad(Q)
Svalues <- substrate_function(xQ,yQ)
# poisson model
fit1 <- ppm(ppp_sub~ y)
fit_sub0 <- ppm(ppp_sub~Sub, data = data.frame(Sub=Svalues))
#anova( fit1,fit_sub0, test="LR")
fit_sub1 <- ppm(ppp_sub~y+Sub, data = data.frame(Sub=Svalues))
AIC(fit1);AIC(fit_sub0);AIC(fit_sub1)
goodfit <- step(ppm(ppp_sub ~ polynom(x,y,3)), trace=0)
formula(goodfit);AIC(goodfit)
formula <- as.formula(paste("ppp_sub",as.character(goodfit$trend)[1],
                            as.character(goodfit$trend)[2],
                            " + Sub",sep=""))
fit_sub2 <- ppm(formula, data = data.frame(Sub=Svalues))
formula <- as.formula(paste("ppp_sub",as.character(goodfit$trend)[1],
                            as.character(goodfit$trend)[2],
                            " + substrate_function",sep=""))
fit_sub2 <-  ppm(formula)
AIC(fit_sub2);AIC(goodfit)
summary(fit_sub2)


X_inhom <- simulate(goodfit)
X_sub <- simulate(fit_sub2)
Show_result <- rbind(
  as.data.frame(cbind(X_sub$`Simulation 1`$x,X_sub$`Simulation 1`$y,
                      rep("X_sub",X_sub$`Simulation 1`$n))),
  as.data.frame(cbind(X_inhom$`Simulation 1`$x,X_inhom$`Simulation 1`$y,
                      rep("X_inhom",X_inhom$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())


# Locally_weighted Poisson point process model
pp <- PPP
# 1) find an estimate of lambda Poisson
fit <- fit_sub2
# 2) calculate the local version of the K-function
localKfunction <- localK(pp, verbose = TRUE)
# 3) find an estimate of the interaction term
theo <- localKfunction$theo
nX <- npoints(pp)
r <- localKfunction$r
dr <- c(diff(r), diff(r)[nX - 1])
l1 <- vector(l = nX)
for(i in 1:nX){
  progressreport(i,nX)
  l1[i] <- exp(sum((sign(localKfunction[[i]] - theo) * (
    (localKfunction[[i]] - theo) ^ 2)) / theo * dr, 
    na.rm = TRUE))
}
pp1 <- pp
pp1$marks <- l1
# 4) smooth the surface with the Nadaraya-Watson smoother
im0_kernel <- Smooth(pp1, sigma = bw.diggle(pp1), positive = TRUE)
#bw.CvL=Cronie and van Lieshout's 
#Criterion for Bandwidth Selection for Kernel Density

#Log-transform the kernel (if all values are positive)
im0_kernel_log <- log(im0_kernel + 1e-6)
# Standardizing to avoid extreme values 
im0_kernel_log_scaled <- im0_kernel_log
im0_kernel_log_scaled[["v"]] <- scale(im0_kernel_log[["v"]])
summary(im0_kernel_log_scaled)

# 5) Fit the final model, maintaining the same covariate dependence specification
formula <- as.formula(paste("pp",as.character(fit$trend)[1],
                            as.character(fit$trend)[2],
                            " + offset(im0_kernel_log_scaled)",sep="")) 
lwppp <- ppm(formula)
AIC(lwppp)

X_inhom_sub <- simulate(fit)
X_lwppp <- simulate(lwppp)
Show_result <- rbind(
  as.data.frame(cbind(X_lwppp$`Simulation 1`$x,X_lwppp$`Simulation 1`$y,
                      rep("X_lwppp",X_lwppp$`Simulation 1`$n))),
  as.data.frame(cbind(X_inhom_sub$`Simulation 1`$x,X_inhom_sub$`Simulation 1`$y,
                      rep("X_inhom_sub",X_inhom_sub$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())



# Show tot #####
Show_result <- rbind(
as.data.frame(cbind(X_lwppp$`Simulation 1`$x,X_lwppp$`Simulation 1`$y,
                                  rep("X_lwppp",X_lwppp$`Simulation 1`$n))),
as.data.frame(cbind(X_inhom$`Simulation 1`$x,X_inhom$`Simulation 1`$y,
                                  rep("X_inhom",X_inhom$`Simulation 1`$n))),
as.data.frame(cbind(X_LGCP$`Simulation 1`$x,X_LGCP$`Simulation 1`$y,
                                  rep("X_LGCP",X_LGCP$`Simulation 1`$n))),
as.data.frame(cbind(X_sub$`Simulation 1`$x,X_sub$`Simulation 1`$y,
                    rep("X_sub",X_sub$`Simulation 1`$n))),
as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
