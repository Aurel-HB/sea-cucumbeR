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
library(nleqslv)

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

list_PPP_epsg4461 <- readRDS(paste(here(),
                    "/point_process/Data/list_PPP_2025_epsg4461.rds",sep=""))

########################
# Start with one PPP
########################

stn <- "96" # # stn="149" stn="179" stn="143"
PPP <- list_PPP[[stn]] # PPP<-list_PPP[["179"]] PPP<-list_PPP[["149"]]
# PPP<-list_PPP[["143"]]
PPP[["window"]][["units"]] <- list("metre","metres")
#PPP[["window"]][["yrange"]] <- PPP[["window"]][["yrange"]] - PPP[["window"]][["yrange"]][[1]]

# or PPP <- list_PPP_epsg4461[[stn]]

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
#Q–Q plot of smoothed raw residuals Baddeley 2005 ####
# 1. Fit model
fit <- goodfit
# 2. Smoothed residual field for observed data
r_obs <- Smooth(residuals(fit, type="raw"), sigma=0.05)
v_obs <- as.vector(r_obs$v)
v_obs <- v_obs[!is.na(v_obs)]
# Sort observed values
v_obs <- sort(v_obs)
# 3. Monte Carlo simulations
nsim <- 39
sim_patterns <- simulate(fit, nsim=nsim)
# 4. For each simulation:
sim_sorted <- matrix(NA, nrow=length(v_obs), ncol=nsim)
for(i in 1:nsim){
  # simulate pattern
  Xsim <- sim_patterns[[i]]
  # re-fit model (IMPORTANT step)
  fit_sim <- ppm(Xsim,formula(fit))
  # smoothed residual field
  r_sim <- Smooth(residuals(fit_sim, type="raw"), sigma=0.05)
  
  v_sim <- as.vector(r_sim$v)
  v_sim <- v_sim[!is.na(v_sim)]
  # sort values (order statistics)
  sim_sorted[,i] <- sort(v_sim)
}
# 5. Expected quantiles (mean of order statistics)
q_mean <- rowMeans(sim_sorted)
# Envelopes (pointwise)
q_lo <- apply(sim_sorted, 1, quantile, 0.025)
q_hi <- apply(sim_sorted, 1, quantile, 0.975)
# 6. Plot Q–Q
plot(q_mean, v_obs,
     xlab="Mean quantile of simulations",
     ylab="Data quantile",
     main="Smoothed residuals: raw")
abline(0,1)
lines(q_mean, q_lo, lty=2, col="red")
lines(q_mean, q_hi, lty=2, col="red")

#Envelope Tests (K-Function) ####
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
legend("top", legend = c("Observed", "Simulated"), lty = c(1, 2),
       col = c("black", "gray"))




# Cox model ####
#fitox0 <- kppm(PPP~1, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
#fitox0 <- kppm(PPP~1, clusters = "LGCP", method="waag", model="matern",nu=0.3)
#fitox0 <- kppm(PPP~1, clusters = "LGCP", method="mincon", model="matern",nu=0.3)
#fitox0 <- kppm(PPP~1, clusters = "LGCP", statistic="pcf", model="matern",nu=0.3)
formule <- formula(goodfit)
#fitox1 <- kppm(PPP,formule, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
fitox1 <- kppm(PPP~y, clusters = "LGCP", method="clik2", model="matern",nu=0.3)
#fitox1 <- kppm(PPP,formule, clusters = "LGCP", method="waag", model="matern",nu=0.3)

#simulation
X_LGCP <- simulate(fitox1); plot(X_LGCP[[1]])
Show_result <- rbind(
  as.data.frame(cbind(X_LGCP$`Simulation 1`$x,X_LGCP$`Simulation 1`$y,
                      rep("X_LGCP",X_LGCP$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())

#Q–Q plot of smoothed raw residuals Baddeley 2005 ####
# 1. Fit model
fit <- fitox1
# 2. Smoothed residual field for observed data
r_obs <- Smooth(residuals(fit, type="raw"), sigma=0.05)
v_obs <- as.vector(r_obs$v)
v_obs <- v_obs[!is.na(v_obs)]
# Sort observed values
v_obs <- sort(v_obs)
# 3. Monte Carlo simulations
nsim <- 39
sim_patterns <- simulate(fit, nsim=nsim)
# 4. For each simulation:
sim_sorted <- matrix(NA, nrow=length(v_obs), ncol=nsim)
for(i in 1:nsim){
  # simulate pattern
  Xsim <- sim_patterns[[i]]
  # re-fit model (IMPORTANT step)
  fit_sim <- kppm(Xsim~y, clusters = "LGCP", method="clik2",
                  model="matern",nu=0.3)
  # smoothed residual field
  r_sim <- Smooth(residuals(fit_sim, type="raw"), sigma=0.05)
  
  v_sim <- as.vector(r_sim$v)
  v_sim <- v_sim[!is.na(v_sim)]
  # sort values (order statistics)
  sim_sorted[,i] <- sort(v_sim)
}
# 5. Expected quantiles (mean of order statistics)
q_mean <- rowMeans(sim_sorted)
# Envelopes (pointwise)
q_lo <- apply(sim_sorted, 1, quantile, 0.025)
q_hi <- apply(sim_sorted, 1, quantile, 0.975)
# 6. Plot Q–Q
plot(q_mean, v_obs,
     xlab="Mean quantile of simulations",
     ylab="Data quantile",
     main="Smoothed residuals: raw")
abline(0,1)
lines(q_mean, q_lo, lty=2, col="red")
lines(q_mean, q_hi, lty=2, col="red")

# extract point process to use it on a another surface ####
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
#####

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
im0_kernel <- Smooth(pp1, sigma = bw.ppl(pp1), positive = TRUE)
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

#Q–Q plot of smoothed raw residuals Baddeley 2005 ####
# 1. Fit model
fit <- lwppp
# 2. Smoothed residual field for observed data
r_obs <- Smooth(residuals(fit, type="raw"), sigma=0.05)
v_obs <- as.vector(r_obs$v)
v_obs <- v_obs[!is.na(v_obs)]
# Sort observed values
v_obs <- sort(v_obs)
# 3. Monte Carlo simulations
nsim <- 39
sim_patterns <- simulate(fit, nsim=nsim)
# 4. For each simulation:
sim_sorted <- matrix(NA, nrow=length(v_obs), ncol=nsim)
for(i in 1:nsim){
  # simulate pattern
  Xsim <- sim_patterns[[i]]
  # re-fit model (IMPORTANT step)
  formula <- as.formula(paste("Xsim",as.character(fit$trend)[1],
                              as.character(fit$trend)[2],
                              " + offset(im0_kernel_log_scaled)",sep="")) 
  fit_sim <- ppm(formula)
  # smoothed residual field
  r_sim <- Smooth(residuals(fit_sim, type="raw"), sigma=0.05)
  
  v_sim <- as.vector(r_sim$v)
  v_sim <- v_sim[!is.na(v_sim)]
  # sort values (order statistics)
  sim_sorted[,i] <- sort(v_sim)
}
# 5. Expected quantiles (mean of order statistics)
q_mean <- rowMeans(sim_sorted)
# Envelopes (pointwise)
q_lo <- apply(sim_sorted, 1, quantile, 0.025)
q_hi <- apply(sim_sorted, 1, quantile, 0.975)
# 6. Plot Q–Q
plot(q_mean, v_obs,
     xlab="Mean quantile of simulations",
     ylab="Data quantile",
     main="Smoothed residuals: raw")
abline(0,1)
lines(q_mean, q_lo, lty=2, col="red")
lines(q_mean, q_hi, lty=2, col="red")

# Validation by scoring rules ####










# Show tot #####
Show_result <- rbind(
  as.data.frame(cbind(X_lwppp$`Simulation 1`$x,X_lwppp$`Simulation 1`$y,
                      rep("X_lwppp",X_lwppp$`Simulation 1`$n))),
  as.data.frame(cbind(X_inhom$`Simulation 1`$x,X_inhom$`Simulation 1`$y,
                      rep("X_inhom",X_inhom$`Simulation 1`$n))),
  as.data.frame(cbind(X_LGCP$`Simulation 1`$x,X_LGCP$`Simulation 1`$y,
                      rep("X_LGCP",X_LGCP$`Simulation 1`$n))),
  as.data.frame(cbind(PPP$x,PPP$y,rep("Initial Point Pattern",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme()


#####

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
as.data.frame(cbind(PPP$x,PPP$y,rep("Initial Point Pattern",PPP$n))))
Show_result$V1 <- as.numeric(Show_result$V1)
Show_result$V2 <- as.numeric(Show_result$V2)

ggplot(data = Show_result)+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme()
