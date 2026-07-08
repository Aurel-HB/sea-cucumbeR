# this code realised the adjustement of 3 types point process models on the 61
# point patterns from the 2025 survey

### ### ### ### ###
#Load packages ####
### ### ### ### ###
library(here)
library(sf)
library(spatstat)
library(spatstat.model)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(fields)

### ### ### ###
#Load data ####
### ### ### ###
list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

check_warnings <- function(expr) {
  warnings_list <- character(0)
  result <- withCallingHandlers(
    expr,
    warning = function(w) {
      warnings_list <<- c(warnings_list, as.character(w))
      invokeRestart("muffleWarning")
    }
  )
  list(result = result, warnings = warnings_list)
}

### ### ### ###
#Analysis  ####
### ### ### ###
## Definition of the variable ####
scoring_tot <- data.frame()
residual_tot <- data.frame()
warnings_list <- list()
for (index in 1:length(list_PPP)){
  stn <- names(list_PPP)[index]
  PPP <- list_PPP[[stn]]
  PPP[["window"]][["units"]] <- list("metre","metres")
## Fit the 3 models ####
### Inhomogeneous Poisson Point Process Model ####
  #select best model for intensity function
  fit_ihP <- step(ppm(PPP ~ polynom(x,y,2)), trace=0)
  # show the potential warnings
  print(paste(stn,"Inhomogeneous Poisson Point Process Model", sep="_"))
  print(fit_ihP)
  
### Inhomogeneous Log-Gaussian Cox Point Process Model ####
  fit_LGCP <- kppm(PPP~y, clusters = "LGCP",
                 method="clik2", model="matern",nu=0.3)
  # show the potential warnings
  print(paste(stn,"Inhomogeneous Log-Gaussian Cox Point Process Model",
              sep="_"))
  print(fit_LGCP)
### Locally Weighted Poisson Point Process Model ####
  # 1) find an estimate of lambda Poisson using the previous Poisson model
  fit <- fit_ihP
  # 2) calculate the local version of the K-function
  localKfunction <- localK(PPP, verbose = TRUE)
  # 3) find an estimate of the interaction term
  theo <- localKfunction$theo
  nX <- npoints(PPP)
  r <- localKfunction$r
  dr <- c(diff(r), diff(r)[nX - 1])
  l1 <- vector(l = nX)
  for(i in 1:nX){
    progressreport(i,nX)
    term  <- sign(localKfunction[[i]] - theo) * 
      (localKfunction[[i]] - theo)^2 / theo * dr
    l1[i] <- exp(pmin(sum(term, na.rm = TRUE),700)) # avoid inf value
  }
  pp1 <- PPP
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
  formula <- as.formula(paste("PPP",as.character(fit$trend)[1],
                              as.character(fit$trend)[2],
                              " + offset(im0_kernel_log_scaled)",sep=""))
  fit_lwppp <- ppm(formula)
  
  # show the potential warnings
  print(paste(stn,"Locally-weighted Inhomogeneous Poisson Point Process Model",
              sep="_"))
  print(fit_lwppp)
  
## Check the resiudal ####
  ##Q–Q plot of smoothed raw residuals Baddeley 2005
  residual <- data.frame(
    stn = rep(stn,3),
    model = c("fit_ihP","fit_LGCP","fit_lwppp"),
    cor.coef = NA,
    lm.coef = NA  
  )
  for (index in 1:3){
    # 1. Fit model
    fit <- list(fit_ihP,fit_LGCP,fit_lwppp)[[index]]
    name <- residual$model[index]
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
      formula <- formula(fit)
      if (name == "fit_LGCP"){
        fit_sim <- kppm(Xsim~y, clusters = "LGCP", method="clik2",
                        model="matern",nu=0.3)
      } else {
        fit_sim <- ppm(Xsim,formula) 
      }
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
    # 6. Extract correlation from plot Q–Q
    residual$cor.coef[index] <- as.numeric(cor.test(q_mean,v_obs)$estimate)
    residual$lm.coef[index] <- as.numeric(
      lm(q_mean~v_obs)$coefficients["v_obs"])
  }
## Apply the scorin of rules of the intensity and K'Ripley function ####
  ### Setup scoring rule functions ####
  source(paste(here(),"/point_process/setup_scoring_rules.R",sep=""))
 
}

# Show Point Pattern ####
ggplot(data = as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))+
  geom_point(aes(x=V1,y=V2,color=V3))+
  facet_wrap(~V3,nrow=1)+
  ylab("")+xlab("")+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank())
