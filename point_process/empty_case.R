# This code is a the version of point process analysis when the locally weighted
# Poisson point process can't be perform because there is no interaction term

## Check the resiudal ####
##Q–Q plot of smoothed raw residuals Baddeley 2005
residual <- data.frame(
  stn = rep(stn,3),
  model = c("fit_ihP","fit_LGCP","fit_lwppp"),
  cor.coef = NA,
  lm.coef = NA,
  p_value_envelope = NA
)

for (indic in 1:2){
  # 1. Fit model
  fit <- list(fit_ihP,fit_LGCP,fit_lwppp)[[indic]]
  name <- residual$model[indic]
  
  # 2. Smoothed residual field for observed data
  r_obs <- Smooth(residuals(fit, type="raw"), sigma=0.05)
  v_obs <- as.vector(r_obs$v)
  v_obs <- v_obs[!is.na(v_obs)]
  # Sort observed values
  v_obs <- sort(v_obs)
  
  # 3. Monte Carlo simulations
  nsim <- 39
  set.seed(1000 + indic)  # <-- seed specific to each model, for independent reproducibility
  sim_patterns <- simulate(fit, nsim=nsim)
  
  # 4. For each simulation:
  sim_sorted <- matrix(NA, nrow=length(v_obs), ncol=nsim)
  for(i in 1:nsim){
    # simulate pattern
    Xsim <- sim_patterns[[i]]
    # re-fit model (IMPORTANT step)
    formula <- formula(fit)
    fit_sim <- NULL  # <-- reset on each iteration (to fix a potential bug)
    if (name == "fit_LGCP"){
      # be careful to simulation with very few point
      if (Xsim$n > 4){
        fit_sim <- kppm(Xsim~y, clusters = "LGCP", #method="clik2",
                        model="matern", nu=0.3)
      }
    } else {
      fit_sim <- ppm(Xsim,formula) 
    }
    
    if (is.null(fit_sim)){
      # Simulation with too few data points: this column is ignored (remains NA)
      next
    }
    
    # smoothed residual field
    r_sim <- Smooth(residuals(fit_sim, type="raw"), sigma=0.05)
    v_sim <- as.vector(r_sim$v)
    v_sim <- v_sim[!is.na(v_sim)]
    # sort values (order statistics)
    sim_sorted[,i] <- sort(v_sim)
  }
  
  # We remove the columns that have not been calculated (LGCP cases with too few data points)
  valid_cols <- colSums(is.na(sim_sorted)) < nrow(sim_sorted)
  sim_sorted <- sim_sorted[, valid_cols, drop = FALSE]
  nsim_valid <- ncol(sim_sorted)
  
  # 5. Expected quantiles (mean of order statistics)
  q_mean <- rowMeans(sim_sorted)
  # Envelopes (pointwise)
  q_lo <- apply(sim_sorted, 1, quantile, 0.025)
  q_hi <- apply(sim_sorted, 1, quantile, 0.975)
  
  # 6. Extract correlation from plot Q–Q
  residual$cor.coef[indic] <- as.numeric(cor.test(q_mean,v_obs)$estimate)
  residual$lm.coef[indic] <- as.numeric(
    lm(q_mean~v_obs)$coefficients["v_obs"])
  
  ## --- Envelope Indicator ---
  # Option 4: Global envelope test (MAD type, Myllymäki/Baddeley)
  #We combine the observed data (column 1) and the simulations (columns 2:n) into a single matrix
  all_curves <- cbind(v_obs, sim_sorted)
  
  # Point wise standard deviation calculated based solely on simulations (test benchmark)
  sd_row <- apply(sim_sorted, 1, sd)
  sd_row[sd_row == 0] <- NA  # évite division par 0
  
  # MAD statistic: standardized maximum absolute deviation, for each curve (observed + simulated)
  mad_stat <- apply(all_curves, 2, function(v) max(abs(v - q_mean) / sd_row, na.rm = TRUE))
  
  # Rank of the observed curve among the (1 + nsim_valid) curves
  # rank 1 = the most extreme (least consistent with the model)
  rank_obs <- rank(-mad_stat, ties.method = "average")[1]
  
  # p-value Monte-Carlo (resolution = 1/(nsim_valid+1))
  residual$p_value_envelope[indic] <- rank_obs / (nsim_valid + 1)
  
  # 7. Q-Q plot with envelope
  smoothed_residuals_raw[[paste(
    c("fit_ihP","fit_LGCP","fit_lwppp")[[indic]],
    stn, sep="_"
  )]] <- ggplot(data.frame(
    q_mean = q_mean,
    v_obs = v_obs,
    q_lo = q_lo,
    q_hi = q_hi), 
    aes(x = q_mean, y = v_obs)) +
    geom_point() + 
    geom_abline(slope = 1, intercept = 0, color = "black") +
    geom_line(aes(y = q_lo), linetype = "dashed", color = "red") +
    geom_line(aes(y = q_hi), linetype = "dashed", color = "red") +
    labs(
      x = "Mean quantile of simulations",
      y = "Data quantile",
      title = paste0("Smoothed residuals: raw (p = ",
                     round(residual$p_value_envelope[indic], 3), ")")
    ) +
    theme_minimal()
}

residual_tot <- rbind(residual_tot,residual)
count_Xsim_LGCP <- rbind(count_Xsim_LGCP,(data.frame(stn, nsim_valid)))

## Apply the scorin of rules of the intensity and K'Ripley function ####
### Setup scoring rule functions ####
#source(paste(here(),"/point_process/setup_scoring_rules.R",sep=""))
N = 100 # samples to compute for CRPS

##### K-score ####
eval_points = seq(0.1,0.375,by = 0.005) # points where the K estimator is to be evaluated. 
#Needs to take very specific form due to the way K_hat is coded in spatstat

K_hat = function(dat){return(as.function(Kinhom(dat,
                                                lambda = density.ppp(
                                                  PPP,bw.scott(PPP)),
                                                nlarge = 2500))
                             (eval_points))}

NEst_Kscore = computeNEst(mod = models[1:2],est = K_hat,N = N)

#K_scores = get_crps(dt = NEst_Kscore,models = models[1:2])

K_PPP <- as.function(
    Kinhom(PPP,lambda = density.ppp(PPP,bw.scott(PPP)))
)(eval_points)
K_scores_PPP <- get_crps_PPP(dt = NEst_Kscore,models = models[1:2],est_PPP = K_PPP)

##### Intensity score ####

lambda_hat = function(dat){return(as.matrix(density(dat)))}

NEst_lambdascore = computeNEst(mod = models[1:2],est = lambda_hat,N = N,
                               silence = FALSE)

#lambda_scores = get_crps(dt = NEst_lambdascore,models = models[1:2])

lambda_PPP <- as.matrix(density(PPP))
lambda_scores_PPP <- get_crps_PPP(dt = NEst_lambdascore,
                                  models = models[1:2],est_PPP = lambda_PPP)

scoring_tot <- rbind(scoring_tot,
                     data.frame(
                       STN = stn,
                       model = models,
                       lambda_score = as.numeric(c(lambda_scores_PPP[1,1],
                                                   lambda_scores_PPP[1,2],
                                                   NA)),
                       K_score = as.numeric(c(K_scores_PPP[1,1],
                                              K_scores_PPP[1,2],
                                              NA))
                     ))
