# this script compile all the work tested to compare the Gamma and Lognormal
# distribution in the context of sea cucumber biomass densities distribution
# model. This code produce the analysis and result for the paper associated to
# this study.

# load packages ####
library(here)
library(sf)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sdmTMB)
library(INLA)
library(RColorBrewer)
library(viridis)
library(ape)
library(ggpubr)
library(cowplot)
library(future)

# load data ####
calcul_area <- readRDS(paste(here(),
                             "/SIG/SIG_Data/study_calcul_area.rds",sep=""))

data_holotv <- readRDS(
  paste(here(),"/process_spatiotemp_data/Data/processed/data_abun_tot_cov.rds",
        sep = "")
)

grid_bathy <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_bathy.rds",
        sep = "")
)

grid_tmp <- readRDS(
  paste(here(),
      "/environment_exploration/Environment_Data/processed/grid_bottom_tmp.rds",
      sep = "")
)

grid_chl <- readRDS(
  paste(here(),
        "/environment_exploration/Environment_Data/processed/grid_chl.rds",
        sep = "")
)

# prepare data ####
## Prepare the grid for prediction ####
# The survey is only in may so we keep the may data
grid_tmp <- grid_tmp[,grep(pattern = "May",names(grid_tmp))]
grid_chl <- grid_chl[,grep(pattern = "May",names(grid_chl))]

# The grid have to be replicated for each year with the environmental covariate values
grid_tot <- data.frame()
for (i in 1:5){
  data <- st_join(grid_tmp[,i],grid_bathy)
  data <- st_join(data,grid_chl[,i])
  data <- data %>% mutate(year = 2020+i)
  names(data) <- c("bottomT","long","lat","bathy","chla","geometry","year")
  grid_tot <- rbind(grid_tot,data[,c(1,2,3,4,5,7,6)])
} # the lat and long are in km for modelling
rm(data)
## Prepare the date for modelling ####
# synchronize variable name
names(grid_tot)[2:3] <- c("X","Y")
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
grid_proj <- as.data.frame(grid_tot)[,1:6]

# extract the surface of the study area
stock_surface <- as.numeric(st_area(calcul_area)/1e6)


###-###-###-###-###-###-###-###
# Mesh sensibility analysis ####
###-###-###-###-###-###-###-###

## creation of the mesh with R-INLA ####
cutoff <- c(1,1.852,3,5)
list_mesh_regular <- list()
list_mesh_data_fit <- list()
for (cut in cutoff){
  #using the grid of the area
  bnd <- INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                   convex = -0.04)
  bnd2 = INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                   convex = -0.15)
  
  mesh_inla <- INLA::inla.mesh.2d(
    loc = as.matrix(as.data.frame(grid_bathy)[,c(1,2)]),
    boundary = list(bnd,bnd2),
    cutoff = cut, # minimum triangle edge length
    max.edge = c(3*cut, 20*cut), # inner and outer max triangle lengths
  ) # 1.852 is the distance in meter of a mile nautic
  
  list_mesh_regular[[paste("mesh_regular_",cut,sep="")]] <- make_mesh(
    data_estimate,c("X", "Y"), mesh = mesh_inla)
  
  mesh_inla <- INLA::inla.mesh.2d(
    loc = as.matrix(data_estimate[,c("X","Y")]),
    boundary = list(bnd,bnd2),
    cutoff = cut, # minimum triangle edge length
    max.edge = c(3*cut, 20*cut), # inner and outer max triangle lengths
  ) # 1.852 is the distance in meter of a mile nautic
  
  list_mesh_data_fit[[paste("mesh_data_fit_",cut,sep="")]] <- make_mesh(
    data_estimate,c("X", "Y"), mesh = mesh_inla)
}
rm(mesh_data_fit,mesh_inla,mesh_regular)

for (i in 1:length(cutoff)){
  plot(list_mesh_regular[[i]]$mesh, main = NA, edge.color = "grey60", asp = 1)
  points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)
  
  plot(list_mesh_data_fit[[i]]$mesh, main = NA, edge.color = "grey60", asp = 1)
  points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)
}


## fit model for the 2 meshes' type ####
data_estimate$time <- data_estimate$year
grid_proj$time <- grid_proj$year

list_fit <- list()
for(i in 1:length(cutoff)){
  cut <- cutoff[i]
  mesh_regular <- list_mesh_regular[[i]]
  mesh_data_fit <- list_mesh_data_fit[[i]]
  
  fit_gamma_regular <- sdmTMB(
    density.t_km2 ~ 1+as.factor(year),
    data = data_estimate,
    mesh = mesh_regular,
    family = Gamma(link = "log"),
    spatial = "on",
    time = "time",
    spatiotemporal = "IID"
  )
  
  fit_gamma_data_fit <- sdmTMB(
    density.t_km2 ~ 1+as.factor(year),
    data = data_estimate,
    mesh = mesh_data_fit,
    family = Gamma(link = "log"),
    spatial = "on",
    time = "time",
    spatiotemporal = "IID"
  )
  
  fit_logn_regular <- sdmTMB(
    density.t_km2 ~ 1+as.factor(year),
    data = data_estimate,
    mesh = mesh_regular,
    family = lognormal(link = "log"),
    spatial = "on",
    time = "time",
    spatiotemporal = "IID"
  )
  
  fit_logn_data_fit <- sdmTMB(
    density.t_km2 ~ 1+as.factor(year),
    data = data_estimate,
    mesh = mesh_data_fit,
    family = lognormal(link = "log"),
    spatial = "on",
    time = "time",
    spatiotemporal = "IID"
  )
  
  list_fit[[paste("fit_gamma_regular_",cut,sep="")]] <- fit_gamma_regular
  list_fit[[paste("fit_gamma_data_fit_",cut,sep="")]] <- fit_gamma_data_fit
  list_fit[[paste("fit_logn_regular_",cut,sep="")]] <- fit_logn_regular
  list_fit[[paste("fit_logn_data_fit_",cut,sep="")]] <- fit_logn_data_fit
}
rm(fit_logn_data_fit,fit_logn_regular,fit_gamma_data_fit,fit_gamma_regular)

## check the result ####
for (i in 1:length(list_fit)){
  sanity(list_fit[[i]])
  print(paste("Check number",i,sep = " "))
}
# model number 12 and 15 have not pass the sanity check
#list_fit[[12]] <- "sanity_error"
#list_fit[[15]] <- "sanity_error"

for (i in 1:length(cutoff)){
  par(mfrow = c(2,2))
  resids_gamma_regular <- residuals(list_fit[[1+4*(i-1)]]) 
  # randomized quantile residuals
  qqnorm(resids_gamma_regular, main = "Normal Q-Q Plot Gamma Regular Mesh")
  qqline(resids_gamma_regular)
  
  resids_gamma_data_fit <- residuals(list_fit[[2+4*(i-1)]]) 
  # randomized quantile residuals
  qqnorm(resids_gamma_data_fit, main = "Normal Q-Q Plot Gamma Data-fit Mesh")
  qqline(resids_gamma_data_fit)
  
  resids_logn_regular <- residuals(list_fit[[3+4*(i-1)]]) 
  # randomized quantile residuals
  qqnorm(resids_logn_regular, main = "Normal Q-Q Plot Lognormal Regular Mesh")
  qqline(resids_logn_regular)
  
  resids_logn_data_fit <- residuals(list_fit[[4+4*(i-1)]]) 
  # randomized quantile residuals
  qqnorm(resids_logn_data_fit, main = "Normal Q-Q Plot Lognormal Data-fit Mesh")
  qqline(resids_logn_data_fit)
  par(mfrow = c(1,1))
}


### Moran index per year ###
moran <- data.frame()
sim <- c()
for(c in 1:length(cutoff)){
  for (i in 1:1000){
    data_estimate$resids_gamma_regular <- residuals(list_fit[[1+4*(c-1)]]) 
    data_estimate$resids_gamma_data_fit <- residuals(list_fit[[2+4*(c-1)]])
    data_estimate$resids_logn_regular <- residuals(list_fit[[3+4*(c-1)]]) 
    data_estimate$resids_logn_data_fit <- residuals(list_fit[[4+4*(c-1)]])
    
    for (time in unique(data_estimate$year)){
      # calcul weight
      dists <- as.matrix(dist(data_estimate[
        grep(time,data_estimate$year), c("X", "Y")]))
      inv_dists <- 1 / dists
      diag(inv_dists) <- 0
      inv_dists[is.infinite(inv_dists)] <- 0
      
      #gamma regular
      data <- Moran.I(data_estimate$resids_gamma_regular[data_estimate$year==time],
                      inv_dists, scaled = TRUE)
      data$year <- time
      data$mesh <- "Regular"
      data$model <- "Gamma"
      data$cutoff <- cutoff[c]
      moran <- rbind(moran, data)
      
      #gamma data-fit
      data <- Moran.I(data_estimate$resids_gamma_data_fit[data_estimate$year==time],
                      inv_dists, scaled = TRUE)
      data$year <- time
      data$mesh <- "Data-fit"
      data$model <- "Gamma"
      data$cutoff <- cutoff[c]
      moran <- rbind(moran, data)
      
      #logn regular
      data <- Moran.I(data_estimate$resids_logn_regular[data_estimate$year==time],
                      inv_dists, scaled = TRUE)
      data$year <- time
      data$mesh <- "Regular"
      data$model <- "Lognormal"
      data$cutoff <- cutoff[c]
      moran <- rbind(moran, data)
      
      #logn data-fit
      data <- Moran.I(data_estimate$resids_logn_data_fit[data_estimate$year==time],
                      inv_dists, scaled = TRUE)
      data$year <- time
      data$mesh <- "Data-fit"
      data$model <- "Lognormal"
      data$cutoff <- cutoff[c]
      moran <- rbind(moran, data)
    }
    sim <- c(sim, rep(paste("sim",i,sep="_"),16))
  }
  
}

moran$sim <- sim
moran$signif <- FALSE
moran$signif[moran$p.value<0.05] <- TRUE
#saveRDS(moran,paste(here(),
#                    "/comparing_biomass_method/output/moran_simulation.rds",
#                    sep=""))
#moran <- readRDS(paste(here(),
#                    "/comparing_biomass_method/output/moran_simulation.rds",
#                    sep=""))

summary_moran <- data.frame(moran[grep("sim_1000",moran$sim),
                                  c("year","mesh","model","cutoff")],
                            percent_ns=0)
for (i in 1:nrow(summary_moran)){
  time = summary_moran$year[i]
  type = summary_moran$mesh[i]
  distrib = summary_moran$model[i]
  cut = summary_moran$cutoff[i]
  data <- moran %>% filter(year==time) %>% filter(mesh==type) %>%
    filter(model==distrib) %>% filter(cut==cutoff)
  summary_moran$percent_ns[i] <- length(grep(TRUE,data$signif))*100/nrow(data)
}
rm(time,type,data,distrib,cut)

ggplot(data = summary_moran, aes(x = as.factor(year), y = as.factor(mesh),
                                 fill = percent_ns))+
  geom_tile(color = "white")+
  scale_fill_distiller(palette = "RdYlBu")+
  geom_text(aes(label = percent_ns))+
  facet_wrap(cutoff~model, ncol=2, labeller = labeller(cutoff = label_both, 
                                                       model = label_value))+
  labs(fill = "% significant Moran's index", x = "Year", y = "Mesh type")

# It seems that the gamma models capt better the residual spatial autocorrelation



### Observed vs predicted ###
list_ggscatter <- list()
R_vect <- c()
for (i in 1:16){
  predicted <- exp(predict(list_fit[[i]])$est)
  observed <- list_fit[[i]]$data$density.t_km2
  
  list_ggscatter[[names(list_fit)[i]]] <- ggscatter(
    data.frame(observed,predicted),x = "observed", y = "predicted",
            add = "reg.line", title = names(list_fit)[i]) +
    stat_cor(label.x = 1000, label.y = 3000) +
    stat_regline_equation(label.x = 2000, label.y = 6.5)
  
  R_vect <- c(R_vect, as.numeric(cor.test(predicted,observed)$estimate))
}

plot_grid(
  list_ggscatter[[1]],
  list_ggscatter[[2]],
  list_ggscatter[[3]],
  list_ggscatter[[4]],
  list_ggscatter[[5]],
  list_ggscatter[[6]],
  list_ggscatter[[7]],
  list_ggscatter[[8]],
  list_ggscatter[[9]],
  list_ggscatter[[10]],
  list_ggscatter[[11]],
  list_ggscatter[[12]],
  list_ggscatter[[13]],
  list_ggscatter[[14]],
  list_ggscatter[[15]],
  list_ggscatter[[16]]
)

### big ggplot that summarize Moran and R² results ###
AFH_theme <-  theme(panel.background = element_rect(fill = '#dddddd', 
                                                    color = '#dddddd',
                                                    linewidth = 1),
                     panel.grid.major = element_line(color = '#00000000', 
                                                     linetype = 'dotted'),
                     panel.grid.minor = element_line(color = '#00000000', 
                                                     linewidth = 2),
                     plot.background = element_rect(fill = "#b4c7dc"),
                     legend.background = element_rect(fill = "#b4c7dc"),
                    strip.background = element_rect(fill ="#dddddd"))

# add a column configuration with mesh type x cutoff
summary_moran <- summary_moran %>%
  mutate(config = paste(mesh,"_cutoff: ",cutoff,sep="")) %>%
  mutate(year = as.factor(year)) %>%
  mutate(config = as.factor(config))

ggplot(summary_moran, aes(year, config, fill = percent_ns)) +
  geom_tile(color = "#00000000")+
  geom_text(aes(label = percent_ns))+
  scale_fill_distiller(palette = "RdYlBu")+
  facet_wrap(~ model, nrow = 1) +
  AFH_theme+
  labs(y = "Mesh configurations",
       x = "Year",
       fill = "% significant Moran's index")

# data.frame for the coefficient correlation
R.coef_data <- data.frame(
  model = names(list_fit),
  R.coef = R_vect
)
R.coef_data <- R.coef_data %>%
  mutate(distrib = ifelse(grepl("gamma", model),"Gamma","Lognormal")) %>%
  mutate(mesh = ifelse(grepl("regular", model),"Regular","Data-fit")) %>%
  mutate(cutoff = as.factor(c(1,1,1,1,1.852,1.852,1.852,1.852,3,3,3,3,5,5,5,5)))
R.coef_data$cutoff <- ordered(R.coef_data$cutoff,c(5,3,1.852,1))
R.coef_data$mesh <- ordered(R.coef_data$mesh,c("Regular","Data-fit"))

ggplot(R.coef_data %>% filter(distrib=="Gamma"), 
       aes(mesh, cutoff, fill = R.coef)) +
  geom_tile(color = "#00000000")+
  geom_text(aes(label = round(R.coef,digits = 2)*100))+
  scale_fill_distiller(palette = "RdYlBu", direction = 1)+
  facet_wrap(~distrib, ncol=4, labeller = labeller(distrib = label_value))+
  AFH_theme+
  labs(y = "Mesh cutoff",
       x = "Mesh type")

## cross validation ####
# the previous result enable to choose the gamma models with 1, 1.852, 3 cutoff

### cross validation random ###
n_folds = 20

set.seed(44)
plan(multisession, workers = 8)
m_cv_gamma_regular_1 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_1,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_gamma_regular_1$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(44)
plan(multisession, workers = 8)
m_cv_gamma_data_fit_1 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_data_fit$mesh_data_fit_1,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_gamma_data_fit_1$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(44)
plan(multisession, workers = 8)
m_cv_gamma_regular_1.852 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_1.852,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_gamma_regular_1.852$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(44)
plan(multisession, workers = 8)
m_cv_gamma_data_fit_1.852 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_data_fit$mesh_data_fit_1.852,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_gamma_data_fit_1.852$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(44)
plan(multisession, workers = 8)
m_cv_gamma_regular_3 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_3,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_gamma_regular_3$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

m_cv_gamma_regular_1$sum_loglik # total log-likelihood
m_cv_gamma_data_fit_1$sum_loglik # total log-likelihood
m_cv_gamma_regular_1.852$sum_loglik # total log-likelihood
m_cv_gamma_data_fit_1.852$sum_loglik # total log-likelihood
m_cv_gamma_regular_3$sum_loglik # total log-likelihood

### cross validation with spatial cluster ###

#cluster by time and space 
k <- 8
clust_2021 <- kmeans(data_estimate %>% filter(year==2021)%>%select(X,Y), k)$cluster
clust_2022 <- kmeans(data_estimate %>% filter(year==2022)%>%select(X,Y), k)$cluster
clust_2023 <- kmeans(data_estimate %>% filter(year==2023)%>%select(X,Y), k)$cluster
clust_2025 <- kmeans(data_estimate %>% filter(year==2025)%>%select(X,Y), k)$cluster
clust <- c(clust_2021,clust_2022+k,clust_2023+(k*2),clust_2025+(k*3))
data_estimate$clust <- clust


set.seed(13)
plan(multisession, workers = 8)
s_cv_gamma_regular_1 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_1,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_gamma_regular_1$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(13)
plan(multisession, workers = 8)
s_cv_gamma_data_fit_1 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_data_fit$mesh_data_fit_1,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_gamma_data_fit_1$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(13)
plan(multisession, workers = 8)
s_cv_gamma_regular_1.852 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_1.852,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_gamma_regular_1.852$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(13)
plan(multisession, workers = 8)
s_cv_gamma_data_fit_1.852 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_data_fit$mesh_data_fit_1.852,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_gamma_data_fit_1.852$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(13)
plan(multisession, workers = 8)
s_cv_gamma_regular_3 <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = list_mesh_regular$mesh_regular_3,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_gamma_regular_3$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

s_cv_gamma_regular_1$sum_loglik # total log-likelihood
s_cv_gamma_data_fit_1$sum_loglik # total log-likelihood
s_cv_gamma_regular_1.852$sum_loglik # total log-likelihood
s_cv_gamma_data_fit_1.852$sum_loglik # total log-likelihood
s_cv_gamma_regular_3$sum_loglik # total log-likelihood

# all the model from the CV had converge
# we stack the sum_loglik in a dataframe
cv_sum_loglik <- data.frame(
  model = rep(c("gamma_regular_1","gamma_data_fit_1",
                "gamma_regular_1.852","gamma_data_fit_1.852",
                "gamma_regular_3"),2),
  CV_type = c(rep("Random",5),rep("Spatial_cluster",5)),
  sum_loglik = c(
    m_cv_gamma_regular_1$sum_loglik,
    m_cv_gamma_data_fit_1$sum_loglik, 
    m_cv_gamma_regular_1.852$sum_loglik, 
    m_cv_gamma_data_fit_1.852$sum_loglik,
    m_cv_gamma_regular_3$sum_loglik,
    s_cv_gamma_regular_1$sum_loglik, 
    s_cv_gamma_data_fit_1$sum_loglik, 
    s_cv_gamma_regular_1.852$sum_loglik, 
    s_cv_gamma_data_fit_1.852$sum_loglik,
    s_cv_gamma_regular_3$sum_loglik
  )
)
cv_sum_loglik <- cv_sum_loglik %>%
  mutate(distrib = ifelse(grepl("gamma", model),"Gamma","Lognormal")) %>%
  mutate(mesh = ifelse(grepl("regular", model),"Regular","Data-fit")) %>%
  mutate(cutoff = as.factor(rep(c(1,1,1.852,1.852,3),2)))
cv_sum_loglik$cutoff <- ordered(cv_sum_loglik$cutoff,c(3,1.852,1))
cv_sum_loglik$mesh <- ordered(cv_sum_loglik$mesh,c("Regular","Data-fit"))


ggplot(cv_sum_loglik %>% filter(CV_type=="Random"),
       aes(mesh, cutoff, fill = sum_loglik)) +
  geom_tile(color = "#00000000")+
  geom_text(aes(label = round(sum_loglik,digits = 0)))+
  scale_fill_distiller(palette = "RdYlBu", direction = 1)+
  facet_wrap(~CV_type, ncol=4, labeller = labeller(CV_type = label_both))+
  AFH_theme+
  labs(y = "Mesh cutoff",
       x = "Mesh type")

ggplot(cv_sum_loglik%>% filter(CV_type=="Spatial_cluster"),
       aes(mesh, cutoff, fill = sum_loglik)) +
  geom_tile(color = "#00000000")+
  geom_text(aes(label = round(sum_loglik,digits = 0)))+
  scale_fill_distiller(palette = "RdYlBu", direction = 1)+
  facet_wrap(~CV_type, ncol=4, labeller = labeller(CV_type = label_both))+
  AFH_theme+
  labs(y = "Mesh cutoff",
       x = "Mesh type")

## show prediction ####
### compare the total index of biomass ###
index_comparison <- data.frame()
for (name in c("fit_gamma_regular_1","fit_gamma_data_fit_1",
               "fit_gamma_regular_1.852","fit_gamma_data_fit_1.852",
               "fit_gamma_regular_3")){
  prediction <- predict(list_fit[[name]], newdata = grid_proj %>%
                          select(X,Y,year,time) %>%
                          filter(time!=2024),
                        return_tmb_object = TRUE,se_fit = TRUE)
  
  index <- get_index(prediction, area = 0.25, bias_correct = TRUE)
  index$model <- name
  
  index_comparison <- rbind(index_comparison,index)
}
rm(prediction,index)

index_time <- c(0,0,0,0,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,
                0.3,0.3,0.3,0.3,0.4,0.4,0.4,0.4)
index_comparison$time <- index_comparison$time + index_time
index_comparison$model <- ordered(index_comparison$model,
                                  unique(index_comparison$model))

ggplot(index_comparison,aes(x = time, y = est, colour = model))+
  geom_point(show.legend = T)+
  scale_color_manual(values = as.character(
    c("#B15928","#6A3D9A","#DC902A","#C2090a","#158466"))
    )+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2,show.legend = F)+
  labs(x = "", y = "Biomass (tonnes)",
       title = "Biomass from survey", colour = "Model")+
  theme(plot.background = element_rect(fill = "#b4c7dc"),
        legend.background = element_rect(fill = "#b4c7dc"),
        panel.background = element_rect(fill = '#eeeeee', 
                                        color = '#000000',
                                        linewidth = 0.5),
        panel.grid.major = element_line(color = '#000000', 
                                        linetype = 'dotted'),
        panel.grid.minor = element_line(color = '#00000000', 
                                        linetype = 'dotted'))

### Map the prediction from the chosen model ###
# take the adapted column from the grid in prediction
predictions_gamma_regular <- predict(list_fit[[9]], newdata = grid_proj %>%
                                       select(X,Y,year,time) %>%
                                       filter(time!=2024),
                                     return_tmb_object = TRUE,se_fit = TRUE)

#Let’s make a small function to make maps
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year, nrow = 1) +
    coord_fixed()+
    theme(aspect.ratio = 3)
}

### show prediction density ###
plot_map(predictions_gamma_regular$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt")+
  ggtitle("Prediction (fixed effects + all random effects)")

### show spatial random effects ###
plot_map(predictions_gamma_regular$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

### show spatiotemporal random effects ###
plot_map(predictions_gamma_regular$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

### plot the density as sequential data ###
cut_min <- trunc(predictions_gamma_regular$data$est)
#cut_min <- trunc(predictions_gamma_data_fit$data$est)
cut_max <- cut_min+1
cuts <- paste(c("("), cut_min, c("-"), cut_max, c("]"), sep = "")

for (i in 1:length(cuts)){
  if (cuts[i]=="(8-9]"){
    cuts[i] <- ">8"
  }
  if (cuts[i]=="(9-10]"){
    cuts[i] <- ">8"
  }
  if (cuts[i]=="(10-11]"){
    cuts[i] <- ">8"
  }
}

predictions_gamma_regular$data$cuts <- as.factor(cuts)

ggplot(predictions_gamma_regular$data, aes(X, Y, fill = cuts)) +
  geom_raster() +
  scale_fill_brewer("log(density)", type = "seq", palette = "YlOrRd")+
  facet_wrap(~year, nrow = 1) +
  coord_fixed()+
  theme(aspect.ratio = 3)+
  ggtitle("Prediction (fixed effects + all random effects)")

###-###-###-###-###-###-###
# Covariates analysis ####
###-###-###-###-###-###-###

## check the correlation between the covariates and the density
ggplot(data=data_holotv, aes(x = bathy, y = density.t_km2)) +
  geom_point() +
  facet_wrap(~ year) +
  labs(x = "Depth (m)", y = "Density")

ggplot(data=data_holotv, aes(x = temp, y = density.t_km2)) +
  geom_point() +
  facet_wrap(~ year) +
  labs(x = "Bottom temp (°C)", y = "Density")

ggplot(data=data_holotv, aes(x = chla, y = density.t_km2)) +
  geom_point() +
  facet_wrap(~ year) +
  labs(x = "Mass concentration of chlorophyll a (mg/m3)", y = "Density")

## creation of the mesh with R-INLA ####
bnd <- INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                   convex = -0.04)
bnd2 = INLA::inla.nonconvex.hull(cbind(grid_proj$X, grid_proj$Y),
                                   convex = -0.15)
  
mesh_inla <- INLA::inla.mesh.2d(
    loc = as.matrix(as.data.frame(grid_bathy)[,c(1,2)]),
    boundary = list(bnd,bnd2),
    cutoff = 3, # minimum triangle edge length
    max.edge = c(3*3, 20*3), # inner and outer max triangle lengths
  ) # 3 is the distance chosen from the previous analysis
  
mesh <- make_mesh(data_estimate,c("X", "Y"), mesh = mesh_inla)

plot(mesh$mesh, main = NA, edge.color = "grey60", asp = 1)
points(data_holotv$X, data_holotv$Y, pch = 19, col = "red",cex = 0.3)


## fit model with different covariates ####

### prepare the data ###
# keep only the column needed
data_estimate <- as.data.frame(data_holotv) %>% 
  select(station,X,Y,density.t_km2,year,temp,bathy,chla)
grid_proj <- as.data.frame(grid_tot)[,1:6]
data_estimate$time <- data_estimate$year
grid_proj$time <- grid_proj$year

fit_none <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)
sanity(fit_none)

fit_bathy <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)
sanity(fit_bathy)

fit_temp <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)
sanity(fit_temp)

fit_temp_bathy <- sdmTMB(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID"
)
sanity(fit_temp_bathy)

## check the result ####

### residuals ###
par(mfrow = c(2,2))
resids_none <- residuals(fit_none) 
# randomized quantile residuals
qqnorm(resids_none, main = "Normal Q-Q Plot Gamma")
qqline(resids_none)

resids_bathy <- residuals(fit_bathy) 
# randomized quantile residuals
qqnorm(resids_bathy, main = "Normal Q-Q Plot Gamma Bathymetry")
qqline(resids_bathy)

resids_temp <- residuals(fit_temp) 
# randomized quantile residuals
qqnorm(resids_temp, main = "Normal Q-Q Plot Gamma Bottom temperature")
qqline(resids_temp)

resids_temp_bathy <- residuals(fit_temp_bathy) 
# randomized quantile residuals
qqnorm(resids_temp_bathy, main = "Normal Q-Q Plot Gamma Bathymetry and Bottom Temperature")
qqline(resids_temp_bathy)
par(mfrow = c(1,1))

### Observed vs predicted ###
observed <- fit_none$data$density.t_km2
predicted_none <- exp(predict(fit_none)$est)
predicted_bathy <- exp(predict(fit_bathy)$est)
predicted_temp <- exp(predict(fit_temp)$est)
predicted_temp_bathy <- exp(predict(fit_temp_bathy)$est)

AFH_theme <-  theme(panel.background = element_rect(fill = '#b4c7dc', 
                                                    color = '#b4c7dc',
                                                    linewidth = 1),
                    panel.grid.major = element_line(color = '#00000000', 
                                                    linetype = 'dotted'),
                    panel.grid.minor = element_line(color = '#00000000', 
                                                    linewidth = 2),
                    plot.background = element_rect(fill = "#b4c7dc"),
                    legend.background = element_rect(fill = "#b4c7dc"),
                    strip.background = element_rect(fill ="#dddddd"))

plot_grid(
  ggscatter(
    data.frame(observed,predicted_none),x = "observed", y = "predicted_none",
    add = "reg.line", title = "No covariates") +
    stat_cor(aes(label = paste(..r.label..)),label.x = 1000, label.y = 3000) +
    AFH_theme,
  ggscatter(
    data.frame(observed,predicted_bathy),x = "observed", y = "predicted_bathy",
    add = "reg.line", title = "Bathymetry") +
    stat_cor(aes(label = paste(..r.label..)),label.x = 1000, label.y = 3000) +
    AFH_theme,
  ggscatter(
    data.frame(observed,predicted_temp),x = "observed", y = "predicted_temp",
    add = "reg.line", title = "Bottom temperature") +
    stat_cor(aes(label = paste(..r.label..)),label.x = 1000, label.y = 3000) +
    AFH_theme,
  ggscatter(
    data.frame(observed,predicted_temp_bathy),x = "observed", y = "predicted_temp_bathy",
    add = "reg.line", title = "Bathymetry and Bottom temperature") +
    stat_cor(aes(label = paste(..r.label..)),label.x = 1000, label.y = 3000) +
    AFH_theme
)

## cross validation ####
# the previous result enable to choose the gamma models with 1 or 1.852 cutoff

### cross validation random ###
n_folds = 20

set.seed(22)
plan(multisession, workers = 8)
m_cv_none <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_none$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(22)
plan(multisession, workers = 8)
m_cv_bathy <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_bathy$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(22)
plan(multisession, workers = 8)
m_cv_temp <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_temp$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(22)
plan(multisession, workers = 8)
m_cv_temp_bathy <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  k_folds = n_folds
)
for (i in 1:n_folds){
  model <- m_cv_temp_bathy$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

m_cv_none$sum_loglik # total log-likelihood
m_cv_bathy$sum_loglik # total log-likelihood
m_cv_temp$sum_loglik # total log-likelihood
m_cv_temp_bathy$sum_loglik # total log-likelihood


### cross validation with spatial cluster ###

#cluster by time and space 
k <- 8
clust_2021 <- kmeans(data_estimate %>% filter(year==2021)%>%select(X,Y), k)$cluster
clust_2022 <- kmeans(data_estimate %>% filter(year==2022)%>%select(X,Y), k)$cluster
clust_2023 <- kmeans(data_estimate %>% filter(year==2023)%>%select(X,Y), k)$cluster
clust_2025 <- kmeans(data_estimate %>% filter(year==2025)%>%select(X,Y), k)$cluster
clust <- c(clust_2021,clust_2022+k,clust_2023+(k*2),clust_2025+(k*3))
data_estimate$clust <- clust


set.seed(11)
plan(multisession, workers = 8)
s_cv_none <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_none$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(11)
plan(multisession, workers = 8)
s_cv_bathy <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_bathy$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(11)
plan(multisession, workers = 8)
s_cv_temp <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_temp$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

set.seed(11)
plan(multisession, workers = 8)
s_cv_temp_bathy <- sdmTMB_cv(
  density.t_km2 ~ 1+as.factor(year)+smooth(temp)+smooth(bathy),
  data = data_estimate,
  mesh = mesh,
  family = Gamma(link = "log"),
  spatial = "on",
  time = "time",
  spatiotemporal = "IID",
  fold_ids = "clust"
)
for (i in 1:max(unique(clust))){
  model <- s_cv_temp_bathy$models[[i]]
  print(paste("Check model ",i,sep=""))
  sanity(model)
}

s_cv_none$sum_loglik # total log-likelihood
s_cv_bathy$sum_loglik # total log-likelihood
s_cv_temp$sum_loglik # total log-likelihood
s_cv_temp_bathy$sum_loglik # total log-likelihood

#Cross-validation is one of the best approaches that can be used to quantify model performance
#it seems that the model with best performance is without covariates

### check deviance and covariates effect ###
#deviance
1-deviance(fit_bathy)/deviance(fit_none)
1-deviance(fit_temp)/deviance(fit_none)
1-deviance(fit_temp_bathy)/deviance(fit_none)

# This version is in link space and the residuals are partial 
#randomized quantile residuals.
visreg::visreg(fit_bathy, "bathy")
visreg::visreg(fit_temp, "temp")
visreg::visreg(fit_temp_bathy, "bathy")
visreg::visreg(fit_temp_bathy, "temp")

## show prediction ####
### compare the total index of biomass ###
predicted_none <- predict(fit_none, newdata = grid_proj %>%
                            select(X,Y,year,time) %>%
                            filter(time!=2024),
                          return_tmb_object = TRUE,se_fit = TRUE)
predicted_bathy <- predict(fit_bathy, newdata = grid_proj %>%
                             select(X,Y,year,time,bathy) %>%
                             filter(time!=2024),
                           return_tmb_object = TRUE,se_fit = TRUE)
predicted_temp <- predict(fit_temp, newdata = grid_proj %>%
                            mutate(temp = bottomT) %>%
                            select(X,Y,year,time,temp) %>%
                            filter(time!=2024),
                          return_tmb_object = TRUE,se_fit = TRUE)
predicted_temp_bathy <- predict(fit_temp_bathy, newdata = grid_proj %>%
                                  mutate(temp = bottomT) %>%
                                  select(X,Y,year,time,temp,bathy) %>%
                                  filter(time!=2024),
                                return_tmb_object = TRUE,se_fit = TRUE)
index_comparison <- rbind(
  get_index(predicted_none, area = 0.25, bias_correct = TRUE),
  get_index(predicted_bathy, area = 0.25, bias_correct = TRUE),
  get_index(predicted_temp, area = 0.25, bias_correct = TRUE),
  get_index(predicted_temp_bathy, area = 0.25, bias_correct = TRUE)
)

index_comparison$model <- c(
  rep("fit_none",4),
  rep("fit_bathy",4),
  rep("fit_temp",4),
  rep("fit_temp_bathy",4)
)

index_time <- c(0,0,0,0,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3)
index_comparison$time <- index_comparison$time + index_time

index_comparison$model <- ordered(index_comparison$model, c("fit_none",
                                                            "fit_bathy",
                                                            "fit_temp",
                                                            "fit_temp_bathy"))

ggplot(index_comparison,aes(x = time, y = est, colour = model))+
  geom_point(show.legend = T)+
  scale_color_manual(values = as.character(
    c("#158466","#B15928","#6A3D9A","#DC902A"))
  )+
  geom_errorbar(aes( ymin = lwr, ymax = upr), width = 0.2,show.legend = F)+
  labs(x = "", y = "Biomass (tonnes)", title = "Biomass from survey")+
  theme(panel.background = element_rect(fill = '#eeeeee', 
                                                      color = '#000000',
                                                      linewidth = 0.5),
                      panel.grid.major = element_line(color = '#000000', 
                                                      linetype = 'dotted'),
                      panel.grid.minor = element_line(color = '#00000000', 
                                                      linetype = 'dotted'),
                      plot.background = element_rect(fill = "#b4c7dc"),
                      legend.background = element_rect(fill = "#b4c7dc"),
                      strip.background = element_rect(fill ="#dddddd"))

### Map the prediction from the chosen model ###
# take the adapted column from the grid in prediction

#Let’s make a small function to make maps
plot_map <- function(dat, column) {
  ggplot(dat, aes(X, Y, fill = {{ column }})) +
    geom_raster() +
    facet_wrap(~year, nrow = 1) +
    coord_fixed()+
    theme(aspect.ratio = 3)
}

### show prediction density ###
plot_map(predicted_none$data, exp(est)) +
  scale_fill_viridis_c(trans = "sqrt")+
  ggtitle("Prediction (fixed effects + all random effects)")

### show spatial random effects ###
plot_map(predicted_none$data, omega_s) +
  ggtitle("Spatial random effects only") +
  scale_fill_gradient2()

### show spatiotemporal random effects ###
plot_map(predicted_none$data, epsilon_st) +
  ggtitle("Spatiotemporal random effects only") +
  scale_fill_gradient2()

### plot the density as sequential data ###
cut_min <- trunc(predicted_none$data$est)
#cut_min <- trunc(predictions_gamma_data_fit$data$est)
cut_max <- cut_min+1
cuts <- paste(c("("), cut_min, c("-"), cut_max, c("]"), sep = "")

for (i in 1:length(cuts)){
  if (cuts[i]=="(8-9]"){
    cuts[i] <- ">8"
  }
  if (cuts[i]=="(9-10]"){
    cuts[i] <- ">8"
  }
  if (cuts[i]=="(10-11]"){
    cuts[i] <- ">8"
  }
}

predicted_none$data$cuts <- as.factor(cuts)

ggplot(predicted_none$data, aes(X, Y, fill = cuts)) +
  geom_raster() +
  scale_fill_brewer("log(density)", type = "seq", palette = "YlOrRd")+
  facet_wrap(~year, nrow = 1) +
  coord_fixed()+
  scale_x_continuous(breaks = seq(555, 568, by = 5))+
  theme(aspect.ratio = 3,
        panel.background = element_rect(fill = '#eeeeee', 
                                        color = '#000000',
                                        linewidth = 0.5),
        panel.grid.major = element_line(color = '#000000', 
                                        linetype = 'dotted'),
        panel.grid.minor = element_line(color = '#00000000', 
                                        linetype = 'dotted'),
        plot.background = element_rect(fill = "#b4c7dc"),
        legend.background = element_rect(fill = "#b4c7dc"),
        strip.background = element_rect(fill ="#dddddd"))+
  ggtitle("Predictions (fixed and random effects) from the model without covariates")

