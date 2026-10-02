# This code analyse the output from the script Point_Process_Modelling

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
library(FactoMineR)
library(factoextra)
library(gridExtra)
library(RColorBrewer)
library(colorspace)
library(purrr)
library(tibble)

### ### ### ###
#Load data ####
### ### ### ###
list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))

scoring_tot <- readRDS(paste(
  here(),"/point_process/Output/scoring_tot.rds",sep=""))

residual_tot <- readRDS(paste(
  here(),"/point_process/Output/residual_tot.rds",sep=""))

smoothed_residuals_raw <- readRDS(
  paste(here(),"/point_process/Output/smoothed_residuals.rds",sep=""))

process_list <- readRDS(paste(
  here(),"/point_process/Output/process_list.rds",sep=""))

nsim_valid_LGCP <- readRDS(paste(
  here(),"/point_process/Output/nsim_valid_LGCP.rds",sep=""))

### ### ### ###
#Analysis  ####
### ### ### ###

## Prepare data ####
metrics <- scoring_tot %>% 
  mutate(cor.coef = as.numeric(residual_tot$cor.coef)) %>%
  mutate(lm.coef = as.numeric(residual_tot$lm.coef)) %>%
  mutate(p_value_envelope = as.numeric(residual_tot$p_value_envelope))

plot_PPP <- function(PPP){
  ggplot(data = data.frame(V1 = as.numeric(PPP$x),
                           V2 = as.numeric(PPP$y),
                           name = rep("PPP",PPP$n)))+
    geom_point(aes(x=V1,y=V2),color="black")+
    ylab("")+xlab("")+
    theme(#axis.ticks = element_blank(),
          #axis.text = element_blank()
    )
}

## Select model ####
data_select <- data.frame(
  STN = unique(metrics$STN),
  lambda_score = NA,
  K_score = NA,
  cor.coef = NA,
  lm.coef = NA#,p_value_envelope = NA
)

for (stn in unique(metrics$STN)){
  data <- metrics %>% filter(STN == stn) %>%
    filter(!is.na(lambda_score)) %>%
    filter(!is.na(K_score))
  
  data_select$lambda_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$lambda_score), data$lambda_score)]
  
  data_select$K_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$K_score), data$K_score)]
  
  data_select$cor.coef[grep(stn, data_select$STN)] <- 
    data$model[grep(min(abs(1-data$cor.coef)), abs(1-data$cor.coef))]
  
  data_select$lm.coef[grep(stn, data_select$STN)] <- 
    data$model[grep(min(abs(1-data$lm.coef)), abs(1-data$lm.coef))]
  
  #data_select$p_value_envelope[grep(stn, data_select$STN)] <- 
  #  data$model[grep(max(data$p_value_envelope), data$p_value_envelope)]
  
}

### Work on the residual to homogeneise the information ####

#### delete all the model that not pass the envelope test
checked_residual <- metrics %>% filter(p_value_envelope > 0.025)
length(unique(checked_residual$STN))==61
# if false identify the stations
warning_residual_station <- c()
for (stn in unique(metrics$STN)){
  if (stn %in% unique(checked_residual$STN)){
  } else {
    warning_residual_station <- c(warning_residual_station, stn)
  }
}

# as some simulations have not be taking account for LGCP,
# we adapt the p_value filter to this particular case
residual_LGCP <- checked_residual %>% filter(model=="LGCP") %>%
  filter(STN %in% nsim_valid_LGCP$stn[nsim_valid_LGCP$nsim_valid<39]) %>%
  mutate(nsim_valid = nsim_valid_LGCP$nsim_valid[nsim_valid_LGCP$nsim_valid<39])

residual_LGCP$test_value <- 1/(residual_LGCP$nsim_valid+1)

residual_LGCP <- residual_LGCP %>% filter(p_value_envelope <= test_value)
# if null there is no problem 
# if not null, use the command below to delete the stations that finally not pass the test
# not pass the test:
### checked_residual[
###   grep(residual_LGCP$lambda_score[1],checked_residual$lambda_score),] <- NA
checked_residual <- checked_residual[!is.na(checked_residual$STN),]

# test again the missing station
warning_residual_station <- c()
for (stn in unique(metrics$STN)){
  if (stn %in% unique(checked_residual$STN)){
  } else {
    warning_residual_station <- c(warning_residual_station, stn)
  }
}

#### separate case when cor.coef and lm.coef are agree and the others
data_select_2 <- data.frame(
  STN = unique(checked_residual$STN),
  cor.coef = NA,
  lm.coef = NA
)

for (stn in unique(checked_residual$STN)){
  
  data <- checked_residual %>% filter(STN == stn)

  data_select_2$cor.coef[grep(stn, data_select_2$STN)] <- 
    data$model[grep(min(abs(1-data$cor.coef)), abs(1-data$cor.coef))]
  
  data_select_2$lm.coef[grep(stn, data_select_2$STN)] <- 
    data$model[grep(min(abs(1-data$lm.coef)), abs(1-data$lm.coef))]
  
}


check_residual <- data_select_2 %>% filter(cor.coef!=lm.coef)
data_select_3 <- data_select_2 %>% filter(cor.coef==lm.coef) %>%
  select(STN,cor.coef)
names(data_select_3)[2] <- "residual_validation"

# visual check of the corresponding station
for (i in 1:nrow(check_residual)){
  stn <- check_residual$STN[i]
  plot1 <- smoothed_residuals_raw[[paste("fit",check_residual$cor.coef[i],
                                         stn,sep="_")]]
  plot2 <- smoothed_residuals_raw[[paste("fit",check_residual$lm.coef[i],
                                         stn,sep="_")]]
  grid.arrange(plot1, plot2, ncol=2)
}
# 129 -> lwppp  
# 151 -> ihP 
# 183 -> ihP
# 185 -> LGCP  
# 190 -> LGCP 
# 193 -> ihP
# 214 -> LGCP
# 226 -> lwppp
data_select_3 <- rbind(data_select_3,
                       data.frame(
                         STN = c(129,151,190,
                                 193,226,214,
                                 185,183),
                         residual_validation = c("lwppp","ihP","LGCP",
                                                 "ihP","lwppp","LGCP",
                                                 "LGCP","ihP")
                       ))
#data_select_3 <- rbind(data_select_3,
#                    data.frame(
#                          STN = c(112,143,226,185,190,193),
#                          residual_validation = c("ihP","ihP","lwppp",
#                                                   "LGCP","ihP","ihP")
#                      ))

# add the station that not pass the envelope test
data_select_3 <- rbind(data_select_3,
                       data.frame(
                         STN = c(123,133,141,155,157,176,216),
                         residual_validation = c(NA,NA,NA,NA,NA,NA,NA)
                       ))

#data_select_3 <- rbind(data_select_3,
#                     data.frame(
#                           STN = c(133,137,187),
#                           residual_validation = c(NA,NA,NA)
#                       ))

# merge data_select with the result of residual validation
data_select_final <- merge(
  data_select[,(1:3)], data_select_3
)

# compare the model choose by score under previous residual filter
# this part was a test but not keep for instance ###
#for (stn in unique(checked_residual$STN)){
#  data <- checked_residual %>% filter(STN == stn) %>%
#    filter(!is.na(lambda_score))
#  
#  data_select$lambda_score[grep(stn, data_select$STN)] <- 
#    data$model[grep(min(data$lambda_score), data$lambda_score)]
#  
#  data_select$K_score[grep(stn, data_select$STN)] <- 
#    data$model[grep(min(data$K_score), data$K_score)]
#}
#names(data_select)[2:3] <- c("lambda_score_res","K_score_res") 
#data_select_compare <- merge(
#  data_select[,(1:3)],data_select_final
#)


## Select one model per station ####
data_select_final$selected_model <- NA
for (i in 1:nrow(data_select_final)){
  IHP <- 0
  LGCP <- 0
  LWPPP <- 0
  for (j in 2:4){
    if(!is.na(data_select_final[i,j])){
      if(data_select_final[i,j]=="ihP"){
        IHP <-IHP + 1
      }
      if(data_select_final[i,j]=="LGCP"){
        LGCP <-LGCP + 1
      }
      if(data_select_final[i,j]=="lwppp"){
        LWPPP <-LWPPP + 1
      }
    }
  }
  data <- data.frame(
    model = c("ihP","LGCP","lwppp"),
    occurence = c(IHP,LGCP,LWPPP)
  )
  if (max(data$occurence)==1){
    data_select_final$selected_model[i] <- "equality"
  } else {
    data_select_final$selected_model[i] <- data$model[
      data$occurence==max(data$occurence)]
  }
}
rm(data,i,j,IHP,LGCP,LWPPP)

# select each selected model in list
ihP_process_list <- list()
LGCP_process_list <- list()
lwppp_process_list <- list()
for (i in 1:nrow(data_select_final)){
  if (data_select_final$selected_model[i]=="ihP"){
    name <- paste("fit","ihP",data_select_final$STN[i],sep = "_")
    ihP_process_list[[name]] <- process_list[[name]]
  }
  if (data_select_final$selected_model[i]=="LGCP"){
    name <- paste("fit","LGCP",data_select_final$STN[i],sep = "_")
    LGCP_process_list[[name]] <- process_list[[name]]
  }
  if (data_select_final$selected_model[i]=="lwppp"){
    name <- paste("fit","lwppp",data_select_final$STN[i],sep = "_")
    lwppp_process_list[[name]] <- process_list[[name]]
  }
}
rm(name)

#saveRDS(data_select_final, paste(
#  here(),"/point_process/Output/model_selection.rds",sep="") )

data_select_final <- readRDS(paste(
   here(),"/point_process/Output/model_selection.rds",sep=""))

## descriptve statistics ####
### Principal Component Analysis ####
row.names(metrics) <- paste(metrics$STN,metrics$model,sep="_")
res.acp <- PCA(metrics[,(3:7)] %>% filter(!is.na(metrics$K_score)),
               scale.unit=F, ncp=5, graph=F)
fviz_eig(res.acp, addlabels = TRUE)
fviz_pca_biplot(res.acp,
                addEllipses = T,      # Add ellipses for categories
                repel = T,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")


### Multiple Correspondence Analysis ####
#Work with the data_select_final
row.names(data_select_final) <- data_select_final$STN
res.mca <- MCA(data_select_final[,(2:4)])
fviz_eig(res.mca, addlabels = TRUE)
fviz_mca_biplot(res.mca,
                addEllipses = T,      # Add ellipses for categories
                repel = T,            # Avoid label overlap
                title = "MCA Biplot",
                col.var = "darkblue",
                col.ind = "#5577aa")+
  theme() +
  labs(title = "Customized MCA Biplot")

data_mod_num <- data_select_final
data_mod_num <- data_mod_num %>% filter(!is.na(residual_validation))
for (i in 1:nrow(data_mod_num)){ 
  for (j in 1:ncol(data_mod_num)){
    if(data_mod_num[i,j]=="ihP"){
      data_mod_num[i,j] <- 0
    }
    if(data_mod_num[i,j]=="LGCP"){
      data_mod_num[i,j] <- 0.5
    }
    if(data_mod_num[i,j]=="lwppp"){
      data_mod_num[i,j] <- 1
    }
  }
}

data_mod_num <- data_mod_num %>% 
  mutate(lambda_score = as.numeric(lambda_score)) %>% 
  mutate(K_score = as.numeric(K_score)) %>% 
  mutate(residual_validation = as.numeric(residual_validation))
res.acp <- PCA(data_mod_num[,(2:4)],
               scale.unit=F, ncp=5, graph=F)
fviz_eig(res.acp, addlabels = TRUE)
fviz_pca_biplot(res.acp,
                addEllipses = T,      # Add ellipses for categories
                repel = T,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")

### Hierarchical Clustering the data-table ####
library(cluster)

# df = your data frame with categorical (factor) columns,
# and optionally numeric ones too
df <- data_select_final
df$lambda_score <- as.factor(df$lambda_score)
df$K_score <- as.factor(df$K_score)
df$residual_validation <- as.factor(df$residual_validation)

# Step 1: compute Gower dissimilarity matrix
gower_dist <- daisy(df[,(2:4)], metric = "gower")

# Step 2: hierarchical clustering on the dissimilarity matrix
hc <- hclust(gower_dist, method = "average")  # or "complete", "single"

# Step 3: plot dendrogram
plot(hc, labels = FALSE, main = "Hierarchical Clustering (Gower distance)")

# Step 4: cut the tree into k clusters
clusters <- cutree(hc, k = 16)
df$cluster <- clusters

def_cluster <- data.frame(df[1,], count = length(grep(TRUE,df$cluster==1))) 
clust <- 1
for (i in 2:nrow(df)){
  if (clust < df$cluster[i]){
    clust <- clust + 1
    def_cluster <- rbind(def_cluster,
                         data.frame(df[i,], 
                                    count=length(grep(TRUE,df$cluster==clust))))
  }
}
def_cluster <- def_cluster %>% select(cluster,lambda_score,K_score,
                                 residual_validation, count) 

### Extract the combinations ####
# Extract all the possible combination and the number of time they appear

data_combi <- as.data.frame(data_map) %>% 
  select(station,lambda_score,K_score,residual_validation)

data_combi <- data_combi %>% 
  mutate(combine = paste(lambda_score,K_score,
                         residual_validation,sep = "/")) %>%
  mutate(combine = as.factor(combine))

data_combi <- as.data.frame(summary(data_combi$combine))
data_combi$combine = row.names(data_combi)

data <- data.frame( lambda_score = NA,
                    K_score = NA,
                    residual_validation = NA,
                    count = data_combi$`summary(data_combi$combine)`)
for (i in 1:nrow(data)){
  data$lambda_score[i] <- gsub("(.*)\\/.*/.*","\\1", data_combi$combine[i], 
                               perl=T)
  data$K_score[i] <- gsub(".*/(.*)\\/.*","\\1", data_combi$combine[i], 
                          perl=T)
  data$residual_validation[i] <- gsub(".*/.*/(.*)","\\1", 
                                      data_combi$combine[i], perl=T)
}
data_combi <- data
rm(data)

write.csv(data_combi,
          paste(here(),"/point_process/Output/data_combinaison.csv",
                sep = ""),
          row.names = FALSE)

# Map the result ####
# load an prepare data
data_2025 <- readRDS(
  paste(here(),"/process_spatiotemp_data/Data/processed/data_abun_tot.rds",
        sep = "")
)
data_2025 <- data_2025 %>% filter(year==2025)
data_map <- merge(data_2025,data_select_final,by.x="station",by.y="STN")
data_map <- data_map %>% 
  mutate(residual_validation = ifelse(is.na(residual_validation),"NA",
                                      residual_validation))

calcul_area <- readRDS(paste(here(),
                             "/workshop_spatiotemp/Data/study_calcul_area.rds",
                             sep=""))

ggplot(data_map)+
  geom_sf(aes(color=lambda_score), size = 5)+
  scale_color_brewer("lambda_score", type = "qua", palette = "Dark2")+
  geom_sf(data=calcul_area, fill = "#11111111")+
  theme(aspect.ratio = 2,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  labs(title = "Model chosen by lambda score by station")

ggplot(data_map)+
  geom_sf(aes(color=K_score), size = 5)+
  scale_color_brewer("K_score", type = "qua", palette = "Dark2")+
  geom_sf(data=calcul_area, fill = "#11111111")+
  theme(aspect.ratio = 2,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  labs(title = "Model chosen by K score by station")

ggplot(data_map)+
  geom_sf(aes(color=residual_validation), size = 5)+
  scale_color_brewer("residual_validation", type = "qua", palette = "Dark2")+
  geom_sf(data=calcul_area, fill = "#11111111")+
  theme(aspect.ratio = 2,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  labs(title = "Model chosen by residual validation by station")

ggplot(data_map)+
  geom_sf(aes(color=selected_model), size = 5)+
  scale_color_brewer("selected_model", type = "qua", palette = "Dark2")+
  geom_sf(data=calcul_area, fill = "#11111111")+
  theme(aspect.ratio = 2,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "white"))+
  labs(title = "Model chosen by station")

# Show the point pattern by model fitted ####
data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

# all
ggplot(data = data_position)+
geom_point(aes(x=X,y=Y))+
facet_wrap(~station, nrow = 4,scales="free_y")+
ylab("")

# only ihP
ggplot(data = data_position %>% 
         filter(station %in% data_select_final$STN[
           data_select_final$selected_model=="ihP"]))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 3,scales="free_y")+
  ylab("")

# only LGCP
ggplot(data = data_position %>% 
         filter(station %in% data_select_final$STN[
           data_select_final$selected_model=="LGCP"]))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 2,scales="free_y")+
  ylab("")

# only lwppp
ggplot(data = data_position %>% 
         filter(station %in% data_select_final$STN[
           data_select_final$selected_model=="lwppp"]))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 1,scales="free_y")+
  ylab("")

# only equality
ggplot(data = data_position %>% 
         filter(station %in% data_select_final$STN[
           data_select_final$selected_model=="equality"]))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 1,scales="free_y")+
  ylab("")

# Work on the result from the model ####

## LGCP ####

### extract effect for one LGCP ####
fit <- LGCP_process_list[["fit_LGCP_104"]]
tr <- predict(LGCP_process_list[["fit_LGCP_104"]], type = "trend")
data_test <- data.frame(
  x = sort(rep(tr$xcol, 128)),
  y = rep(tr$yrow, 128),
  value = as.numeric(tr$v)
)
ggplot(data_test)+
  geom_point(aes(x = x, y = y, colour=value))+
  scale_color_gradient(low = "#eeee44",
                       #mid = "#bb1166",
                       high = "blue")+
  ggtitle("Tendance (partie déterministe)")+
  theme(aspect.ratio=6)

sims <- simulate(fit, nsim = 4, saveLambda = TRUE)
Lam  <- attr(sims[[1]], "Lambda")     # realized random intensity
data_test <- data.frame(
  x = sort(rep(Lam$xcol, 128)),
  y = rep(Lam$yrow, 128),
  value = as.numeric(Lam$v)
)
ggplot(data_test)+
  geom_point(aes(x = x, y = y, colour=value))+
  scale_color_gradient(low = "#eeee44",
                       #mid = "#bb1166",
                       high = "blue")+
  ggtitle("Intensité aléatoire Lambda(u)")+
  theme(aspect.ratio=6)

# the gaussian field only, without the trend :
Z <- log(Lam / tr)                   # = Z(u) - sigma2/2 pour un LGCP
data_test <- data.frame(
  x = sort(rep(Z$xcol, 128)),
  y = rep(Z$yrow, 128),
  value = as.numeric(Z$v)
)
ggplot(data_test)+
  geom_point(aes(x = x, y = y, colour=value))+
  scale_color_gradient(low = "#eeee44",
                       #mid = "#bb1166",
                       high = "blue")+
  ggtitle("Réalisation du champ gaussien latent")+
  theme(aspect.ratio=6)

### comparison between all LGCP ####
extract_kppm <- function(f, nm) {
  s <- summary(f)
  mp <- f$modelpar             # var, scale (sometimes nu if estimated)
  cm <- f$covmodel              # list: model, margs (contain nu fixed)
  
  # because mu is somtimes image and not numeric
  mu_val <- f$mu
  mean_logint <- if (is.im(mu_val)) {
    mean(mu_val, na.rm = TRUE)      # spatial mean of log-intensity
  } else {
    as.numeric(mu_val)
  }
  return(
    tibble(
      pattern      = nm,
      n_points     = npoints(f$X),
      intercept    = coef(f)[["(Intercept)"]],
      slope_y      = coef(f)[["y"]],
      var          = unname(mp["sigma2"]),
      scale        = unname(mp["alpha"]),
      nu           = if (!is.null(cm$margs$nu)) cm$margs$nu else NA_real_,
      mean_logint  = mean_logint,                     # mean of log-field
      method       = f$Fit$method,              # "mincon", "clik2", "palm"...
      covmodel     = cm$model,
      long_window  = diameter(f$X$window),
      scale_ratio  = scale / long_window
    ))
}

tab <- data.frame()
for (i in 1:length(LGCP_process_list)){
  tab <- rbind(tab,
               extract_kppm(
                 LGCP_process_list[[i]],
                 names(LGCP_process_list[i])
               ))
}
rm(i)
# check the variance and scale parameter to detect impossible numeric value
# when var~0 = the LGCP degenerates into a simple Poisson process
# when scale~inf = confusion large-scale structure with deterministic trends
tab_verified <- tab %>% filter(scale_ratio < 1) %>% filter(var > 0.001)
LGCP_no_check <- substr(tab$pattern[tab$scale_ratio >=1 ],10,12)

data_select_final[data_select_final$STN %in% LGCP_no_check,]
metrics[metrics$STN %in% LGCP_no_check,]
# only for stn 127 all the others have one parameter choosing ihP
# the selected model has therefore been replaced by ihP for this stn
data_select_final_corrected <- data_select_final %>% 
  mutate(selected_model = ifelse(STN %in% LGCP_no_check, "ihP", selected_model))

plot_PPP(list_PPP[["127"]])
data_select_final_corrected$selected_model[data_select_final_corrected$STN==
                                             "127"] <- NA
data_map <- merge(data_2025,data_select_final_corrected,
                  by.x="station",by.y="STN")
data_map <- data_map %>% 
  mutate(residual_validation = ifelse(is.na(residual_validation),"NA",
                                      residual_validation)) %>%
  mutate(selected_model = ifelse(is.na(selected_model),"NA",
                                 selected_model))
mypalette <- c("#ff71a2","#ffb069","#bdd2f9","#97dda9","#c3c3c3")
ggplot(data_map)+
  geom_sf(data=calcul_area, fill = "#11111100")+
  geom_sf(aes(fill=selected_model), size = 5,shape = 21)+
  scale_fill_manual(values = mypalette) +
  #scale_color_brewer("selected_model", type = "qua", palette = "Dark2")+
  theme(aspect.ratio = 2,
        legend.title = element_blank(),
        title = element_text(color = "black",face = "bold"),
        plot.title = element_text( size = 12, hjust = 0.5),
        plot.subtitle = element_text(size = 8,hjust = 0.5),
        panel.border = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey90"),
        panel.background = element_rect(fill = "lightblue"),
        panel.grid.minor = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "grey90"))+
  labs(title = "Model chosen by station")

#fit_lgcp_y_bounded <- kppm(X ~ y,
#                           clusters = "LGCP",
#                           model = "matern",
#                           statistic = "pcf",
#                           covfunargs = list(nu = 0.3),
#                           startpar = c(var = 1, scale = 2),   # sensible local-scale start
#                           control = list(
#                             method = "L-BFGS-B",
#                             lower  = c(var = 1e-4, scale = 0.1),   # e.g. > pixel/measurement error
#                             upper  = c(var = 50,   scale = 50)     # well below your 600 m domain
#                           ))

#for more transparence
#fit <- lgcp.estpcf(X, 
#                   startpar = c(var = 1, scale = 2),
#                   covmodel = list(model = "matern", nu = 0.3),
#                   q = 1/4, p = 2,
#                   rmin = NULL, rmax = 20,   # restrict the r-range used in the contrast — 
#                   # forces the fit to focus on local scales
#                   control = list(method = "L-BFGS-B",
#                                  lower = c(1e-4, 0.1),
#                                  upper = c(50, 50)))


## ihP ####
ihP_model_caract <- data.frame()
for (i in 1:length(ihP_process_list)){
  ihP_model_caract <- rbind(ihP_model_caract,
                            data.frame(
                              stn = substr(names(ihP_process_list[i]),start = 9,
                                           stop = 11),
                              model = "ihP",
                              formula = deparse(
                                ihP_process_list[[i]][["trend"]][[2]])
                            )
  )
}
rm(i)

ihP_model_caract <- merge(ihP_model_caract, data_select_final_corrected,
                          by.x = "stn", by.y = "STN")

## lwppp ####
lwppp_model_caract <- data.frame()
for (i in 1:length(lwppp_process_list)){
  lwppp_model_caract <- rbind(lwppp_model_caract,
                            data.frame(
                              stn = substr(names(lwppp_process_list[i]),
                                           start = 11,
                                           stop = 13),
                              model = "lwppp",
                              formula = deparse(
                                lwppp_process_list[[i]][["trend"]][[2]])
                            )
  )
}
rm(i)

