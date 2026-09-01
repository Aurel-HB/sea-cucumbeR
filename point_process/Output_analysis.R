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
  ggplot(data = as.data.frame(cbind(PPP$x,PPP$y,rep("PPP",PPP$n))))+
    geom_point(aes(x=V1,y=V2),color="black")+
    ylab("")+xlab("")+
    theme(axis.ticks = element_blank(),
          axis.text = element_blank())
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

# the stations that finally not pass the test are delete from the table
checked_residual[grep(
  residual_LGCP$lambda_score[1],checked_residual$lambda_score),
  ] <- NA

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
# 125 ->  ihP  
# 127 -> LGCP 
# 159 -> lwppp
# 173 ->  ihP  
# 190 -> LGCP 
# 193 -> LGCP 
data_select_3 <- rbind(data_select_3,
                       data.frame(
                         STN = c(125,127,159,173,190,193),
                         residual_validation = c("ihP","LGCP","lwppp",
                                                 "ihP","LGCP","LGCP")
                       ))

# add the station that not pass the envelope test
data_select_3 <- rbind(data_select_3,
                       data.frame(
                         STN = c(157,216,179),
                         residual_validation = c(NA,NA,NA)
                       ))

# merge data_select with the result of residual validation
data_select_final <- merge(
  data_select[,(1:3)], data_select_3
)

# compare the model choose by score under previous residual filter
for (stn in unique(checked_residual$STN)){
  data <- checked_residual %>% filter(STN == stn) %>%
    filter(!is.na(lambda_score))
  
  data_select$lambda_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$lambda_score), data$lambda_score)]
  
  data_select$K_score[grep(stn, data_select$STN)] <- 
    data$model[grep(min(data$K_score), data$K_score)]
}
names(data_select)[2:3] <- c("lambda_score_res","K_score_res") 
data_select_compare <- merge(
  data_select[,(1:3)],data_select_final
)

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
clusters <- cutree(hc, k = 18)
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

ggplot(merge(data_2025,df,by.x="station",by.y="STN") %>%
         mutate(cluster = as.factor(cluster)))+
  geom_sf(aes(color=cluster), size = 5)+
  scale_color_manual(values = pal)+
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
  labs(title = "Cluster for each onfiguration by station")
