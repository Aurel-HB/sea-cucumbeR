# this code is a synthetic version of the first exploratory_PPP_analysis.R
# it used the chapter 2 of Spatial Point Process from Baddeley 2015 to test 
# if the distribution of the sea cucumber follow complete randomness

#############
#load packages
#############
library(here)
library(sf)
library(spatstat)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(FactoMineR)
library(factoextra)
library(cluster)
library(NbClust)

#############
#load data
#############

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
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


########################################
# test of Complete Spatial Randomness
########################################

# homogeneity of the intensity -> Quadrat counting test of homogeneity
# quadrat.test(point_pattern, nx,ny ) 
# where the expected counts greater than 5 for all quadrats

# measuring dependence between points in a point pattern, 
# using the concept of correlation -> DCLF (Diggle-Cressie-Loosmore-Ford) test
# based on the mean squared deviation between the empirical summary function and
# the theoretical function over a range of distance values
# significant if p-value < 1/(nsim+1)

# summary statistics based on the spacing between points in a point pattern
# Hopkins-Skellam test of Complete Spatial Randomness
# based on comparing nearest-neighbour distances with point-event distances
# random:A=1, cluster:A<1, regularity:A>1

CSR_PP_2025 <- HOLOTVSPM2025[,c("STN","X.centroid_track","Y.centroid_track")]
CSR_PP_2025 <- CSR_PP_2025 %>% filter(STN %in% names(list_PPP))

CSR_PP_2025 <- CSR_PP_2025 %>% 
  mutate(Quadrat = NA) %>% 
  mutate(DCLF = NA) %>% 
  mutate(Hopkins = NA)

for (stn in names(list_PPP)){
  X <- list_PPP[[stn]]
  # Quadrat counting test of homogeneity
  tS <- quadrat.test(X,1,10)
  tS
  ## Diggle-Cressie-Loosmore-Ford
  dclf <- dclf.test(X, Lest, nsim=99, use.theo=TRUE)
  dclf
  # Hopkins-Skellam index
  hopkins_index <- hopskel(X)
  hopkins_index
  
  CSR_PP_2025[CSR_PP_2025$STN==stn,]$Quadrat <- tS$p.value
  CSR_PP_2025[CSR_PP_2025$STN==stn,]$DCLF <- dclf$p.value
  CSR_PP_2025[CSR_PP_2025$STN==stn,]$Hopkins <- hopkins_index
}

# ACP ####
CSR_acp <- CSR_PP_2025 %>% 
  select(Quadrat, DCLF, Hopkins)
row.names(CSR_acp) <- CSR_PP_2025$STN


acp3 = FactoMineR::PCA(CSR_acp,scale.unit=F, ncp=5, graph=F)
acp3
fviz_pca_biplot(acp3,
                addEllipses = T,      # Add ellipses for categories
                repel = F,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")

#check the ACP
fviz_eig(acp3, addlabels = TRUE)
# Graph of the variables
fviz_pca_var(acp3, col.var = "black")
#contribution of each variable
fviz_cos2(acp3, choice = "var", axes = 1:2)
fviz_pca_var(acp3, col.var = "cos2",
             gradient.cols = c("black", "orange", "green"),
             repel = TRUE)

# Hierarchical clustering ####

### Scale the data
CSR_hc <- as.data.frame(scale(CSR_acp))

distance_matrix <- dist(CSR_hc)
cluster_result <- hclust(distance_matrix, method = "ward.D2")

### be careful of highly correlated variables
#### Calculate correlation matrix
cor_matrix <- cor(CSR_hc)
####  Find highly correlated pairs (e.g., |cor| > 0.9)
high_cor <- which(abs(cor_matrix) > 0.9 & cor_matrix < 1, arr.ind = TRUE)
####  Remove one variable from each highly correlated pair
#CSR_hc <- CSR_hc[, -unique(high_cor[,2])[c(1,3)]]

# Determine the Optimal Number of Clusters
### Plot WSS for different k and look for the "elbow" point
wss <- sapply(1:20, function(k) {
  km <- kmeans(CSR_hc, centers = k, nstart = 10)
  km$tot.withinss
})
#plot(1:20, wss, type = "b", xlab = "Number of Clusters",
#     ylab = "Within-cluster Sum of Squares")
#### First differences
first_diff <- diff(wss)
#### Second differences
second_diff <- diff(first_diff)
#### Optimal k: where second_diff is maximized
optimal_k_wss <- which.max(second_diff) + 2  # +2 because of diff offset

### Maximize the average silhouette width.
silhouette_scores <- sapply(2:20, function(k) {
  km <- kmeans(CSR_hc, centers = k, nstart = 10)
  sil <- silhouette(km$cluster, dist(CSR_hc))
  mean(sil[, "sil_width"])
})
#plot(2:20, silhouette_scores, type = "b", xlab = "Number of Clusters",
#     ylab = "Average Silhouette Width")
#### Optimal k: where silhouette score is maximized
optimal_k_silhouette <- which.max(
  silhouette_scores)+1 #+1 because sapply starts at 2
#print(paste("Optimal k (Silhouette):", optimal_k_silhouette))

### Example code for NbClust
index <- c("kl", "ch", "hartigan",
           "gap","dunn","pseudot2","duda")
nbclust_result <- list()
nbclust_best.nc <- c()
for (i in 1:length(index)){
  nbclust_result[[index[i]]] <- NbClust(CSR_hc, distance = "euclidean", min.nc = 2,
                                        max.nc = 20, method = "ward.D2",
                                        index = index[i])
  nbclust_best.nc <- c(nbclust_best.nc,
                       nbclust_result[[index[i]]][["Best.nc"]][1])
}
index <- c(index, "wss", "silhouette")
nbclust_best.nc <- c(nbclust_best.nc, optimal_k_wss, optimal_k_silhouette)
nbclust_data <- data.frame(index,nbclust_best.nc)


# plot the dendogramme
nb_cluster <- 2 #or nb_cluster <- 3
plot(cluster_result, main = "Dendrogramme des processus")
rect.hclust(cluster_result, k = nb_cluster, border = 2:4) 

# explore cluster
# c(100,210,161,106,220,143,110,193,129,222,127,144)
cluster2 <- as.character(c(100,210,161,106,220,143,110,193,129,222,127,144))
print(CSR_PP_2025 %>% filter(STN %in% cluster2))

CSR_PP_2025$STN[CSR_PP_2025$DCLF>0.1]
CSR_PP_2025$STN[CSR_PP_2025$DCLF>0.1] %in% cluster2
print(CSR_PP_2025 %>% filter(STN %in% CSR_PP_2025$STN[CSR_PP_2025$DCLF>0.05]))

CSR_PP_2025$STN[CSR_PP_2025$Quadrat>0.05]
CSR_PP_2025$STN[CSR_PP_2025$Quadrat>0.05] %in% cluster2
print(CSR_PP_2025 %>% filter(STN %in% CSR_PP_2025$STN[CSR_PP_2025$Quadrat>0.05]))

ggplot()+
  geom_point(aes(list_PPP[["222"]]$x, list_PPP[["222"]]$y))
