# this code use the chapter 2 of Spatial Point Process from Baddeley
# the aim is to explore the different haul to see the PPP that define the 
# sea cucumber population.

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

#############
#load data
#############
data_position <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_position_2025.rds",
  sep=""))

data_abun_2025 <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_abun_2025.rds",
  sep=""))

count_tot <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/HOLOTVSPM2025.rds",
  sep=""))


# creation of ppp for one haul
##From sf to ppp
stn <- 100
stn_position <- st_as_sf(data_position %>% filter(station==stn),
                         coords = c("X","Y"))
X <- as.ppp(st_coordinates(stn_position), st_bbox(stn_position))
plot(X)

# creation of object of class owin then a ppp
stn_abun <- data_abun_2025 %>% filter(station==stn)

win <- owin(xrange = c(0, 150), yrange = c(0, stn_abun$area[1]/1.5))
X2 <- ppp(st_coordinates(stn_position)[,1]*100,
          st_coordinates(stn_position)[,2], window = win)
plot(X2)

win <- owin(xrange = c(0, 1.5), yrange = c(0, stn_abun$area[1]/1.5))
X3 <- ppp(st_coordinates(stn_position)[,1],
          st_coordinates(stn_position)[,2], window = win)
plot(X3)


#############
#Intensity
#############

# test on ppp for one haul ####
intensity(X2)
intensity(X3)
sqrt(intensity(X2)/area(Window(X2)))
quadratcount(X2,3,461)
intensity(quadratcount(X2,3,46.1))
plot(intensity(quadratcount(X2,3,46.1), image=TRUE))

den <- density(X2, sigma=1)
persp.im(den)
contour.im(den)

b <- bw.ppl(X2)

den <- density(X2, sigma=b)
persp.im(den)
contour.im(den)

#####

#############
#Correlation
#############

# test on ppp for one haul ####

# Morisita index
miplot(X2)
miplot(X3)# reality but not enough point

# Fry plot
fryplot(X2)
fryplot(X3)

# K- and L- Riley functions
K3 <- Kest(X3)
L3 <- Lest(X3)

K2 <- Kest(X2)
L2 <- Lest(X2)

#plot(K,L,K2,L2)
fct_L3 <- as.function(Lest(X3))
x <- seq(0,0.37,0.0001)
ggplot()+
  geom_point(aes(x=x, y=x), shape = 3, color = '#999999')+
  geom_point(aes(x=x, y=fct_L3(x)), color = "black")+
  ylab("L3(X)")

# pair correlation function
g3 <- pcf(X3)
g2 <- pcf(X2)
plot(g3)
g3.1 <- pcf(X3, divisor = "d")
plot(g3.1)
g3.2 <- pcf(K3, method = "b")
plot(g3.2)
g3.3 <- pcf(K3, spar = 0.5)
plot(g3.3)

fct_G3.1 <- as.function(g3.1)
fct_G3.1(0.375)

# variance under CSR
varK3 <- Kest(X3, var.approx = TRUE)
plot(varK3)

# block bootstrap
Kvb3 <- varblock(X3,Kest,nx=3,ny=20)
plot(Kvb3)

# Loh's bootstrap
Kloh3 <- lohboot(X3, Kest)
plot(Kloh3)

# pointwise envelopes
E3 <- envelope(X3, Kest, nsim = 39, fix.n = TRUE)
plot(E3)
E3.1 <- envelope(X3, Kest, nsim = 99)
plot(E3.1)

# global envelopes
Eg3 <- envelope(X3, Kest, nsim = 19, rank=1, global = TRUE)
plot(Eg3)

# Non-graphical tests
alpha <- 1/(99+1)
## maximum absolute deviation
mad.test(X3, Lest, nsim=99,rmax=0.375, use.theo=TRUE)
## Diggle-Cressie-Loosmore-Ford
dclf.test(X3, Lest, nsim=99, rmax=0.375, use.theo=TRUE)#$p.value

# test in case of inhomogeneity
lambda3 <- density(X3, bw.ppl)
inhK3 <- Kinhom(X3, lambda3)
plot(inhK3)

Einh3 <- envelope(X3, Linhom, sigma=bw.ppl,
              simulate=expression(rpoispp(lambda3)),
              use.theory=TRUE, nsim=19, global=TRUE)
plot(Einh3)

# calculate difference between the Poisson function and the ppp
max(L3$r) - integrate(fct_G3.1, lower = 0, upper = max(L3$r))$value
integrate(function(x){x}, lower = 0, upper = max(L3$r))$value - 
  integrate(as.function(L3), lower = 0, upper = max(L3$r))$value

#####

# test all ppp ####
data_PPP_2025 <- data_abun_2025 %>% select(station,intensity) %>%
  mutate(L_area_dif = 0) %>%
  mutate(g_area_dif = 0) %>%
  mutate(MAD_test = 0) %>%
  mutate(DCLF_test = 0)

test <- c()
for (indice in 1:nrow(data_PPP_2025)){
  stn <- data_PPP_2025$station[indice]
  view_field <- 1.5 # view of the GoPro in meter
  stn_position <- st_as_sf(data_position %>% filter(station==stn),
                           coords = c("X","Y"))
  stn_abun <- data_abun_2025 %>% filter(station==stn)
  stn_count <- count_tot %>% filter(STN == as.character(stn))
  win <- owin(xrange = c(0, view_field),
              yrange = c(stn_count$start_distance[1], 
                      (stn_abun$area[1]/view_field) +
                        stn_count$start_distance[1] +#to adapt cucumber position
                        1.5))# to capt last that aren't in owin due to process
  #ratio for the first part of the video 380 pixel = 1m and video 540 pixels 
  #height so 1/380*540 = 1.421~1.5
  X <- ppp(st_coordinates(stn_position)[,1],
            st_coordinates(stn_position)[,2], window = win)
  
  if (nrow(stn_position)!=X$n){
    test <- c(test,stn)
  }# check if all the point are in the ppp for each station
  
  # calculate L for testing the Ripley's K function of the PPP
  L <- Lest(X) #max(L$r) define the interval of definition of the L function
  # calculate the pair correlation function g of the ppp
  g <- pcf(X, divisor = "d") #using d divisor for better feat at small r
  # calculate difference between the Poisson function and the ppp
  # the reference for a Poisson ppp is L(r)=r
  L_dif <- integrate(function(x){x}, 
                     lower = 0, upper = max(L$r))$value - 
    integrate(as.function(L), lower = 0, upper = max(L$r),
              subdivisions=2000)$value # error due to not enougth subdivision
  # the reference for a Poisson ppp is g(r)=1
  g_dif <- max(L$r) - integrate(as.function(g), 
                                 lower = 0, upper = max(L$r))$value
  #g_dif negative = cluster and g_dif positive = inhibition
  
  # Non-graphical tests
  alpha <- 1/(99+1) # so the threshold value is 0.01
  ## maximum absolute deviation
  mad <- mad.test(X, Lest, nsim=99,rmax=max(L$r), use.theo=TRUE)$p.value
  ## Diggle-Cressie-Loosmore-Ford
  dclf <- dclf.test(X, Lest, nsim=99, rmax=max(L$r), use.theo=TRUE)$p.value
  
  data_PPP_2025$L_area_dif[indice] <- L_dif
  data_PPP_2025$g_area_dif [indice] <- g_dif
  data_PPP_2025$MAD_test [indice] <- mad
  data_PPP_2025$DCLF_test [indice] <- dclf
}
rm(indice,stn,stn_position,stn_abun,win,X,L,g,L_dif,g_dif,mad,dclf)

#####
saveRDS(data_PPP_2025,
        paste(here(),
              "/via3_data_exploration/Data/processed/data_PPP_2025.rds",
              sep=""))

data_PPP_2025 <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_PPP_2025.rds",
  sep=""))
row.names(data_PPP_2025) <- data_PPP_2025$station

# visualize PPP ####
stn_test <- data_position %>% filter(station==stn)
ggplot()+
  geom_point(data = stn_test, aes(x=X,y=Y))

ggplot(data = data_position)+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 4,scales="free_y")+
  ylab("")

#acp to try to create category of ppp
x <- as.data.frame(data_PPP_2025) %>% 
  select(intensity, L_area_dif, g_area_dif, MAD_test,
         DCLF_test)
row.names(x) <- data_PPP_2025$station
acp1 <- prcomp(x,  scale = F)
acp1
plot(acp1)

acp3 = FactoMineR::PCA(x,scale.unit=F, ncp=5, graph=F)
acp3
fviz_pca_ind(acp3,
             addEllipses = TRUE,      # Add ellipses for categories
             repel = TRUE,            # Avoid label overlap
             title = "PCA with Categories")
fviz_pca_biplot(acp3,
                addEllipses = T,      # Add ellipses for categories
                repel = F,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")

# Avec des données ACP issues de prcomp() de R par défaut
#axeX <- acp1$x[,1] ; axeY <- acp1$x[,2] 
#
## Avec des données ACP issues de PCA() de FactoMineR
#axeX <- acp3$ind$coord[,1] ; axeY <- acp3$ind$coord[,2]
#
#plot(axeX,axeY,pch=16);grid()
#text(axeX,axeY,rownames(x))
ggplot(data = data_position %>% filter(station %in% c(139,169,201,135,212,185)))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y")+
  ylab("")
ggplot(data = data_position %>% filter(station %in% c(157,123,183,
                                                      127,155,174,104)))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow=1)+
  ylab("")

#acp to try to create category of ppp without the intensity
x2 <- as.data.frame(data_PPP_2025) %>% 
  select(L_area_dif, g_area_dif, MAD_test,
         DCLF_test)
row.names(x2) <- data_PPP_2025$station
acp1.2 <- prcomp(x2,  scale = F)
acp1.2
plot(acp1.2)

acp3.2 = FactoMineR::PCA(x2,scale.unit=F, ncp=5, graph=F)
acp3.2
fviz_pca_ind(acp3.2,
             addEllipses = TRUE,      # Add ellipses for categories
             repel = TRUE,            # Avoid label overlap
             title = "PCA with Categories")
fviz_pca_biplot(acp3.2,
                addEllipses = T,      # Add ellipses for categories
                repel = F,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")

#groups of high MAD_test and DCLF_test
ggplot(data = data_position %>% filter(station %in% c(161,143,190,
                                                      106,222,144,
                                                      129,210,220,
                                                      100,193)))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow=2)+
  ylab("")
#groups of low g_area_dif
ggplot(data = data_position %>% filter(station %in% c(157,123,183,
                                                      155,174,104)))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow=1)+
  ylab("")
# stn of high tests and low g_area_dif
ggplot(data = data_position %>% filter(station==127))+
  geom_point(aes(x=X,y=Y))+
  ylab("")

#####

#############
#Spacing 
#############

# test on ppp for one haul ####
stn = 102
stn_position <- st_as_sf(data_position %>% filter(station==stn),
                         coords = c("X","Y"))

# creation of object of class owin then a ppp
stn_abun <- data_abun_2025 %>% filter(station==stn)

win <- owin(xrange = c(0, 1.5), yrange = c(0, stn_abun$area[1]/1.5))
X <- ppp(st_coordinates(stn_position)[,1],
          st_coordinates(stn_position)[,2], window = win)
plot(X)
ggplot(data = data_position %>% filter(station==stn))+
  geom_point(aes(x=X,y=Y))+
  ylab("")

#pairwise distances
M <- pairdist(X)
summary(M)[,1:3]
#nearest-neighbour distances
v <- nndist(X)
v[1:3]
#empty-space distance
Z <- distmap(X)
plot(Z)
U <- runifpoint(3, Window(X))
Z[U]

#Clark-Evans index
clarkevans(X)
clarkevans.test(X,correction="donnelly",
                alternative="clustered")
###random:R=1, cluster:R<1, regularity:R>1 biased by inhomogeneity
#Hopkins-Skellam index
hopskel(X)
hopskel.test(X, alternative="clustered")
###random:A=1, cluster:A<1, regularity:A>1

#nearest-neighbour function G
Gs <- Gest(X)
Gs
plot(Gs)
### G under ref=regular and G above ref=cluster
#empty-space function F
Fs <- Fest(X)
Fs
plot(Fs)
### F above ref=regular and F under ref=cluster p266
#hazard rate
plot(Fest(X), cbind(hazard, theohaz) ~ r)
ggplot(as.data.frame(Fs))+
  geom_point(aes(x=r, y=theohaz), shape = 3, color = '#999999')+
  geom_point(aes(x=r, y=hazard), color = "black")+
  xlim(c(0,1.5))+
  ylim(c(0,2))+
  ylab("h(r)")

## Extract the values of F
#r <- Fs$r
#Fvals <- Fs$theo
## Compute the hazard function
#hazard <- -diff(log(1 - Fvals)) / diff(r)
## The 'r' values for the hazard function are midpoints
#r_mid <- (r[-1] + r[-length(r)]) / 2
## Combine into a data frame or 'fv' object
#hazard_fv <- fv(x = data.frame(r_mid,hazard),"r_mid",quote(A(r_mid)),"hazard",
#                hazard ~ r_mid)

#the diagnostic Dmax
DD <- eval.fv(Gs-Fs)
plot(DD)
DigDif <- function(X, ..., r=NULL) {
  FX <- Fest(X, ..., r=r)
  GX <- Gest(X, ..., r=FX$r)
  eval.fv(GX-FX)
}
DE <- envelope(X, DigDif, nsim=19, fix.n=TRUE, global=TRUE)
plot(DE)

#the J-function
### Poisson:J(r)=1, regular:J(r)>1, cluster:J(r)<1
J <- (1-Gs)/(1-Fs)
plot(J)

#all functions
#plot(allstats(X)) error
#####

# test all ppp ####
data_PPP_2025 <- data_PPP_2025 %>% 
  mutate(hopkins_index = 0) %>%
  mutate(J_area_dif = 0) %>%
  mutate(J_area_percent = 0)

list_PPP <- list()
for (indice in 1:nrow(data_PPP_2025)){
  stn <- data_PPP_2025$station[indice]
  view_field <- 1.5 # view of the GoPro in meter
  stn_position <- st_as_sf(data_position %>% filter(station==stn),
                           coords = c("X","Y"))
  stn_abun <- data_abun_2025 %>% filter(station==stn)
  stn_count <- count_tot %>% filter(STN == as.character(stn))
  win <- owin(xrange = c(0, view_field),
              yrange = c(stn_count$start_distance[1], 
                         (stn_abun$area[1]/view_field) +
                           stn_count$start_distance[1] +#to adapt cucumber position
                           1.5))# to capt last that aren't in owin due to process
  #ratio for the first part of the video 380 pixel = 1m and video 540 pixels 
  #height so 1/380*540 = 1.421~1.5
  X <- ppp(st_coordinates(stn_position)[,1],
           st_coordinates(stn_position)[,2], window = win)
  
  list_PPP[[indice]] <- X

  #Hopkins-Skellam index
  hopkins_index <- hopskel(X)
  ###random:A=1, cluster:A<1, regularity:A>1
  
  #nearest-neighbour function G
  Gs <- Gest(X)
  Gs <- Gs %>% filter(!is.na(rs))
  ### G under ref=regular and G above ref=cluster
  #empty-space function F
  Fs <- Fest(X)
  #the J-function
  ### Poisson:J(r)=1, regular:J(r)>1, cluster:J(r)<1
  # Define a finely spaced sequence of r values because error if > 0.00293
  r <- seq(0, max(Gs$r), by = 0.001)  # Adjust the range and step size as needed
  J <- Jest(X, r=r)
  #problem with creation of NA that error the integrate
  J <- J %>% filter(!is.na(rs))
  # the reference for a Poisson ppp is J(r)=1
  J_dif <- max(J$r) - integrate(as.function(J), 
                                lower = 0, upper = max(J$r),
                                subdivisions=2000)$value
  #J_dif negative = regular and J_dif positive = cluster
  J_ratio <- integrate(as.function(J), 
                                  lower = 0, upper = max(J$r),
                                  subdivisions=2000)$value / max(J$r)  
  #J_ratio > 1  = regular and J_ratio < 1 = cluster
    
  data_PPP_2025$hopkins_index[indice] <- hopkins_index
  data_PPP_2025$J_area_dif [indice] <- J_dif
  data_PPP_2025$J_area_percent [indice] <- J_ratio
}
rm(indice,stn,stn_position,stn_abun,stn_count,win,X,J,Fs,Gs,J_dif,hopkins_index)

#####
saveRDS(data_PPP_2025,
        paste(here(),
              "/via3_data_exploration/Data/processed/data_PPP_2025.rds",
              sep=""))

data_PPP_2025 <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/data_PPP_2025.rds",
  sep=""))

#acp to try to create category of ppp ####
x <- as.data.frame(data_PPP_2025) %>% 
  select(hopkins_index,J_area_dif,J_area_percent)
row.names(x) <- data_PPP_2025$station
acp1 <- prcomp(x,  scale = F)
acp1
plot(acp1)

acp3 = FactoMineR::PCA(x,scale.unit=F, ncp=5, graph=F)
acp3
fviz_pca_ind(acp3,
             addEllipses = TRUE,      # Add ellipses for categories
             repel = TRUE,            # Avoid label overlap
             title = "PCA with Categories")
fviz_pca_biplot(acp3,
                addEllipses = T,      # Add ellipses for categories
                repel = F,            # Avoid label overlap
                title = "PCA Biplot")+
  theme() +
  labs(title = "Customized PCA Biplot")

coord_dim12 <- as.data.frame(acp3[["ind"]][["coord"]][,c(1,2)])

group1 <- row.names(coord_dim12 %>% filter(Dim.1<0.5) %>% filter(Dim.2<(-0.15)))
# high J_area_dif, low J_area_percent, hopkins<1 = cluster
group2 <- row.names(coord_dim12 %>% filter(Dim.2>(0.15)))
# 
group3 <- row.names(coord_dim12 %>% filter(Dim.1>1))
#high hpkins index = regularity

ggplot(data = data_position %>% filter(station %in% group1))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow = 3)+
  ylab("")

ggplot(data = data_position %>% filter(station %in% group2))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow = 3)+
  ylab("")

ggplot(data = data_position %>% filter(station %in% group3))+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station,scales="free_y",nrow = 1)+
  ylab("")


#####


#############
#Summary of all PPP
#############
x <- as.data.frame(data_PPP_2025) %>% 
  select(intensity, L_area_dif, g_area_dif, MAD_test,
         DCLF_test,hopkins_index,J_area_dif,J_area_percent)
row.names(x) <- data_PPP_2025$station
acp1 <- prcomp(x,  scale = F)
acp1
plot(acp1)

acp3 = FactoMineR::PCA(x,scale.unit=F, ncp=5, graph=F)
acp3
fviz_pca_ind(acp3,
             addEllipses = TRUE,      # Add ellipses for categories
             repel = TRUE,            # Avoid label overlap
             title = "PCA with Categories")
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

# ACP without the intensity but just test of spacing and correlation ####
x <- as.data.frame(data_PPP_2025) %>% 
  select(L_area_dif, g_area_dif, MAD_test,
         DCLF_test,hopkins_index,J_area_dif,J_area_percent)
row.names(x) <- data_PPP_2025$station
acp1 <- prcomp(x,  scale = F)
acp1
plot(acp1)

acp3 = FactoMineR::PCA(x,scale.unit=F, ncp=5, graph=F)
acp3
fviz_pca_ind(acp3,
             addEllipses = TRUE,      # Add ellipses for categories
             repel = TRUE,            # Avoid label overlap
             title = "PCA with Categories")
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

# plot the PPP ####
ggplot(data = data_position)+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~station, nrow = 4,scales="free_y")+
  ylab("")

data_test <- data_position %>% mutate(PPP_intensity=0)
for (indice in 1:nrow(data_test)){
  stn <- data_test$station[indice]
  if (stn %in% as.character(data_PPP_2025$station)){
    data_test$PPP_intensity[indice]<- data_PPP_2025[stn,]$intensity
  }
}
rm(indice,stn)

ggplot(data = data_test)+
  geom_point(aes(x=X,y=Y))+
  facet_wrap(~PPP_intensity*station, nrow = 4,scales="free_y")+
  ylab("")

#####

saveRDS(list_PPP,
        paste(here(),
              "/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
              sep=""))

list_PPP <- readRDS(paste(
  here(),"/via3_data_exploration/Data/processed/list_PPP_tuyau_2025.rds",
  sep=""))
