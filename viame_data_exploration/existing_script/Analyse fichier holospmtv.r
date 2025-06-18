#############################################################
##
## Script d'analyses des données de la campagne HOLOSPMTV
## Elaboration d'un TAC à partir des données de la campagne
## Structure du document
##    Partie A - Importation des données
##    Partie B - 1eres inferences stats
##    Partie C - Estimation des Variances d'échantillonnage
##    Partie D - Regle de contrôle de gestion (HCR)
##    Partie E - Analyse en inter annuel
##
##  Joël Vigneau (HMMN/LRHPB) - Juillet 2022
##
#############################################################


library(readxl)
library(dplyr)
library(ggplot2)

##############
##############     Partie A - Importation des données
##############

# Choisir un chemin où se trouvent les tables excel à vérifier (le nom du pays doit être dans le nom du fichier Excel) :
# Paramétrage
# wd_path <- 'C:/_PROGRAMME/NAFO/concombre'
wd_path <- './'
AN      <- 22
####

myxls <- list.files(wd_path)[grepl('.xls',list.files(wd_path))]
myxls <- myxls[substring(myxls,1,1) %in% 'H']

(myxls <- myxls[grepl(AN, myxls)])


# myxls <- myxls[2]


tab  <- read_xlsx(path = paste(wd_path,myxls,sep='/'), sheet = 'Feuil1')
rangeRows <- range(which(!(is.na(tab[,1])))) +1
param <- read_xlsx(path = paste(wd_path,myxls,sep='/'), range = 'A1:Z10')
clns <- which(apply(!apply(param,2, is.na),2,sum) >0)
info <- apply(param,2,is.na)[,clns]
dimnames(info) <- list(2:10,gsub('[nomunitevalré.]','',dimnames(info)[[2]]))
info <- param[as.numeric(dimnames(info)[[1]]), as.numeric(dimnames(info)[[2]])]
param <- list(lTransect=as.numeric(info[which(grepl('transect',c(info[,1])[[1]])),2]),
              champVision=as.numeric(info[which(grepl('vision',c(info[,1])[[1]])),2]),
              pdsAdultes=as.numeric(info[which(grepl('adult',c(info[,1])[[1]])),2]),
              pdsJuveniles=as.numeric(info[which(grepl('juvenile',c(info[,1])[[1]])),2]),
              tauxExp=as.numeric(info[which(grepl('Quota',c(info[,ncol(info)-2])[[1]]))[1],ncol(info)-1]),
              surfCarre=as.numeric(info[which(grepl('carré stat',c(info[,1])[[1]])),2]))

tabl  <- read_xlsx(path = paste(wd_path,myxls,sep='/'), range = cell_rows(rangeRows), sheet = 'Feuil1')
tabl <- as.data.frame(tabl)

table(tabl$STN)
tabl$area <- ''
tabl$area[substring(tabl$STN,nchar(tabl$STN),nchar(tabl$STN)) %in% 'B'] <- 'Boite'
tabl$area[substring(tabl$STN,nchar(tabl$STN),nchar(tabl$STN)) %in% 'C'] <- 'Ouest'
tabl$area[substring(tabl$STN,nchar(tabl$STN),nchar(tabl$STN)) %in% 'N'] <- 'Tuyau Nord'
tabl$area[substring(tabl$STN,nchar(tabl$STN),nchar(tabl$STN)) %in% 'S'] <- 'Tuyau Sud'

if (AN %in% 21) tabl <- tabl[!(tabl$STN %in% c('169S','171S','173S','179S')),]  #Suppression des traines ajoutées par les pros (en 2021)

##############
##############    Partie B - 1eres inferences stats
##############

## nbJuvStation = nombre de juvéniles comptés pour chaque station
nbJuvStation <- (tabl[,which(grepl('Nbre juv',names(tabl)))] / tabl$'taux echant')
nbJuvStation[is.na(nbJuvStation)] <- 0
## nbAdultStation = nombre d'adultes comptés pour chaque station
nbAdultStation <- (tabl[,which(grepl('Nbre adult',names(tabl)))] / tabl$'taux echant')
nbAdultStation[is.na(nbAdultStation)] <- 0
control <- (nbJuvStation + nbAdultStation)-tabl$'Nombre total' == 0
## intervalTemps = durée de la vidéo relative au comptage
intervalTemps <- tabl$Stop - tabl$Start
## surfCompte = surface relative au comptage (postulat = traîne à 3 knots)
vitesseNavire <- 1.4
surfCompte <- as.numeric(intervalTemps*vitesseNavire*1852/60)

## abJuv et abAdult = nombre de juvéniles et adultes rapporté à la surface totale
nbCarreStat=2*nrow(tabl)
abJuv <- param$surfCarre * (sum(nbJuvStation, na.rm=TRUE)/sum(surfCompte,na.rm=TRUE)) * nbCarreStat
abAdult <- param$surfCarre * (sum(nbAdultStation, na.rm=TRUE)/sum(surfCompte, na.rm=TRUE)) * nbCarreStat
## bioAdult = biomasse d'adultes en tonnes rapporté à la surface totale
bioAdult <- param$pdsAdultes *abAdult/1000

TAC <- bioAdult*param$tauxExp
TAC

## Analyse par zone
nbCarreStat=2*table(tabl$area)
## abJuv et abAdult = nombre de juvéniles et adultes rapporté à la surface totale
juvRatio <- tapply(nbJuvStation,tabl$area,sum,na.rm=TRUE)/tapply(surfCompte,tabl$area,sum,na.rm=TRUE)
adultRatio <- tapply(nbAdultStation,tabl$area,sum,na.rm=TRUE)/tapply(surfCompte,tabl$area,sum,na.rm=TRUE)
abJuv <- juvRatio * nbCarreStat * as.numeric(param$surfCarre)
abAdult <-  adultRatio * nbCarreStat * as.numeric(param$surfCarre)
## bioAdult = biomasse d'adultes en tonnes rapporté à la surface totale

## RESULTAT BIOMASSE --------------------------------------
bioAdult <- as.numeric(param$pdsAdultes) *(abAdult/1000)
## --------------------------------------------------------


TAC <- bioAdult*as.numeric(param$tauxExp)
TAC


bioAdult
# 2023 Tuyau Nord  Tuyau Sud 
# Tuyau Nord  Tuyau Sud
#    23271.7   115489.7
  
#    Ouest Tuyau Nord  Tuyau Sud
# 23506.64   22255.75  199354.94

##############
##############    Partie C - Estimation des Variances d'échantillonnage
##############

### Variances par bootstrap
# N = nombre de réplicats à définir
N <- 10000

nZones <- length(table(tabl$area))
vj_a <- numeric(nZones)
va_a <- numeric(nZones)
abJu <- list()
abAd <- list()

for (i in 1:nZones) {
  abjuv <- numeric(N)
  abAdult <- numeric(N)
  for (j in 1:N) {
    y <- nbJuvStation[tabl$area %in% names(table(tabl$area))[i]]
    z <- nbAdultStation[tabl$area %in% names(table(tabl$area))[i]]
    x <- surfCompte[tabl$area %in% names(table(tabl$area))[i]]
    lesNonNA <- which(!(is.na(x)))
    y <- y[lesNonNA]; z = z[lesNonNA]; x <- x[lesNonNA]
    idx <- sample(1:length(y), length(y), replace=TRUE)
    yy <- y[idx]
    zz <- z[idx]
    xx <- x[idx]
    juvRatio <- sum(yy)/sum(xx)
    adultRatio <- sum(zz)/sum(xx)
    abJuv[j] <- juvRatio * nbCarreStat[i] * as.numeric(param$surfCarre)
    abAdult[j] <-  adultRatio * nbCarreStat[i] * as.numeric(param$surfCarre)
  }
  abJu[[i]] <- abJuv
  names(abJu[i]) <- names(nbCarreStat[i])
  abAd[[i]] <- abAdult
  names(abAdult[i]) <- names(nbCarreStat[i])
}

ic_biomasse <- array(dim=c(nZones,7))
for (i in 1:nZones) {
  pipo <- abAd[[i]][order(abAd[[i]])]
  ic_biomasse[i,] <- c(pipo[c(N*0.025,N*0.050,N*0.100)]* as.numeric(param$pdsAdultes)/1000, bioAdult[i], pipo[c(N*0.900,N*0.950,N*0.975)] * as.numeric(param$pdsAdultes)/1000)
}

dimnames(ic_biomasse) <- list(c(names(table(tabl$area))),c('IC 95% min','IC 90% min','IC 80% min','Moyenne','IC 80% max','IC 90% max','IC 95% max'))
## RESULTAT BIOMASSE ADULTE AVEC IC BOOTSTRAP-------------------------------------
floor(ic_biomasse)
## ------------------------------------------------------------------------

#  Biomasse juvéniles à revoir....
#ic_biomasse <- array(dim=c(nZones,7))
# for (i in 1:nZones) {
#   pipo <- abJu[[i]][order(abJu[[i]])]
#   ic_biomasse[i,] <- c(pipo[c(N*0.025,N*0.050,N*0.100)]* as.numeric(param$pdsJuveniles)/1000, bioAdult[i], pipo[c(N*0.900,N*0.950,N*0.975)] * as.numeric(param$pdsAdultes) /1000)
# }
#
# dimnames(ic_biomasse) <- list(c(names(table(tabl$area))),c('IC 95% min','IC 90% min','IC 80% min','Moyenne','IC 80% max','IC 90% max','IC 95% max'))
# ## RESULTAT BIOMASSE JUVENILE AVEC IC BOOTSTRAP-------------------------------------
# floor(ic_biomasse)
# ## ------------------------------------------------------------------------

TAC <- ic_biomasse*as.numeric(param$tauxExp)
floor(TAC)

if (nZones == 2) {
  dfAd <- data.frame(zone=rep(names(table(tabl$area)),each=N), replicat=c(abAd[[1]],abAd[[2]]))
  dfJu <- data.frame(zone=rep(names(table(tabl$area)),each=N), replicat=c(abJu[[1]],abJu[[2]]))
} else {
  dfAd <- data.frame(zone=rep(names(table(tabl$area)),each=N), replicat=c(abAd[[1]],abAd[[2]],abAd[[3]]))
  dfJu <- data.frame(zone=rep(names(table(tabl$area)),each=N), replicat=c(abJu[[1]],abJu[[2]],abJu[[3]]))
}

## FIGURE : ABONDANCE ADULTES
dfAd2 <- dfAd %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_path,'/Abondance ',AN,' adultes.jpeg'), width=1000, height=400)
ggplot(data=dfAd, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  geom_vline(data=dfAd2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='Abondance (réplicats bootstrap)', y='Fréquence', title=paste0("Abondance adultes holothuries SPM - campagne 20",AN)) +
  facet_wrap(~zone, scales='free')
dev.off()

# FIGURE : BIOMASSE ADULTES
dfBio <- dfAd
dfBio$replicat <- dfAd$replicat*.4*10^-6
dfBio2 <- dfBio %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_path,'/Biomasse ',AN,' adultes.jpeg'), width=1000, height=400)
ggplot(data=dfBio, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  geom_vline(data=dfBio2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='Biomasse (kt)', y='Fréquence', title=paste0("Biomasse adultes holothuries SPM - campagne 20",AN)) +
  facet_wrap(~zone, scales='free')
dev.off()

# FIGURE : BIOMASSE JUVENILES
dfJu2 <- dfJu %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_path,'/Abondance ',AN,' juvéniles.jpeg'), width=1000, height=400)
ggplot(data=dfJu, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  geom_vline(data=dfJu2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='Abondance (réplicats bootstrap)', y='Fréquence', title=paste0("Abondance juvéniles holothuries SPM - campagne 20",AN)) +
  facet_wrap(~zone, scales='free')
dev.off()

dfTot <- rbind.data.frame(data.frame(dfAd, stade='Adultes'), data.frame(dfJu, stade='Juvéniles'))
dfTot$stade <- factor(dfTot$stade, ordered=TRUE, levels = c('Juvéniles','Adultes'))

# FIGURE : BOXPLOT ABONDANCE
jpeg(paste0(wd_path,'/boxplot abondance ',AN,'.jpeg'), width=1000, height=400)
ggplot(data=dfTot, aes(x=stade, y=replicat, fill=stade)) +
  geom_boxplot() +
  labs(x='', y='Nombre total estimé', title=paste0("Abondance holothuries SPM - campagne 20",AN)) +
  facet_wrap(~zone, scales='free') +
  scale_fill_brewer()
dev.off()

save(dfAd,dfAd2, dfJu, dfJu2, dfTot, file=paste0(wd_path,'Boot 20',AN,'.Rdata'))

##############
##############    Partie D - Regle de contrôle de gestion (HCR)
##############


Tmax <- 0.015  #Taux d'exploitation donnant une valeur de biomasse moyenne = Vmax

dfAd22 <- data.frame(annee=2022,dfAd); dfJu22 <- data.frame(annee=2022,dfJu); dfTot22 <- data.frame(annee=2022,dfTot); statAd22 <- dfAd2; statJu22 <- dfJu2
dfBmax <- dfAd22[dfAd22$zone %in% c('Tuyau Sud','Ouest'),]
dfBmax$zone <- as.character(dfBmax$zone)
dfBmax$replicat <- Tmax*dfBmax$replicat*.4*10^-3
tab.Bmax <- dfBmax %>% group_by(zone) %>% summarise(q01 = quantile(replicat, 0.01), q05 = quantile(replicat, 0.05), q35 = quantile(replicat, 0.35),
                                                   moyenne=mean(replicat))
tab.Bmax
Vmax <- tab.Bmax[,c('zone','moyenne')]
Vmax

#QUOTA POSSIBLE (Valeur donnant 95% de chance d'être inférieur à valeur de biomaase (Vmax) correspondant à un taux d'exploitation = 1.5% (Tmax) ..................
tab.Bmax$q05
#---------------------------------

# FIGURE ILLUSTRANT LE PROPOS
jpeg(paste0(wd_path,'/TAC Campagne ',AN,' adultes.jpeg'), width=1000, height=400)
ggplot(data=dfBmax[dfBmax$zone %in% c('Tuyau Sud','Ouest'),], aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='blue') +
  geom_vline(data=tab.Bmax, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=tab.Bmax, aes(xintercept=c(moyenne)), color="red", linetype="dashed") +
  labs(x='Captures (t)', y='Fréquence', title=paste0("campagne HOLOSPMTV 20",AN)) +
  facet_wrap(~zone, scales='free')
dev.off()


## Alternative : distribution centrée sur TAC (non gardée in fine car jugée trop complexe)
## idée : rechercher un Teff dont la queue de ditribution (90eme, 95eme ou 99eme percentile) = Vmax

# QUOTA POSSIBLE : meme idée que ci-dessus mais de l'autre côté; recherche d'un taux effectif (Teff) dont la biomasse correspondante a une distribution < Vmax à 95%
Teff <- 0.00997  #Tuyau Sud
dfBio$replicat <- dfAd22$replicat*.4*10^-3*Teff
tab.Beff <- dfBio %>% group_by(zone) %>% summarise(q35=quantile(replicat,0.35), mediane=quantile(replicat,0.5),moyenne=mean(replicat),
                                               q90=quantile(replicat,0.9), q95=quantile(replicat,.95), q99=quantile(replicat,.99))
tab.Beff
Teff
Vmax

# QUOTAS POSSIBLES POUR TUYAU SUD
tab.Beff[3,'moyenne']
# -------------------------------


zz <- 'Tuyau Sud'
jpeg(paste0(wd_path,'/TAC tuyau Sud.jpeg'), width=300, height=250)
ggplot(data=dfBio[dfBio$zone %in% zz,], aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='blue') +
  geom_vline(data=tab.Beff[tab.Beff$zone %in% zz,], aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  geom_vline(data=tab.Bmax, aes(xintercept=c(moyenne)), color="red", linetype="dashed") +
  labs(x='TAC (réplicats bootstrap)', y='Fréquence', title=paste0("Campagne 20",AN,' - Zone ',zz))
dev.off()



# Meme construction pour la zone Ouest
Teff <- 0.00735  #Ouest
dfBio$replicat <- dfAd22$replicat*.4*10^-3*Teff
tab.Beff <- dfBio %>% group_by(zone) %>% summarise(q35=quantile(replicat,0.35), mediane=quantile(replicat,0.5),moyenne=mean(replicat),
                                                   q90=quantile(replicat,0.9),q95=quantile(replicat,.95), q99=quantile(replicat,.99))
tab.Beff
Teff
Vmax

zz <- 'Ouest'
jpeg(paste0(wd_path,'/TAC Ouest.jpeg'), width=300, height=250)
ggplot(data=dfBio[dfBio$zone %in% zz,], aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='blue') +
  geom_vline(data=tab.Beff[tab.Beff$zone %in% zz,], aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  geom_vline(data=tab.Bmax, aes(xintercept=c(moyenne)), color="red", linetype="dashed") +
  labs(x='TAC (réplicats bootstrap)', y='Fréquence', title=paste0("Campagne 20",AN,' - Zone ',zz))
dev.off()


##############
##############   Partie E - Analyse en inter annuel
##############

load(paste0(wd_path,'/Boot 2021.Rdata'))
dfAd21 <- data.frame(annee=2021,dfAd); dfJu21 <- data.frame(annee=2021,dfJu); dfTot21 <- data.frame(annee=2021,dfTot); statAd21 <- dfAd2; statJu21 <- dfJu2
load(paste0(wd_path,'/Boot 2022.Rdata'))
dfAd22 <- data.frame(annee=2022,dfAd); dfJu22 <- data.frame(annee=2022,dfJu); dfTot22 <- data.frame(annee=2022,dfTot); statAd22 <- dfAd2; statJu22 <- dfJu2

dfTotn <- rbind.data.frame(dfTot21, dfTot22)

jpeg(paste0(wd_path,'/boxplot abondance.jpeg'), width=1000, height=400)
ggplot(data=dfTotn, aes(x=as.character(annee), y=replicat, fill=stade)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='', y='Nombre total estimé', title="Abondance holothuries SPM") +
  facet_wrap(~zone, scales='free') +
  scale_fill_brewer()
dev.off()

jpeg(paste0(wd_path,'/boxplot abondance2.jpeg'), width=1000, height=400)
ggplot(data=dfTotn, aes(x=as.character(annee), y=replicat, fill=stade)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='', y='Nombre total estimé', title="Abondance holothuries SPM") +
  facet_wrap(~stade+zone, scales='free') +
  scale_fill_brewer()
dev.off()


