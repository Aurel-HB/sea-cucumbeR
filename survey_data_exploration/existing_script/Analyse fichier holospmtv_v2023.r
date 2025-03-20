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
##  Joël Vigneau (HMMN/LRHPB) - Juillet 2023
##
#############################################################


library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpol)
library(viridis)



##############
##############     Partie A - Importation des données
##############

# Choisir un chemin où se trouvent les tables excel à vérifier (le nom du pays doit être dans le nom du fichier Excel) :
# Paramétrage
wd_path <- 'C:/_PROGRAMME/NAFO/concombre'
wd_figures <- 'C:/_PROGRAMME/NAFO/concombre/figs'
AN      <- 23

fig.width <- 2000
fig.height <- 2000
fig.res <- 300

#Taille des Zones
Tzone <- c(60, 60, 50, 15)  #Arbitraire à ce stade!
names(Tzone) <- c('Tuyau Nord','Tuyau Sud','Boîte à pétoncles','Ouest Miquelon')
####
strZone <- read_xlsx(path=paste(wd_path,'ZEE SPM strate Zone.xlsx',sep='/'))

myxls <- list.files(wd_path)[grepl('.xls',list.files(wd_path))]
myxls <- myxls[substring(myxls,1,1) %in% 'H']
myxls <- myxls[grepl(AN, myxls)]

tab  <- read_xlsx(path = paste(wd_path,myxls,sep='/'), sheet = 'Feuil1')
#Prise d'info de la feuille principale = paramètres (param) et matrice de données (tabl)
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
tabl <- merge(tabl, strZone, by.x='STN', by.y='Strate', all.x=TRUE)

if (AN %in% 21) tabl <- tabl[!(tabl$STN %in% c('169','171','173','179')),]  #Suppression des traines ajoutées par les pros (en 2021)
if (sum(is.na(tabl$STN)>0)) tabl <- tabl[!is.na(tabl$STN),]

##############
##############    Partie B - Initialisation des paramètres nécessaires aux estimations
##############

## taux echant = 1 à partir de 2022 car toutes les vidéos lues et interprétées entièrement
if (AN %in% 21) {nbJuvStation <- (tabl[,which(grepl('Nbre juv',names(tabl)))] / tabl$'taux echant')
	} else {nbJuvStation <- tabl[,which(grepl('Nbre juv',names(tabl)))]}
nbJuvStation[is.na(nbJuvStation)] <- 0
## nbAdultStation = nombre d'adultes comptés pour chaque station
if (AN %in% 21) {nbAdultStation <- (tabl[,which(grepl('Nbre adult',names(tabl)))] / tabl$'taux echant')
    } else {nbAdultStation <- tabl[,which(grepl('Nbre adult',names(tabl)))]}
nbAdultStation[is.na(nbAdultStation)] <- 0
control <- (nbJuvStation + nbAdultStation)-tabl$'Nombre total' == 0
## intervalTemps = durée de la vidéo relative au comptage
intervalTemps <- tabl$Stop - tabl$Start
## surfCompte = surface relative au comptage (en m2)
vitesseNavire <- 1.4
surfCompte <- as.numeric(intervalTemps*param$champVision*vitesseNavire*1852/60)
control

#Exploitation des infos de la feuille Excel
tt <- tapply(tabl$'biomasse carré stat (tonne)', tabl$Zone, mean, na.rm=TRUE)
tt*Tzone[names(tt)]

##############
##############    Partie C - &eres inférences statistiques
##############

## Analyse par zone
## abJuv et abAdult = nombre de juvéniles et adultes rapporté à la surface totale
## Ratio = somme(x)/somme(y) avec x = nombre d'individus comptés et y = surface balayée au comptage
juvRatio <- tapply(nbJuvStation,tabl$Zone,sum,na.rm=TRUE)/tapply(surfCompte,tabl$Zone,sum,na.rm=TRUE)
adultRatio <- tapply(nbAdultStation,tabl$Zone,sum,na.rm=TRUE)/tapply(surfCompte,tabl$Zone,sum,na.rm=TRUE)

## RESULTAT ABONDANCE --------------------------------------
## abondance en nombre d'individus rapporté à la surface totale de chaque zone
abJuv <- juvRatio * Tzone[names(juvRatio)] * as.numeric(param$surfCarre)
abAdult <-  adultRatio * Tzone[names(adultRatio)] * as.numeric(param$surfCarre)

## RESULTAT BIOMASSE --------------------------------------
## biomasse en tonnes rapporté à la surface totale de chaque zone
bioAdult <- (abAdult/1000) * param$pdsAdultes
bioJuv   <- (abJuv/1000)   * param$pdsJuveniles
## --------------------------------------------------------

myDF <- data.frame(station=tabl$STN, zone=tabl$Zone, type='Juveniles', nbInd=nbJuvStation)
myDF2 <- data.frame(station=tabl$STN, zone=tabl$Zone, type='Adults', nbInd=nbAdultStation)
myDF <- rbind.data.frame(myDF, myDF2)
techDF <- data.frame(station=tabl$STN,zone=tabl$Zone, nbInd=tapply(myDF$nbInd,myDF$station,max), taux=tabl$'taux echant')
techDF <- techDF[techDF$taux < 1,]
techDF <- techDF[!is.na(techDF$taux),]

## FIGURE : ABONDANCE ADULTES PAR STATION
jpeg(paste0(wd_figures,'/Stations ',AN,' abondance.jpeg'), width=fig.width, height=fig.height, res=0.5*fig.res)
ggplot(data=myDF, aes(x=station, y = nbInd)) +
  geom_point(aes(col=type, shape=type), size=3, alpha=0.5) +
  geom_text(data=techDF, aes(label=taux),nudge_y=1000, cex=2.5) +
  facet_wrap(~zone, scales='free', ncol=1) +
  guides(x = guide_axis(angle = 60))+
  labs(x='station', y='Nombre individus', title=paste0('Holothuries SPM - Campagne 20',AN,'\nAbondance par stations'),subtitle='(Taux de lecture indiqués pour station si <100%)',
       caption='Data source: Ifremer') +
  theme(legend.position = "bottom", plot.subtitle=element_text(size=8))
dev.off()

##############
##############    Partie C - Estimation des Variances d'échantillonnage
##############

### Variances par bootstrap
# N = nombre de réplicats à définir
N <- 10000

nZones <- length(table(tabl$Zone))
vj_a <- numeric(nZones)
va_a <- numeric(nZones)
abJu <- list()
abAd <- list()

for (i in 1:nZones) {
  abjuv <- numeric(N)
  abAdult <- numeric(N)
  for (j in 1:N) {
    zone_i <- names(table(tabl$Zone))[i]
    y <- nbJuvStation[tabl$Zone %in% zone_i]
    z <- nbAdultStation[tabl$Zone %in% zone_i]
    x <- surfCompte[tabl$Zone %in% zone_i]
    lesNonNA <- which(!(is.na(x)))
    y <- y[lesNonNA]; z = z[lesNonNA]; x <- x[lesNonNA]
    idx <- sample(1:length(y), length(y), replace=TRUE)
    yy <- y[idx]
    zz <- z[idx]
    xx <- x[idx]
    juvRatio <- sum(yy)/sum(xx)
    adultRatio <- sum(zz)/sum(xx)
    abJuv[j] <- juvRatio * Tzone[zone_i] * as.numeric(param$surfCarre)
    abAdult[j] <-  adultRatio * Tzone[zone_i] * as.numeric(param$surfCarre)
  }
  abJu[[i]] <- abJuv
  names(abJu[i]) <- names(Tzone[zone_i])
  abAd[[i]] <- abAdult
  names(abAdult[i]) <- names(Tzone[zone_i])
}

ic_biomasse <- array(dim=c(nZones,8))
for (i in 1:nZones) {
  pipo <- abAd[[i]][order(abAd[[i]])]
  ic_biomasse[i,] <- c(pipo[c(N*0.025,N*0.050,N*0.100)]* as.numeric(param$pdsAdultes)/1000, bioAdult[i], pipo[c(N*0.5,N*0.900,N*0.950,N*0.975)] * as.numeric(param$pdsAdultes)/1000)
}

dimnames(ic_biomasse) <- list(c(names(table(tabl$Zone))),c('IC_95%_min','IC_90%_min','IC_80%_min','Moyenne','Mediane','IC_80%_max','IC_90%_max','IC_95%_max'))
## RESULTAT BIOMASSE ADULTE AVEC IC BOOTSTRAP-------------------------------------
floor(ic_biomasse)
## ------------------------------------------------------------------------


dfAd <- data.frame(zone=rep(names(table(tabl$Zone)),each=N), replicat=unlist(abAd))
dfJu <- data.frame(zone=rep(names(table(tabl$Zone)),each=N), replicat=unlist(abJu))
quiPas0 <- which(tapply(dfAd$replicat, dfAd$zone, mean)>0)
qui0 <- which(tapply(dfAd$replicat, dfAd$zone, mean)==0)
dfAd <- dfAd[dfAd$zone %in% names(quiPas0),]
if (is.factor(dfAd$Zone) & length(qui0) > 0) dfAd$zone <- droplevels(dfAd$zone)

## FIGURE : ABONDANCE ADULTES
dfAd2 <- dfAd %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_figures,'/Abondance ',AN,' adultes.jpeg'), width=1.5*fig.width, height=fig.height, res=fig.res)
ggplot(data=dfAd, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  facet_wrap(~zone, scales='free') +
  geom_vline(data=dfAd2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfAd2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='Abondance (réplicats bootstrap)', y='Fréquence', title=paste0("Holothuries SPM - campagne 20",AN,'\nDistribution des valeurs estimées d\'abondance des adultes'),
       subtitle='Lignes verticales hachurées = 5ème percentile, médiane, moyenne (en bleu) et 95ème percentile', caption='Data source: Ifremer') +
  theme(plot.subtitle=element_text(size=8))
dev.off()

# FIGURE : BIOMASSE ADULTES
dfBio <- dfAd
dfBio$replicat <- dfAd$replicat*.4*10^-6
dfBio2 <- dfBio %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_figures,'/Biomasse ',AN,' adultes.jpeg'), width=1.5*fig.width, height=fig.height, res=fig.res)
ggplot(data=dfBio, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  facet_wrap(~zone, scales='free') +
  geom_vline(data=dfBio2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfBio2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='biomasse (kt)', y='Fréquence', title=paste0("Holothuries SPM - campagne 20",AN,'\nDistribution des valeurs estimées de biomasse des adultes'),
     subtitle='Lignes verticales hachurées = 5ème percentile, médiane, moyenne (en bleu) et 95ème percentile', caption='Data source: Ifremer') +
  theme(plot.subtitle=element_text(size=8))
dev.off()

# FIGURE : ABONDANCE JUVENILES
qui0 <- which(tapply(dfJu$replicat, dfJu$zone, mean)==0)
dfJu <- dfJu[dfJu$zone %in% names(quiPas0),]
if (is.factor(dfJu$zone) & length(qui0) > 0) dfJu$zone <- droplevels(dfJu$zone)

dfJu2 <- dfJu %>% group_by(zone) %>% summarise(q05=quantile(replicat,0.05), mediane=quantile(replicat,0.5),moyenne=mean(replicat), q95=quantile(replicat,.95))
jpeg(paste0(wd_figures,'/Abondance ',AN,' juvéniles.jpeg'),  width=1.5*fig.width, height=fig.height, res=fig.res)
ggplot(data=dfJu, aes(x=replicat)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='red') +
  facet_wrap(~zone, scales='free') +
  geom_vline(data=dfJu2, aes(xintercept=c(mediane)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(q95)), color="black", linetype="dashed") +
  geom_vline(data=dfJu2, aes(xintercept=c(moyenne)), color="blue", linetype="dashed") +
  labs(x='Abondance (réplicats bootstrap)', y='Fréquence', title=paste0("Holothuries SPM - campagne 20",AN,'\nDistribution des valeurs estimées d\'abondance des juvéniles'),
       subtitle='Lignes verticales hachurées = 5ème percentile, médiane, moyenne (en bleu) et 95ème percentile', caption='Data source: Ifremer') +
  theme(plot.subtitle=element_text(size=8))
dev.off()

dfTot <- rbind.data.frame(data.frame(dfAd, stade='Adultes'), data.frame(dfJu, stade='Juvéniles'))
dfTot$stade <- factor(dfTot$stade, ordered=TRUE, levels = c('Juvéniles','Adultes'))

# FIGURE : BOXPLOT ABONDANCE
jpeg(paste0(wd_figures,'/boxplot abondance ',AN,'.jpeg'), width=1000, height=400)
ggplot(data=dfTot, aes(x=stade, y=replicat, fill=zone)) +
  geom_boxplot() +
  labs(x='', y='Nombre total estimé', title=paste0("Holothuries SPM - campagne 20",AN,'\nabondance comparée Juvéniles/Adultes par zone'),
       caption='Data source: Ifremer') +
  facet_wrap(~stade, scales='free') +
  scale_fill_brewer()
dev.off()

save(dfAd,dfAd2, dfJu, dfJu2, dfTot, tabl, file=paste0(wd_path,'/Boot 20',AN,'.Rdata'))

##############
##############    Partie D - Regle de contrôle de gestion (HCR)
##############


Tmax <- 0.015  #Taux d'exploitation donnant une valeur de biomasse moyenne = Vmax

dfAdult <- data.frame(annee=2000+AN,dfAd); dfJuv <- data.frame(annee=2000+AN,dfJu); dfTotal <- data.frame(annee=2000+AN,dfTot); statAdult <- dfAd2; statJuv <- dfJu2
#dfBmax <- dfAd22[dfAd22$zone %in% c('Tuyau Sud','Ouest'),]
#dfBmax$zone <- as.character(dfBmax$zone)
dfAdult$replicatTAC <- Tmax*dfAdult$replicat*param$pdsAdultes*10^-3
tab.B <- dfAdult %>% group_by(zone) %>% summarise(q01 = quantile(replicatTAC, 0.01), q05 = quantile(replicatTAC, 0.05), q10 = quantile(replicatTAC, 0.10),
                                                   q20 = quantile(replicatTAC, 0.2), Moy_TE015=mean(replicatTAC),
												   Bio_000t = mean(replicat*param$pdsAdultes*10^-3), Ab_Adult_millions=mean(replicat*10^-6))
tab.C <- dfJuv %>% group_by(zone) %>% summarise(Ab_Juv_millions=mean(replicat*10^-6))

tab.B <-  cbind.data.frame(tab.B, tab.C)
tab.B[,-9]
write.table(tab.B,file=paste0(wd_figures,'/chiffres_clés 20',AN,'.csv'))

#QUOTA POSSIBLE (Valeur donnant 95% de chance d'être inférieur à valeur de biomaase (Vmax) correspondant à un taux d'exploitation = 1.5% (Tmax) ..................
quota <- 10*round(tab.B$q05/10,0) #Arrondi à 10 tonnes
names(quota) <- dfAd2$zone
quota
#---------------------------------

# FIGURE ILLUSTRANT LE PROPOS
quotaDF <- data.frame(zone=names(quota), replicatTAC=quota)
jpeg(paste0(wd_figures,'/TAC Campagne ',AN,' adultes.jpeg'), width=1.5*fig.width, height=fig.height, res=fig.res)
ggplot(data=dfAdult, aes(x=replicatTAC)) +
  geom_histogram(aes(y=..density..), color='grey',fill='light blue') +
  geom_density(col='blue') +
#  geom_text(data=quotaDF, aes(y=0.00005, label=replicatTAC)) +
  facet_wrap(~zone, scales='free') +
  geom_vline(data=tab.B[,-9], aes(xintercept=c(q05)), color="black", linetype="dashed") +
  geom_vline(data=tab.B[,-9], aes(xintercept=c(Moy_TE015)), color="red", linetype="dashed") +
  labs(x='Quota', y='Fréquence', title=paste0("Holothuries SPM - campagne 20",AN,'\nDistribution des valeurs estimées de quotas centrée sur la moyenne'),
       subtitle='Lignes verticales hachurées = 5ème percentile (bleu) et médiane (rouge)', caption='Data source: Ifremer') +
  theme(plot.subtitle=element_text(size=8))
dev.off()


##############
##############   Partie E - Analyse en inter annuel
##############

load(paste0(wd_path,'/Boot 2021.Rdata'))
dfAd21 <- data.frame(annee=2021,dfAd); dfJu21 <- data.frame(annee=2021,dfJu); dfTot21 <- data.frame(annee=2021,dfTot); statAd21 <- dfAd2; statJu21 <- dfJu2
tabltot1 <- data.frame(annee=2021, tabl[,c(1,12,13,16,29)])
load(paste0(wd_path,'/Boot 2022.Rdata'))
dfAd22 <- data.frame(annee=2022,dfAd); dfJu22 <- data.frame(annee=2022,dfJu); dfTot22 <- data.frame(annee=2022,dfTot); statAd22 <- dfAd2; statJu22 <- dfJu2
tabltot2 <- data.frame(annee=2022, tabl[,c(1,8,9,13,25)])
load(paste0(wd_path,'/Boot 2023.Rdata'))
dfAd23 <- data.frame(annee=2023,dfAd); dfJu23 <- data.frame(annee=2023,dfJu); dfTot23 <- data.frame(annee=2023,dfTot); statAd22 <- dfAd2; statJu22 <- dfJu2
tabltot3 <- data.frame(annee=2023, tabl[,c(1,8,9,13,25)])

tabletot <- rbind.data.frame(tabltot1, tabltot2)
tabletot <- rbind.data.frame(tabletot, tabltot3)


dfTotn <- rbind.data.frame(dfTot21, dfTot22)
dfTotn <- rbind.data.frame(dfTotn, dfTot23)
dfTotn$replicat <- dfTotn$replicat/1000000
dfTotn$bio <- dfTotn$replicat * param$pdsAdultes
dfTotn$bio[dfTotn$stade %in% 'Juvéniles'] <- dfTotn$replicat[dfTotn$stade %in% 'Juvéniles'] * param$pdsJuveniles

jpeg(paste0(wd_figures,'/boxplot abondance.jpeg'), width=1000, height=400)
ggplot(data=dfTotn, aes(x=as.character(annee), y=replicat, fill=zone)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='Année', y='Nombre total estimé (en millions)', title=paste0('Holothuries SPM - toutes campagnes\nMédiane et répartition interquartile des abondances'),
       caption='Data source: Ifremer') +
  facet_wrap(~stade+zone, scales='free', ncol=4) +
  scale_fill_brewer()
dev.off()

jpeg(paste0(wd_figures,'/boxplot abondance2.jpeg'), width=1000, height=600)
ggplot(data=dfTotn, aes(x=as.character(annee), y=replicat, fill=stade)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='année', y='Nombre total estimé (en millions)', title=paste0('Holothuries SPM - toutes campagnes\nMédiane et répartition interquartile des abondances'),
       caption='Data source: Ifremer') +
  facet_wrap(~stade+zone, scales='free', ncol=4) +
  scale_fill_brewer()
dev.off()


#Biomasses comparées

jpeg(paste0(wd_figures,'/boxplot biomasse.jpeg'), width=1000, height=400)
ggplot(data=dfTotn, aes(x=as.character(annee), y=bio, fill=zone)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='Année', y='Biomasse Adultes (000t.)', title=paste0('Holothuries SPM - toutes campagnes\nMédiane et répartition interquartile des biomasses'),
       caption='Data source: Ifremer') +
  facet_wrap(~stade+zone, scales='free', ncol=4) +
  scale_fill_brewer()
dev.off()

#filtre zone
dfT <- dfTotn[substring(dfTotn$zone,1,1) %in% 'T',]; vncol=2; vwidth=500


jpeg(paste0(wd_figures,'/boxplot biomasse2.jpeg'), width=vwidth, height=400)
ggplot(data=dfT, aes(x=as.character(annee), y=bio, fill=zone)) +
  geom_boxplot(position = position_dodge(width=1)) +
  labs(x='Année', y='Biomasse Adultes (000t.)', title=paste0('Holothuries SPM - toutes campagnes\nMédiane et répartition interquartile des biomasses'),
       caption='Data source: Ifremer') +
  facet_wrap(~stade+zone, scales='free', ncol=vncol) +
  scale_fill_brewer()
dev.off()

#Différences interannuelles

vsd <- tapply(tabletot$Nbre.juveniles+tabletot$Nbre.adultes, list(tabletot$STN), sd, na.rm=TRUE)
vmean <- tapply(tabletot$Nbre.juveniles+tabletot$Nbre.adultes, list(tabletot$STN), mean, na.rm=TRUE)
vvl <- tapply(tabletot$Nbre.juveniles+tabletot$Nbre.adultes, list(tabletot$STN, tabletot$annee), sum)
myval <- vvl[,3]-vmean
vmean <- tapply(tabletot$Nbre.juveniles, list(tabletot$STN), mean, na.rm=TRUE)
vvl <- tapply(tabletot$Nbre.juveniles, list(tabletot$STN, tabletot$annee), sum)
myval.juv <- vvl[,3]-vmean
vmean <- tapply(tabletot$Nbre.adultes, list(tabletot$STN), mean, na.rm=TRUE)
vvl <- tapply(tabletot$Nbre.adultes, list(tabletot$STN, tabletot$annee), sum)
myval.ad<- vvl[,3]-vmean
myval.df <- data.frame(STN=names(myval), diff.tot=myval, diff.juv=myval.juv, diff.ad=myval.ad)


jpeg(paste0(wd_figures,'/Variations annualles.jpeg'), width=vwidth, height=400)
ggplot(data=tabletot, aes(x=as.character(STN), y=Nbre.adultes, group=annee)) +
  geom_line(aes(color=factor(annee))) +
  geom_line(data = filter(tabletot, annee == "2023"), color='blue',size = 1.05) +
  guides(x = guide_axis(angle = 90))+
  theme(text = element_text(size = 9)) +
  labs(x='Station', y='Nombre Adultes', title=paste0('Holothuries SPM - toutes campagnes\nVariations interannuelles nombre adultes par station'),
       caption='Data source: Ifremer') +
  facet_wrap(~Zone, scales='free', ncol=1) +
  scale_fill_brewer()
dev.off()

tabletot$vcol='Negative'
npos <- names(which(myval>0))
tabletot$vcol[tabletot$STN %in% npos] <- 'Positive'

jpeg(paste0(wd_figures,'/Variations annualles 2.jpeg'), width=vwidth, height=400)
ggplot(data=tabletot, aes(x=as.character(STN), y=Nbre.adultes, fill=vcol)) +
  geom_boxplot(position = position_dodge(width=1)) +
  guides(x = guide_axis(angle = 90))+
  theme(text = element_text(size = 9)) +
  labs(x='Station', y='Nombre Adultes', title=paste0('Holothuries SPM - toutes campagnes\nVariations interannuelles nombre adultes par station (2023 vs moyenne 2021-2022)'),
       caption='Data source: Ifremer') +
  facet_wrap(~Zone, scales='free', ncol=1) +
  scale_fill_brewer()
dev.off()


tabletot2 <- with(tabletot,aggregate(x=list(nb.juv=Nbre.juveniles, nb.adults=Nbre.adultes), by=list(STN=STN, Zone=Zone), sum, na.rm=TRUE))
tabletot2 <- merge(tabletot2, myval.df)
tabletot2$diff.tot.sign <- case_when(tabletot2$diff.tot <0 ~'2023 inférieur', tabletot2$diff.tot >= 0 ~ '2023 supérieur ou égal', is.na(tabletot2$diff.tot) ~'2023 inférieur')
tabletot2$diff.juv.sign <- case_when(tabletot2$diff.juv <0 ~'2023 inférieur', tabletot2$diff.juv >= 0 ~ '2023 supérieur ou égal')
tabletot2$diff.ad.sign <- case_when(tabletot2$diff.ad <0 ~'2023 inférieur', tabletot2$diff.ad >= 0 ~ '2023 supérieur ou égal')


ZonePlot <- 'Tuyau Sud'
ggplot(data=tabletot2[!is.na(tabletot2$diff.tot) & tabletot2$Zone %in% ZonePlot,], aes(x = diff.tot, y = STN, fill = diff.tot.sign)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 7)) +
  facet_share(~diff.tot.sign, dir = "h", scales = "free", reverse_num = TRUE) +  # note: scales = "free"
  #coord_flip() +
  #theme_minimal() +
  labs(y = "Difference (nombre)", x = "Station", title = paste0('Différentiel 2023 Juvéniles + Adultes - ',ZonePlot)) +
  scale_fill_manual(values = c("pink", "blue"))

ggplot(data=tabletot2[!is.na(tabletot2$diff.juv) & tabletot2$Zone %in% ZonePlot,], aes(x = diff.juv, y = STN, fill = diff.juv.sign)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 7)) +
  facet_share(~diff.juv.sign, dir = "h", scales = "free", reverse_num = TRUE) +  # note: scales = "free"
  #coord_flip() +
  #theme_minimal() +
  labs(y = "Difference (nombre)", x = "Station", title = paste0('Différentiel 2023 Juvéniles - ',ZonePlot)) +
  scale_fill_manual(values = c("pink", "blue"))

ggplot(data=tabletot2[!is.na(tabletot2$diff.ad),], aes(x = diff.ad, y = STN, fill = diff.ad.sign)) +
  geom_bar(stat = "identity") +
  theme(text = element_text(size = 7)) +
  facet_share(~diff.ad.sign, dir = "h", scales = "free", reverse_num = TRUE) +  # note: scales = "free"
  facet_wrap(~Zone, ncol = 2) +
  #coord_flip() +
  #theme_minimal() +
  labs(y = "Difference (nombre)", x = "Station", title = paste0('Différentiel 2023 Juvéniles - ',ZonePlot)) +
  scale_fill_manual(values = c("pink", "blue"))

#################
tabletot3 <- data.frame(tabletot2[,c(1,2,5,8)],pop='Juv + Adults')
names(tabletot3) <- c('STN','Zone','diff.lastY','diff.sign','pop')
tabletemp <- data.frame(tabletot2[,c(1,2,6,9)],pop='Juveniles')
names(tabletemp) <- c('STN','Zone','diff.lastY','diff.sign','pop')
tabletot3 <- rbind.data.frame(tabletot3, tabletemp)
tabletemp <- data.frame(tabletot2[,c(1,2,7,10)],pop='Adults')
names(tabletemp) <- c('STN','Zone','diff.lastY','diff.sign','pop')
tabletot3 <- rbind.data.frame(tabletot3, tabletemp)
ind <- which(is.na(tabletot3$diff.lastY))
if (length(ind)>0) tabletot3 <- tabletot3[-ind,]

tabletot4 <- tabletot3[tabletot3$pop %in% c('Juveniles','Adults') & tabletot3$Zone %in% c('Tuyau Nord','Tuyau Sud'),]
tabletot4$Zone <- factor(tabletot4$Zone)
tabletot4$pop <- as.character(tabletot4$pop)

jpeg(paste0(wd_figures,'/Variations annuelles 2023 vs moyenne 2021-2022.jpeg'), width=2000, height = 1800, res=300)
ggplot(na.omit(tabletot4), aes(x=diff.lastY)) +
  geom_bar(data=na.omit(tabletot4[tabletot3$diff.sign=="2023 inférieur",]), aes(y=STN, fill=pop), stat="identity") +
  geom_bar(data=na.omit(tabletot4[tabletot3$diff.sign=="2023 supérieur ou égal",]), aes(y=STN, fill=pop), stat="identity") +
  geom_hline(yintercept=0, colour="white", lwd=1) +
  scale_fill_brewer(palette="Paired") +
  facet_wrap(~Zone, scales='free') +
  theme_minimal() +
  #coord_flip(ylim=c(-101,101))
  #scale_y_continuous(breaks=seq(-100,100,50), labels=c(100,50,0,50,100)) +
  labs(y="Station", x="Nombre") +
  ggtitle('Différentiel 2023 vs moyenne 2021-2022')
dev.off()

#names(tabletot2)[5:10] <- c('difn_tot','difn_juv','difn_ad','difs_tot','difs_juv','difs_ad')
#tabletot2 %>% pivot_longer(names_to = c('difn','difs'), values_to = c('n'), cols=difn_tot:difs_ad, names_pattern='/dif(.)_(.*)')




######### SPATIAL MAPPING
tabletot5 <- merge(tabletot, strZone, by.x='STN',by.y='Strate', all.x=TRUE)
ind <- which(is.na(tabletot5$Latitude))
if (length(ind)>0) tabletot5 <- tabletot5[-ind,]
tabletot5$Nbre.juveniles[is.na(tabletot5$Nbre.juveniles)] <- 0
vcut <- cut(tabletot5$Nbre.adultes,breaks=c(0,1,20,100,1000,10000), include.lowest=TRUE, labels=c('0','[1-20[','[20-100[','[100-1000[','>=1000'))
tabletot5$Nbre.adultes.classes <- vcut
vcut <- cut(tabletot5$Nbre.juveniles,breaks=c(0,1,20,100,1000,10000), include.lowest=TRUE, labels=c('0','[1-20[','[20-100[','[100-1000[','>=1000'))
tabletot5$Nbre.juveniles.classes <- vcut


jpeg(paste0(wd_figures,'/Variations spatiales annuelles adultes.jpeg'), width=2000, height = 2000, res=300)
ggplot(tabletot5, aes(x=-Longitude, y=Latitude)) +
  geom_point(aes(shape=Nbre.adultes.classes, color=Nbre.adultes.classes, size=Nbre.adultes.classes)) +
  scale_shape_manual(values=c(3,rep(16,4)), name='classe de nombre d\'adultes') +
  scale_color_manual(values=c('black','blue1','blue2','blue3','blue4'), name='classe de nombre d\'adultes') +
  scale_size_manual(values=c(1,1,1.5,2,4),name='classe de nombre d\'adultes') +
  facet_wrap(~annee) +
  theme_light(base_size=7) +
  theme(legend.position = 'bottom', legend.background = element_rect(colour = 1)) +
    labs(x='Longitude Ouest', y='Latitude Nord', title='Holothuries SPM - toutes campagnes\nVariations spatiales annuelles du nombre d\'adultes par station',
         caption='Data source: Ifremer')
dev.off()

jpeg(paste0(wd_figures,'/Variations spatiales annuelles juvéniles.jpeg'), width=2000, height = 2000, res=300)
ggplot(tabletot5, aes(x=-Longitude, y=Latitude)) +
  geom_point(aes(shape=Nbre.juveniles.classes, color=Nbre.juveniles.classes, size=Nbre.juveniles.classes)) +
  scale_shape_manual(values=c(3,rep(16,4)), name='classe de nombre de juvéniles') +
  scale_color_manual(values=c('black','blue1','blue2','blue3','blue4'), name='classe de nombre de juvéniles') +
  scale_size_manual(values=c(1,1,1.5,2,4), name='classe de nombre de juvéniles') +
  facet_wrap(~annee) +
  theme_light(base_size=7) +
  theme(legend.position = 'bottom', legend.background = element_rect(colour = 1)) +
  labs(x='Longitude Ouest', y='Latitude Nord', title='Holothuries SPM - toutes campagnes\nVariations spatiales annuelles du nombre de juvéniles par station',
       caption='Data source: Ifremer')
dev.off()

jpeg(paste0(wd_figures,'/carte stations.jpeg'), width=2000, height = 2000, res=300)
ggplot(tabletot5, aes(x=-Longitude, y=Latitude)) +
  geom_text(aes(label=STN), size=3) +
  theme(text = element_text(size = 8))
dev.off()



ggplot(tabletot5, aes(x=-Longitude,y=Latitude, weight= Nbre.adultes)) +
  geom_bin2d(bins = 10) +
  facet_wrap(~annee) +
  scale_fill_gradient(low="lightblue1",high="darkblue") +
  labs(title = "Distribution des concombres adultes",
          x = "Longitude Ouest",
          y = "Latitude Nord",
          fill = "Nombre par station")


ggplot(tabletot5, aes(x=-Longitude,y=Latitude, weight= Nbre.juveniles)) +
  geom_bin2d(bins = 10) +
  facet_wrap(~annee) +
  scale_fill_gradient(low="lightblue1",high="darkblue") +
  labs(title = "Distribution des concombres juvéniles",
       x = "Longitude Ouest",
       y = "Latitude Nord",
       fill = "Nombre par station")


library(ggOceanMaps)
library(ggspatial)
jpeg(paste0(wd_figures,'/carto zone stations 2023 zoom out.jpeg'), width=2000, height = 2000, res=300)
dt <- data.frame(lon = c(-55, -55, -58, -58), lat = c(44, 47.5, 47.5, 44))
#dt <- data.frame(lon = c(-56, -56, -56.5, -56.5), lat = c(45, 46, 46, 45))
basemap(data = dt, bathymetry = TRUE, rotate=TRUE) +
  geom_spatial_point(data = tabletot5 %>% filter(annee == 2023),
                     aes(x = -Longitude, y = Latitude))
dev.off()



library(leaflet)
dt <- data.frame(lon = c(-50, -50, -60, -60), lat = c(40, 50, 50, 40))
leaflet(dt) %>%
  addTiles(urlTemplate = paste0("https://server.arcgisonline.com/ArcGIS/",
                                "rest/services/Ocean_Basemap/MapServer/tile/{z}/{y}/{x}"),
           attribution = paste0("Tiles &copy; Esri &mdash; Sources: ",
                                "GEBCO, NOAA, CHS, OSU, UNH, CSUMB, National Geographic, ",
                                "DeLorme, NAVTEQ, and Esri")) %>%
  addCircleMarkers(lat = ~ lat, lng = ~ lon,
                   weight = 1, radius = 5,
                   color = "red", fillOpacity = 0.5
  )



ggplot(tabletot5, aes(x=-Longitude, y=Latitude) ) +
  geom_bin2d(bins = 70) +
  scale_fill_continuous(type = "viridis") +
  theme_bw()



#### SPATIAL ANALYSIS
library(kriging)
library(sp)
library(gstat)

tabletot5$Long_deci <- floor(tabletot5$Longitude) + (tabletot5$Longitude %% floor(tabletot5$Longitude))/0.6
tabletot5$Lat_deci <- floor(tabletot5$Latitude) + (tabletot5$Latitude %% floor(tabletot5$Latitude))/0.6
tabletot6 <- tabletot5[,c(1:4,6,13,14)]
names(tabletot6)[5]<- 'Zone'
tabletot6$Long_deci <- -tabletot6$Long_deci

#Distribution des variables
hist(log10(sqrt(tabletot6$Nbre.adultes)+1))
hist(log10(sqrt(tabletot6$Nbre.juveniles)+1))

hist(tabletot6$Nbre.adultes^0.2)
qqnorm(tabletot6$Nbre.adultes^0.2)


#transformation en classe SpatialPointsdataframe
coordinates(tabletot6) <- ~Long_deci+Lat_deci
class(tabletot6)
str(tabletot6)
plot(tabletot6, pch=20)
plot(tabletot6[tabletot6$annee == 2021,], pch=20)

tuyau <- matrix(c(-56.14,-56.14,-56.41,-56.41,45.3,46, 46, 45.30), byrow=FALSE, ncol=2)
lines(tuyau)

vgm1 <- variogram(Nbre.juveniles^.2 ~ 1, tabletot6[tabletot6$annee == 2023,])
plot(vgm1)
model.1 <- fit.variogram(vgm1,vgm(1,"Sph",300,1))
plot(vgm1, model=model.1)
plot(vgm1, plot.numbers = TRUE, pch = "+")
vgm2 <- variogram(Nbre.adultes ~ 1, tabletot6, alpha=c(0,45,90,135))
plot(vgm2)
# the following demonstrates plotting of directional models:
model.2 <- vgm(.59,"Sph",926,.06,anis=c(0,0.3))
plot(vgm2, model=model.2)

g = gstat(NULL, "Adultes < 200", I(Nbre.adultes<200)~1, tabletot6)
g = gstat(g, "Adultes < 400", I(Nbre.adultes<400)~1, tabletot6)
g = gstat(g, "Adultes < 1000", I(Nbre.adultes<1000)~1, tabletot6)
# calculate multivariable, directional variogram:
v = variogram(g, alpha=c(0,45,90,135))
plot(v, group.id = FALSE, auto.key = TRUE) # id and id pairs panels
plot(v, group.id = TRUE, auto.key = TRUE)  # direction panels

# variogram maps:
plot(variogram(g, cutoff=1000, width=100, map=TRUE),
     main = "(cross) semivariance maps")
plot(variogram(g, cutoff=1000, width=100, map=TRUE), np=TRUE,
     main = "number of point pairs")













data(meuse.grid)

Vario_RDT <- variogram(Nbre.adultes~1, data=tabletot6)
plot(Vario_RDT)
Vario_fit <- fit.variogram(Vario_RDT, model = vgm(psill=2, model='Sph',range=30, nugget=0.5))


Vario_RDT.fit = fit.variogram(Vario_RDT,
                              ## A semi-variogram model needs to be manually proposed but it is only to drive the fit
                              ## The parameters will certainly be slightly changed after the fit
                              model = vgm(psill=2, ## Partial Sill (do not confound with the sill)
                                          model="Sph",  ## A Spherical model seems appropriate
                                          range=1,     ## Portée pratique = Maximal distance of autocorrelation
                                          nugget=0.5))  ## Effet pépite = Small-scale variations

RDT_vario.dir=variogram(Nbre.adultes~1,
                        data=tabletot6,
                        alpha = c(0, 45, 90, 135))
## Semi-variograms will be plotted in the directions 0, 45, 90 and 135°
## The variogram is a symmetric measure so for instance, the directions 45
## and 225° (45+180) will lead to the same results
plot(RDT_vario.dir)


RDT_variodir.fit = vgm(psill=2,  ## A semi-variogram model is estimated by eye as before
                       model="Sph",
                       range=30000,
                       nugget=0.5,
                       ## Set the anisotropy parameters. They also have to be estimated by eye
                       anis = c(90,  ## Direction of larger range regarding all the directional variograms
                                ## Here, the direction 45° seems to be that of larger range
                                0.6)) ## Ratio of the minor range (approximately 25m)
## to the maximal range (approximately 40m)

## Plot the directional variogram and associated model
plot(RDT_vario.dir, RDT_variodir.fit)

