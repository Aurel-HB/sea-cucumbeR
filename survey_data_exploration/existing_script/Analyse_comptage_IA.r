# Analyse comptage IA
#
# Campagnes HoloSPM_TV de 2021 à 2023
# Analyse spatiale des holothuries à partir des comptages
# et positionnnement dans l'espace par IA
#
# J. Vigneau, R. Failletaz, J. Simon, E. Foucher
# 2024
#
#

library(readxl)
library(ggplot2)

# 1 - IMPORT DES FICHIERS
# -----------------------
filePath <- 'C:/_PROGRAMME/NAFO/concombre/Analyse comptage IA/'
figPath <- paste0(filePath,'Figures/')
compte <- read.csv(paste0(filePath,'SPM_IA_predictions.csv'), header=TRUE, sep=',')
names(compte) <- c('id','station','nFrames','species','confidence','positionX','sizeFrame','time','year')
compte$time <- as.POSIXlt(as.character(compte$time), format = '%H:%M:%S')
coord <- read_xlsx(paste0(filePath, 'coordonnées holospm2022.xlsx'))
names(coord) <- c('date','trait','station','timeEndShoot','latEndShoot','longEndShoot','depthEnd','timeBegHaul','latBegHaul','longBegHaul','depthBeg','vesselSpeed','seaStatus','validStation','comment','tempSurf')



# 2 - VISUALISATION
#
# 2.1 - Données de spatialisation empiriques
#
for (yr in 2021:2023) {
  for (st in unique(compte$station[compte$year %in% yr])) {
    dfTemp <- compte[compte$year %in% yr & compte$station %in% st,]
    myfig <- ggplot(dfTemp, aes(x=positionX, y=time$min + time$sec/60, shape=species, colour=species, size=sizeFrame)) +
      geom_point() +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
      labs(x='Position X', y='Position temps', title=paste0('Spatialisation holothurie\n',yr,' - station ',st,' (N=',nrow(dfTemp),')'), caption='Data source: Ifremer') +
      theme(legend.position="bottom") +  #ça ne marche pas!
      theme_light()
    jpeg(paste0(figPath,yr,'-',st,'-individual_position.jpeg'), width=300, height=800)
    print(myfig)
    dev.off()

  }
}


# 2.2 - Essai de distinguer des classes de taille à partir de sizeFrame
#
jpeg(paste0(figPath,'histogramme des frames vus par IA.jpeg'), width=800, height=700)
hist(floor(compte$sizeFrame/10000), breaks=200, xlim=c(0,60), col='cadetblue2', border='cornflowerblue')
dev.off()

compte$sizeCut <- cut(floor(compte$sizeFrame/10000), breaks=c(0,8,20,10000), labels=c('Petits','Moyens','Gros'))
compte$sizeCut2 <- cut(floor(compte$sizeFrame/10000), breaks=c(0,9,10000), labels=c('Juvéniles','Adultes'))
compte$sizeCut2 <- factor(compte$sizeCut2, levels = c('Adultes','Juvéniles'))

#Optimisation du cut avec les comptages manuels
# A faire

for (st in unique(compte$station)) {
    dfTemp <- compte[compte$station %in% st & substring(compte$species,1,4) %in% 'cucu',]
    nCol=length(unique(dfTemp$year))
    nbAdultes <- paste(table(dfTemp$sizeCut2,dfTemp$year)[1,], collapse = ', ')
    nbJuv <- paste(table(dfTemp$sizeCut2,dfTemp$year)[2,], collapse = ', ')
    myfig <- ggplot(dfTemp, aes(x=positionX, y=time$min + time$sec/60, shape=sizeCut2, colour=sizeCut2, size=sizeFrame)) +
      geom_point(alpha=1/3) +
      facet_grid(cols = vars(year)) +
      theme(axis.ticks = element_blank(), axis.text.x = element_blank()) +
      labs(x='Position X', y='Position temps', title=paste0('HOLOSPM-TV - Spatialisation holothurie, Station ',st), subtitle=paste0('\nNb adultes = ',nbAdultes, '\nNb Juvéniles = ',nbJuv), caption='Data source: Ifremer') +
      theme(legend.position="bottom") +  #ça ne marche pas!
      theme_light()
    jpeg(paste0(figPath,'Station ',st,'- sizeCut - individual_position.jpeg'), width=nCol*300, height=900)
    print(myfig)
    dev.off()
}




# 2.2 - Intégration coordonnées latitude et longitude
#



