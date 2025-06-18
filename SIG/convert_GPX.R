library(xml2)
library(methods)
library(sf)

#import the grid and the area
point <- read.csv(
  file = "C:/Users/ahebertb/Documents/mission_spm/point_campagne.csv",
  sep = ";")


# Fonction pour convertir DMS en degrés décimaux
dms_to_decimal <- function(dms) {
  # Extraire les degrés et les minutes
  parts <- strsplit(dms, "[°']")[[1]]
  degrees <- as.numeric(parts[1])
  minutes <- as.numeric(parts[2])
  
  # Calculer les degrés décimaux
  decimal <- degrees + minutes / 60
  
  # Retourner la valeur
  return(decimal)
}

# Appliquer la fonction aux colonnes longitude et latitude
point$longitude_decimal <- sapply(point$Longitude, dms_to_decimal)
point$latitude_decimal <- sapply(point$Latitude, dms_to_decimal)


point_gps <- data.frame(ID = point$Station, Long = point$longitude_decimal,
                        Lat = point$latitude_decimal)

# Créer un nouveau document XML
gpx_doc <- xml_new_document()

# Ajouter la déclaration XML et l'élément racine gpx
gpx_root <- xml_add_child(gpx_doc, "gpx",
                          version = "1.1",
                          creator = "R script",
                          xmlns = "http://www.topografix.com/GPX/1/1")

# Ajouter les métadonnées (optionnel)
metadata <- xml_add_child(gpx_root, "metadata")
xml_add_child(metadata, "name", "Exported Track")

# Ajouter un segment de trace (trkseg)
trk <- xml_add_child(gpx_root, "trk")
trkseg <- xml_add_child(trk, "trkseg")

# Ajouter chaque point de trace (trkpt) avec un identifiant
for (i in 1:nrow(point_gps)) {
  trkpt <- xml_add_child(trkseg, "trkpt",
                         Lat = point_gps$Lat[i],
                         Long = point_gps$Long[i])
  xml_add_child(trkpt, "ID", point_gps$station[i])  # Ajouter un identifiant unique
}


# Sauvegarder le document XML dans un fichier GPX
write_xml(gpx_doc, "track.gpx")
