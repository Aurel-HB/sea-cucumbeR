#' degree_decimal
#'
#' @description A fct function that take in input a vector of degres coordinates
#' of shapes XX°XX,XXX and transform them in decimal XX,XXX WG84
#'
#' @param coord a character
#'
#' @return a vector of decimal coordinates
#'
#'

#transform the coordinates to use in sf
deg_dec <- function(coord){
  # coords is the shape xx°xx°xxx 
  dec <- as.numeric(substr(coord, start = 1, stop = 2))
  min <- as.numeric(substr(coord, start = 4, stop = 5))
  
  # adaptation to XX°XX,XXXX or XX°XX,XXX
  min_virgule <- substr(coord, start = 7, stop = nchar(coord))
  min_virgule <- as.numeric(min_virgule)/(10^nchar(min_virgule)) 

  min <- (min+min_virgule)/60
  #sec <- as.numeric(substr(coord, start = 7, stop = nchar(coord)))/3600
  return(dec+min)
}
