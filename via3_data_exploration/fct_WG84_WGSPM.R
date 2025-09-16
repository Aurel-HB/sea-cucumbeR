#' transfo_WG84_WGSPM
#'
#' @description A fct function that take in input a dataframe with latitude 
#' and longitude in decimal shape in the WG84 (epsg:4326) referential and 
#' translate them in Saint Pierre and Miquelon projection (epsg:4467)
#'
#' @param data a dataframe (id,lon,lat)
#'
#' @return a dataframe with the longitude and latitude in epsg:4326 and 
#' epsg:4467
#'
#'

transfo_WG84_WGSPM = function(data,
                              src.proj = CRS("+init=epsg:4326"),
                              dst.proj = CRS("+init=epsg:4467")) {
  require(sp)
  #data is the shape (id,lon,lat) in this order
  names(data) <- c("id","lon","lat")
  as.data.frame(
    spTransform(
      SpatialPointsDataFrame(
        coords = data.frame(Xbng = data$lon,
                            Ybng = data$lat),
        data = data.frame(id = data$id,
                          Xlon = data$lon,
                          Ylat = data$lat),
        proj4string = src.proj), dst.proj))
  
}