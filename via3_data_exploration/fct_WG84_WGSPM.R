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
                              src.proj = 4326,
                              dst.proj = 4467) {
  require(sf)
  #data is the shape (id,lon,lat) in this order
  names(data) <- c("id","x","y")
  sf_transform_xy(data, dst.proj, src.proj)
}