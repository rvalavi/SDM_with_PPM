# Note: to generalise the maptools::as.im.RasterLayer to terra::SpatRaster
as.im.SpatRaster <- function (from, factor.col.name = NULL){
  im <- maptools::as.im.RasterLayer(raster::raster(from))
  return(im)
}
