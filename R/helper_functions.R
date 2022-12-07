# Note: to generalise the maptools::as.im.RasterLayer to terra::SpatRaster
as.im.SpatRaster <- function (from, factor.col.name = NULL){
  im <- maptools::as.im.RasterLayer(raster::raster(from))
  return(im)
}

points <- cleaned_data
r <- raster_covars[[1]]

clean_occ <- function(
    points, # dataframe of coordinates of occurrences (e.g., logitude and latitude)
    r, # one raster covariate to be used as a mask
    # jitter = 0,
    remove_outside = TRUE, # remove points outside of window
    remove_duplicates = TRUE #remove duplicate points
){
  # for considering jitter look at sf::st_jitter()
  r <- r[[1]]
  
  # revise this; as in PPM the cell duplicates doesn't matter
  if(remove_duplicates){
    cells <- terra::cellFromXY(r, as.matrix(points))
    points <- points[!duplicated(cells),]
  }
  
  if(remove_outside){
    imr <- as.im(r)
    win <- spatstat.geom::as.owin(imr)
    is_inside <- inside.owin(x = points[,1], 
                             y = points[,2],
                             w = win)
    points <- points[is_inside, ]
  }
  
  return(points)
}

