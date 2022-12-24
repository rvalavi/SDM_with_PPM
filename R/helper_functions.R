library(sf)
library(spatstat)
library(terra)

.raster_to_im <- function(x){
  r <- as.data.frame(x, xy = TRUE) 
  im <- spatstat.geom::as.im(r)
  return(im)
}

# Note: to generalise the spatstat.geom::as.im to SpatRaster and RasterLayer
as.im.SpatRaster <- function(x){ .raster_to_im(x) }
as.im.RasterLayer <- function(x){ .raster_to_im(x) }



# CRS transform function for xy dataframe
crs_transform <- function(xy, # a dataframe with two column: x/longitude, and y/latitude
                          in_crs = 4326,
                          out_crs = 4326){
  
  xy <- as.data.frame(xy)
  
  ppt <- sf::st_as_sf(xy, coords = 1:2, crs = in_crs) %>% 
    sf::st_transform(crs = out_crs) %>% 
    sf::st_coordinates() %>% 
    as.data.frame()

  return(ppt)
}



# cleaning occurrence points
clean_occ <- function(
    points, # dataframe of coordinates of occurrences (e.g., logitude and latitude)
    r, # one raster covariate to be used as a mask
    # jitter = 0,
    remove_outside = TRUE, # remove points outside of window
    remove_dup_dist = 0, # remove duplicate points with some distance
    crs = 4326
){
  # for considering jitter look at sf::st_jitter()
  r <- r[[1]]
  
  # # revise this; as in PPM the cell duplicates doesn't matter
  # if(remove_duplicates){
  #   cells <- terra::cellFromXY(r, as.matrix(points))
  #   points <- points[!duplicated(cells),]
  # }
  
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

