# library(maptools)
library(raster)

grid_dir <- "/Users/rvalavi/Dropbox/MyProjects/SDM_with_PPM/data/grids"
vars <- list.files(grid_dir, pattern = ".tif$", full.names = TRUE)

# read the raster layers as a raster stack
r <- stack(vars)
# set the coordinate system to mercator / geographic coordinates
# the epsg code for geographic coordinates system is 4326
crs(r) <- CRS("+init=epsg:4283")
# transfer a metric coordinate system, here UTM South zone 56 is suitable
# the epsg code for this crs is 32756
nsw_stack <- projectRaster(r, 
                           res = c(250,250), 
                           crs = CRS("+init=epsg:32756"))

setwd("/Users/rvalavi/Dropbox/MyProjects/SDM_with_PPM/data/grids_utm")

writeRaster(nsw_stack, 
            bylayer = TRUE, 
            filename = paste0(names(nsw_stack), ".tif"))