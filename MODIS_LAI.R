
require(raster)
library(ncdf)

source("ModisDownload.R")
source("WRF_functions.R")

CRS.WGS84  = CRS("+proj=longlat +datum=WGS84")


Sys.setenv(MRT_HOME="/Users/guido/local/MRT/",
           MRT_DATA_DIR="/Users/guido/local/MRT//data")


ls = ModisDownloadGuido(x="MCD15A2",dates=c('2015.05.1','2015.05.30'), h=10:14, v=3:5)

# Remove those files that do not exist
#
valid = as.vector(sapply(ls[,2], file.exists))
ls    = ls[valid,]

ModisMosaicGuido(ls,mosaic=T, 
              MRTpath="~/local/MRT/bin",pixel_size=.01, proj=T, proj_type="GEO")



# Read the MODIS grid
#
WRF    = open.ncdf("~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc",write=T)
WRF.grid.polygons  = WRF.grid2polygons(WRF)
WRF.grid.sp        = SpatialPolygons(WRF.grid.polygons, proj4string = CRS.WGS84)



# Read the data
#
filename = "merged.tif"
dat      = raster( filename )
dat.crop = crop(dat, WRF.grid.sp)



# Values from http://cybele.bu.edu/modismisr/products/modis/userguide.pdf
#
offset = 0
gain   = 0.10
valid  = c(0,100)


values           = values(dat.crop)

values[ values<valid[1] | values>valid[2] ] = NA
values = values*gain + offset

values(dat.crop)  = values


writeRaster(dat.crop, file="MODIS_LAI.tif",format="GTiff",overwrite=T)
