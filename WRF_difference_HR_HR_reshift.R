#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
# File: WRF_functions.R
#
#Author: Guido Cervone (cervone@polygonsu.edu) and Yanni Cao (yvc5268@polygonsu.edu)
#        Geoinformatics and Earth Observation Laboratory (http://geolab.polygonsu.edu)
#        Department of Geography and Institute for CyberScience
#        The Pennsylvania State University
#
#



# This file is used to write the output difference once WRF has completed the runs

require(rgdal)
library(ncdf)
require(tools)
require(raster)

#using the new source code "WRF_function_result_Yanni"since WRF output files uses XLAT, XLAT
#describe CLAT,CLONG in WRF input file
#
source("~/geolab/projects/Marcellus/Reprojection/WRF_functions.R")
source("~/geolab/projects/Rbin/RGIS_functions.R")

CRS.WGS84  = CRS("+proj=longlat +datum=WGS84")

processWRF = function(dir1, dir2, dir3, id, WRF.layer, layer ) {

  # Check if we have already created this file... and if so, just skip it
  #
  ofile = paste("raster.HR.crop.",WRF.layer,"_",id,"_",layer,".tif",sep="")
  if (file.exists(ofile)) {
    return
  }
  
  print(paste(dir1, dir2, dir3, id, WRF.layer, layer))
  
  files                  = list.files(path=dir1, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR              = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  ##files                  = list.files(path=dir2, pattern="\\.nc$",full.names=T)
  ##filename.wrf           = files[id]
  ##raster.HR.SHIFT        = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  files                  = list.files(path=dir3, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.DEFAULT         = getRasterWRF( filename.wrf, WRF.layer, layer)
  
  # Shift up (reshift)
  #
  raster.HR.RESHIFT      = raster.HR.SHIFT
  extent(raster.HR.RESHIFT)[3:4] = transform2spheroid( extent(raster.HR.RESHIFT)[3:4])
  
  
  extent1 = extent(raster.HR)
  extent2 = extent(raster.HR.RESHIFT)
  
  extent3 = extent(extent1[1], 
                   extent1[2], 
                   max(extent1[3],extent2[3]), 
                   min(extent1[4],extent2[4]) )
  spoly   = extent2SpatialPolygons( extent3, CRS.WGS84 )
  
  #raster.HR.RESHIFT.crop = crop(raster.HR.RESHIFT, spoly, snap="out" )
  #raster.HR.SHIFT.crop   = crop(raster.HR.SHIFT, spoly, snap="out" )
  raster.HR.crop         = crop(raster.HR, spoly, snap="out" )
  raster.DEFAULT.crop    = crop(raster.DEFAULT, spoly, snap="out" )
  
  diff1 = raster.HR.crop - raster.DEFAULT.crop
  #diff2 = raster.HR.crop - raster.HR.SHIFT.crop
  #diff3 = raster.HR.crop - raster.HR.RESHIFT.crop
  
  
  #writeRaster(raster.HR.crop, file=paste("raster.HR.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  #writeRaster(raster.HR.SHIFT.crop, file=paste("raster.HR.SHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  #writeRaster(raster.HR.RESHIFT.crop, file=paste("raster.HR.RESHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(raster.DEFAULT.crop, file=paste(dir3,"raster.DEFAULT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  
  writeRaster(diff1, file=paste(outdir,"diff.HR-DEFAULT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  #writeRaster(diff2, file=paste("diff.HR-HR.SHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  #writeRaster(diff3, file=paste(,"diff.HR-HR.RESHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
}



files                  = list.files(path=dir2, pattern="\\.nc$",full.names=T)
filename.wrf           = files[1]
raster.HR.SHIFT        = getRasterWRF( filename.wrf, WRF.layer, layer)

processWRFFast = function(dir1, dir2, dir3, id, WRF.layer, layer ) {
  
  print(paste(dir1, dir2, dir3, id, WRF.layer, layer))
  
  files                  = list.files(path=dir1, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR              = getRasterWRF( filename.wrf, WRF.layer, layer,"output")
  
  files                  = list.files(path=dir2, pattern="\\.nc$",full.names=T)
  filename.wrf           = files[id]
  raster.HR.SHIFT        = getRasterWRF( filename.wrf, WRF.layer, layer,"output")
  
  # Shift up (reshift)
  #
  raster.HR.RESHIFT      = raster.HR.SHIFT
  extent(raster.HR.RESHIFT)[3:4] = transform2spheroid( extent(raster.HR.RESHIFT)[3:4])
  
  
  extent1 = extent(raster.HR)
  extent2 = extent(raster.HR.RESHIFT)
  
  extent3 = extent(extent1[1], 
                   extent1[2], 
                   max(extent1[3],extent2[3]), 
                   min(extent1[4],extent2[4]) )
  spoly   = extent2SpatialPolygons( extent3, CRS.WGS84 )
  
  raster.HR.crop         = crop(raster.HR, spoly, snap="out" )
  raster.HR.RESHIFT.crop = crop(raster.HR.RESHIFT, spoly, snap="out" )
  extent(raster.HR.RESHIFT.crop) = extent(raster.HR.crop)
  
  diff3 = raster.HR.crop - raster.HR.RESHIFT.crop

  writeRaster(raster.HR.RESHIFT.crop, file=paste("d03.raster.HR.RESHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(raster.HR.crop, file=paste("d03.raster.HR.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  #writeRaster(raster.HR.SHIFT.crop, file=paste("raster.HR.SHIFT.crop.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
  writeRaster(diff3, file=paste(outdir, "d03.diff.HR-HR.RESHIFT.",WRF.layer,"_",id,"_",layer,".tif",sep=""),format="GTiff",overwrite=T)
}


#dir1="~/geolab/data/Marcellus/Reprojection/WRF_result/HR"
dir1="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-May2016/d03_HR/"
#dir2="~/geolab/data/Marcellus/Reprojection/WRF_result/HR_shift"
dir2="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-May2016/d03_HR_SHIFT/"
dir3="~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Dec2015/Default/"
outdir = "~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Comparisons-May2016/"

#layers    = c(1,7,13,19,25)    # vertical layer for the parameter (2 = second layer from surface)
layers      = c(1)
WRF.layers = c("tracer_1","tracer_5")
ids        = seq(1,25)          # The hour for the result (e.g. 30 = 30th hour from beginning of simulation)


for ( layer in layers ) {
  for (WRF.layer in WRF.layers) {
    for (id in ids) {
      processWRFFast(dir1, dir2,  id, WRF.layer, layer)      
    }
  }
}




