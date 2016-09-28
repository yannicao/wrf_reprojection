#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
# File: WRF_preprocess.R
#
#Author: Guido Cervone (cervone@psu.edu) and Yanni Cao (yvc5268@psu.edu)
#        Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#        Department of Geography and Institute for CyberScience
#        The Pennsylvania State University
#
#

require(rgdal)
require(raster)
library(ncdf)
library(parallel)        # to make things in parallel (multicore)
library(tools)           # for filename operations


source("WRF_functions.R")
source("../../RBin/RGIS_functions.R")
source("Modis_functions.R")


cores                  = 6

write.shapefile        = TRUE   # Shapefile with original WRF data 
shift.to.sphere        = TRUE  # Perform the correction of the raster data from WGS84 to WRF coordinates
layer.id               = 1      # This is used only if the layer has more than 2 dimensions (e.g. LAI) 




##################################
############ MODIS ###############
##################################
if ( F ) {
month=05
year = 2015
h=10:14; v=3:5


Sys.setenv(MRT_HOME="/Users/guido/local/MRT/",
           MRT_DATA_DIR="/Users/guido/local/MRT//data")


# Values from http://cybele.bu.edu/modismisr/products/modis/userguide.pdf
#
offset   = 0 ; gain   = 0.10 ; valid  = c(0,100)
product  = "MCD15A2"
filename = paste("~/geolab_storage/data/Marcellus/Reprojection/MODIS_LAI12M_2015_",month,".tif",sep="")

MODIS.data( product, tmp.dir = "~/tmp/", month, year, h=h,v=v)
MODIS.calibrate(pattern=".Lai_1km.tif", filename, offset, gain, valid)


res         = raster(filename)
a           = values(res)
a[is.na(a)] = 0
values(res) = a

writeRaster(res, file=filename,format="GTiff",overwrite=T)



# Albedo values from https://www.umb.edu/spectralmass/terra_aqua_modis/v006/mcd43a1_brdif_albedo_model_parameters_product
#
offset = 0 ; gain   = 0.001 ; valid  = c(0,32766)
offset = 0 ; gain   = 0.1 ; valid  = c(0,32766)
product = "MCD43B3"
filename = paste("~/geolab_storage/data/Marcellus/Reprojection/MODIS_Albedo_2015_",month,".tif",sep="")

filename.WSA = "AlbedoMeanWSA.tif"
filename.BSA = "AlbedoMeanBSA.tif"

MODIS.data( product, tmp.dir = "~/tmp/", month, year, h=h,v=v)
MODIS.calibrate(pattern=".WSA_shortwave.tif", filename.WSA, offset, gain, valid)
MODIS.calibrate(pattern=".BSA_shortwave.tif", filename.BSA,  offset, gain, valid)

wsa = raster(filename.WSA)
bsa = raster(filename.BSA)
tmp = stack(wsa,bsa)
res = calc(tmp, mean, na.rm=T)

a           = values(res)
a[is.na(a)] = 5  # WRF uses 5 as minimum
values(res) = a

writeRaster(res, file=filename,format="GTiff",overwrite=T)
}




##################################
############## ALB ###############
##################################

for ( shift.to.sphere in c(T,F)) {
filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_Albedo_2015_5.tif"
WRF.layer        = "ALBEDO12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
               filename.raster = filename.raster,
               WRF.layer = WRF.layer,
               layer.id = layer.id,
               shift.to.sphere = shift.to.sphere,
               write.shapefile = write.shapefile,
               cores = cores) 

filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_Albedo_2015_5.tif"
WRF.layer        = "ALBEDO12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                layer.id = layer.id,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_Albedo_2015_5.tif"
WRF.layer        = "ALBEDO12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                layer.id = layer.id,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 





##################################
############## LAI ###############
##################################

filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_LAI12M_2015_5.tif"
WRF.layer        = "LAI12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                layer.id = layer.id,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 

filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_LAI12M_2015_5.tif"
WRF.layer        = "LAI12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                layer.id = layer.id,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/MODIS_LAI12M_2015_5.tif"
WRF.layer        = "LAI12M"
layer.id         = 5
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                layer.id = layer.id,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


}




##################################
############## D03 ###############
##################################


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/elevation_merged.tif"
WRF.layer        = "HGT_M"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
filename.raster  = "~/geolab_storage/data/NLCD2011/NLCD_merged_wgs84.tif"
filename.raster2 =  "~/geolab_storage/data/Landcover_CanadaUSMexico/Land_Cover_2010_TIFF/LandCover_2010/data/NA_LandCover_2010_25haMMU.tif_WGS84_WRF_NLCDreclass.tif"
WRF.layer        =   "LU_INDEX"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                filename.raster2 = filename.raster2,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
WRF.layer        = "F"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d03.nc"
WRF.layer        = "E"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


##################################
############## D02 ###############
##################################

filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/elevation_merged_1KM.tif"
WRF.layer        = "HGT_M"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
filename.raster  = "~/geolab_storage/data/NLCD2011/NLCD_merged_wgs84.tif"
filename.raster2 =  "~/geolab_storage/data/Landcover_CanadaUSMexico/Land_Cover_2010_TIFF/LandCover_2010/data/NA_LandCover_2010_25haMMU.tif_WGS84_WRF_NLCDreclass.tif"
WRF.layer        =   "LU_INDEX"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                filename.raster2 = filename.raster2,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
WRF.layer        = "F"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d02.nc"
WRF.layer        = "E"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 

##################################
############## D01 ###############
##################################

filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
filename.raster  = "~/geolab_storage/data/Marcellus/Reprojection/elevation_merged_1KM.tif"
WRF.layer        = "HGT_M"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
filename.raster  = "~/geolab_storage/data/NLCD2011/NLCD_merged_wgs84.tif"
filename.raster2 =  "~/geolab_storage/data/Landcover_CanadaUSMexico/Land_Cover_2010_TIFF/LandCover_2010/data/NA_LandCover_2010_25haMMU.tif_WGS84_WRF_NLCDreclass.tif"
WRF.layer        =   "LU_INDEX"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                filename.raster2 = filename.raster2,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
WRF.layer        = "F"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 


filename.wrf     = "~/geolab_storage/data/Marcellus/Reprojection/geo_em.d01.nc"
WRF.layer        = "E"
WRF.preprocess( filename.wrf = filename.wrf, 
                filename.raster = filename.raster,
                WRF.layer = WRF.layer,
                shift.to.sphere = shift.to.sphere,
                write.shapefile = write.shapefile,
                cores = cores) 










# Thisis the code to preprocess the North America data so it can be used.  It should be run only once to generate the
# land cover file that it is used to crop data for Canada and outside what it is covered by the NLCD
#
if ( FALSE ) {
  
  filename.CA = "../../../../data/Landcover_CanadaUSMexico/Land_Cover_2010_TIFF/LandCover_2010/data/NA_LandCover_2010_25haMMU.tif_WGS84.tif"
  filename.CA.recl = paste(file_path_sans_ext(filename.CA),"_WRF_NLCDreclass.tif",sep="")
  filename.NLCD = "../../../../data/NLCD2011/NLCD_merged_wgs84.tif"
  
  
  
  CA   = raster(filename.CA)
  NLCD = raster(filename.NLCD)
  
  WRF                = open.ncdf("../../../../data/Marcellus/Reprojection/geo_em.d01.nc",write=T)
  lats   = range(get.var.ncdf(WRF,"CLAT"))
  longs  = range(get.var.ncdf(WRF,"CLONG"))
  
  
  lats[1] = lats[1] -.05
  lats[2] = lats[2] + 1
  
  longs[1] = longs[1] - 0.5
  longs[2] = longs[2] + 0.5
  
  spoly.lats  = c(lats[1],lats[2],lats[2],lats[1],lats[1])
  spoly.longs = c(longs[1],longs[1],longs[2],longs[2],longs[1]) 
  
  poly   = Polygons( list( Polygon( cbind(spoly.longs,spoly.lats) ) ), 1)
  spoly  = SpatialPolygons( list(poly), proj4string = CRS(proj4string(CA) ))
  
  CA.crop = crop(CA, spoly, snap="out")
  
  
  conv = matrix(c(1:19,
                  42,42,42,41,41,43,52,52,71,71,52,71,31,90,82,31,23,11,12),ncol=2)
  
  for (row in 1:nrow(conv) ) {
    print(row)
    CA.crop[ CA.crop==conv[row,1] ] = conv[row,2]
  }    
  
  
  CA.crop@legend = NLCD@legend
  # Let's crop only the region we need
  #
  writeRaster(CA.crop, filename.CA.recl, format="GTiff" )

}

