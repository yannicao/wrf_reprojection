library(ncdf4)
library(tools)
library(raster)

path       = "~/geolab_storage/data/Marcellus/Reprojection/WRF-Input-Feb2016/"
parameters = c("HGT_M","LU_INDEX","F","E","LAI12M", "ALBEDO12M")
ncfiles    = c("geo_em.d01.nc","geo_em.d02.nc","geo_em.d03.nc")
plot       = TRUE
shift      = TRUE  # Do we generate for the SHIFT case


for ( ncfile in ncfiles ) {
  stub = file_path_sans_ext(ncfile)
  
  # We are assuming that the original nc files are in the parent directory than the
  # one specified in path
  #
  filename.wrf       = paste(dirname(path),"/",ncfile,sep="")    #WRF input NC files
  
  if ( shift == TRUE) {
    filename.wrf.new = paste(path,stub,"-SHIFT.nc",sep="")
    append           = "-SHIFT"
  } else {
    filename.wrf.new = paste(path,stub,"-HR.nc",sep="")
    append           = ""
  }
  
  # First make a copy of the file
  #
  
  # From filename.wrf copy to filename.wrf.new 
  #
  file.copy( filename.wrf, filename.wrf.new,overwrite = TRUE )
  
  WRF.new            = nc_open(filename.wrf.new,write=T)
  WRF                  = nc_open(filename.wrf,    write=F)
  
  for ( WRF.layer in parameters ) {
    print(paste("Working on:", filename.wrf.new, WRF.layer))
    
    # Get the original WRF data layer
    #
    original = ncvar_get(WRF, WRF.layer)
    
    
    if ( WRF.layer == "LAI12M" || WRF.layer == "ALBEDO12M") {
      # 5 is for may... 
      #Because LAI and Albedo are retreived from MODIS. They are varied monthly.
      #Month needs define. 
      filename.data    = paste(path,stub,"_",WRF.layer,"_HR_5",append,".Rdata",sep="")
      load(filename.data)

      temp = WRF.data.HR  # Just replace May
      WRF.data.HR      = original     # Make sure this includes all the months, not only May
      WRF.data.HR[,,5] = temp
      
    } else {
      filename.data    = paste(path,stub,"_",WRF.layer,"_HR",append,".Rdata",sep="")
      load(filename.data)
    }
   
    
    if (plot) {
      layout(t(1:2))
      
      o = original
      n = WRF.data.HR
      
      if ( length(dim(original)) > 2) { 
        o = o[,,5] ; n = n[,,5] 
        }

      plot(flip(raster(t(o)), direction = "y"))   ; title("Original Data")  
      plot(flip(raster(t(n)), direction = "y"))   ; title("New Shifted Data")
    }
    
    # WRF.new is an object of class ncdf
    # WRF.layer is what variable to write the data to
    # WRF.data.HR is the values to be written
    put.var.ncdf(WRF.new, WRF.layer, WRF.data.HR)
    
    
    # Now take care of the landmask
    #
    if ( WRF.layer == "LU_INDEX") {
      WRF.data.HR[WRF.data.HR<=17] = 0
      WRF.data.HR[WRF.data.HR>17]  = 1
      put.var.ncdf(WRF.new, "LANDMASK", WRF.data.HR)
    }

  } 
  close.ncdf(WRF)
  close.ncdf(WRF.new)
}
