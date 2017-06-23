#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
# File: WRF_functions.R
#
#Author: Guido Cervone (cervone@psu.edu) and Yanni Cao (yvc5268@psu.edu)
#        Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#        Department of Geography and Institute for CyberScience
#        The Pennsylvania State University
#
#

library(tools)
library(ncdf4)
library(raster)
library(rgdal)

CRS.WGS84  = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
CRS.WRF    = CRS("+proj=lcc +lat_1=30 +lat_2=60 +lat_0=41.8389129639 +lon_0=-77 +x_0=0 +y_0=0 +a=6370000 +b=6370000 +units=m +no_defs")

WRF.qc.obs = function(file) {
  
  data = scan(file,what="character",sep="?")
  
  coords = NULL
  values = list()
  mat    = NULL
  
  for ( i in 1:length(data) ) 
  {
    temp = unlist(strsplit(unlist(strsplit(data[i],"[-]") )," "))
    temp = temp[temp!=""]
    
    first = as.numeric(temp[1])
    
    if ( !is.na(first) && first < 100 && first > 20) {
      
      values[[length(values)+1]] = mat
      second =  as.numeric(temp[2])
      coords = rbind(coords, cbind(as.numeric(temp[2]), first))
      
      mat = NULL
      
    } else if ( !is.na(first) && first > 1000 & first !=777777) {
      mat = rbind(mat,temp)
    }
  }
  
  values[[length(values)+1]] = mat
  
  #========= generate a table which include lat, long, temp, windspeed,winddir========
  #http://www2.mmm.ucar.edu/wrf/users/docs/user_guide_V3/users_guide_chap7.htm#format
  #values[[10]][,5] #gives a table with the observations.
  #coords[10,]  #gives you long and lat
  #In value
  #1,2,  represent  Pressure (Pa) of observation, and QC
  #3,4   represent  Height (m MSL) of observation, and QC
  #5,6   represent  Temperature (K) and QC
  #7,8   represent  Dewpoint (K) and QC
  #9,10  represent  Wind speed (m/s) and QC
  #11,12 represent  Wind direction (degrees) and QC
  
  coords<-data.frame(coords)               #convert to data frame
  coords$ID<-seq.int(nrow(coords))         #generate ID
  colnames(coords) <- c("long","lat","ID") # rename
  coords$long <- -coords$long              #It is in west, long is negative
  
  #select the data within domain 03
  #
  obs<-coords[which(coords$long < -75.05096 & coords$long > -78.08807 & coords$lat > 41.1605 & coords$lat < 42.69054),] #gives you long and lat
  
  #====generate table===
  #
  vtemp <- NULL
  vwsp  <- NULL
  vwdir <- NULL
  vhgt  <- NULL
  for (i in obs$ID) {
    temp  <- values[[i]][,5]
    wsp   <- values[[i]][,9]
    wdir  <- values[[i]][,11]
    hgt   <- values[[i]][,3]
    vtemp <- rbind(vtemp,temp)
    vwsp  <- rbind(vwsp,wsp)
    vwdir <- rbind(vwdir,wdir)
    vhgt  <- rbind(vhgt,hgt)
  }
  
  obs <- cbind(obs,vtemp,vwsp,vwdir,vhgt)
  obs = obs[!obs[,4] == "888888.00000",]              #remove null column
  obs[,4] <- as.numeric(as.character(obs[,4]))-273.15 #convert from Kevin to celsis degree
  colnames(obs) <- c("long","lat","ID","Temp_C", "WindSpeed_m_s", "WindDir_degree","Height_m")
  
  return(obs)
}


WRF.shift = function( data.raster, filename = "" ) {

  # Check if the converted file exists
  #
  if ( filename == "" || !file.exists(filename) ) {
    
    extent(data.raster)[3:4] = transform2sphere( extent(data.raster)[3:4] )
    
    if ( filename != "") {
      writeRaster(data.raster, filename, format="GTiff" )
    }
  } else {
    data.raster = raster(filename)
  }
   
  return(data.raster)
}




deg2rad <- function(deg) return(deg*pi/180)
rad2deg <- function(rad) return(rad*180/pi)

transform2sphere = function( lat ) 
{
  a=6378137 
  b=6356752.314
  
  lat.rad = deg2rad(lat)
  ret.rad = ( atan(tan(lat.rad) * b^2 / a^2))
  
  return( rad2deg(ret.rad) )
}

transform2spheroid = function( lat ) 
{
  a=6378137 
  b=6356752.314
  
  lat.rad = deg2rad(lat)
  ret.rad = ( atan(tan(lat.rad) * a^2 / b^2))
  
  return( rad2deg(ret.rad) )
}


# The input is the ncdf connection to a WRF file
#
WRF.grid2polygons = function(WRF, type="input") {
  
  # Read the WRF domain definition
  #
  if ( type == "input") {
    lats.c   = ncvar_get(WRF,"CLAT")
    longs.c  = ncvar_get(WRF,"CLONG")
  } else {
    lats.c   = ncvar_get(WRF,"XLAT")
    longs.c  = ncvar_get(WRF,"XLONG")
  }
  
  lats.u   = ncvar_get(WRF,"XLAT_U")
  longs.u  = ncvar_get(WRF,"XLONG_U")
  
  lats.v   = ncvar_get(WRF,"XLAT_V")
  longs.v  = ncvar_get(WRF,"XLONG_V")
 
  #Convert the grid to shapefiles
  #
  polygons = list()
  index    = 1
  
 
  # Let's define the polygons (grid cells)
  # This does some basic interpolation to improve the visualization 
  #
  for ( row in 1:nrow(longs.c) ) {
    print(paste("Analyzing Row",row,"of",nrow(longs.c)))  
    
    for ( col in 1:ncol(longs.c) ) {  
      
      longs.LL = longs.u[row,col]   
      lats.LL  = lats.v[row,col] 
      lats.UL  = lats.v[row,col+1]
      longs.LR = longs.u[row+1,col] 
      
      
      if (col < ncol(longs.c) ) {
        longs.UL = longs.u[row,col+1]
        longs.UR = longs.u[row+1,col+1]    
      } else {
        longs.UL = longs.u[row,col]   + ( longs.u[row,col-1] - longs.u[row,col-2])
        longs.UR = longs.u[row+1,col] + ( longs.u[row+1,col-1] - longs.u[row+1,col-2])
      }
      
      if (row < nrow(longs.c) ) {
        lats.UR = lats.v[row+1,col+1]
        lats.LR = lats.v[row+1,col]
      } else {
        lats.UR = lats.v[row,col+1] + ( lats.v[row-1, col+1] - lats.v[row-2, col+1])
        lats.LR = lats.v[row,col]   + ( lats.v[row-1, col] - lats.v[row-2, col])     
      }
      
      longs = c( longs.LL, longs.UL, longs.UR, longs.LR, longs.LL)
      lats  = c( lats.LL,  lats.UL,  lats.UR,  lats.LR,  lats.LL)
      
      #longs = c( longs.u[row,col], longs.u[row,c2], longs.u[row+1,c2], longs.u[row+1,col], longs.u[row,col])
      #lats = c( lats.v[row,col], lats.v[row,col+1], lats.v[r2,col+1], lats.v[r2,col],lats.v[row,col])
      
      polygons[[ index ]]  = Polygons( list( Polygon( cbind(longs,lats) ) ),index)  
      index=index+1
    }
  }
  
  return(polygons)
}




WRF.preprocess = function(  filename.wrf, filename.raster, filename.raster2="", WRF.layer, layer.id = 1, shift.to.sphere = FALSE, write.shapefile = FALSE, cores = 12 ) {
  
  filename.raster.esri   = paste(basename(file_path_sans_ext(filename.wrf)),"_",WRF.layer,"_HR",sep="")
  filename.wrf.esri      = paste(basename(file_path_sans_ext(filename.wrf)),"_",WRF.layer,"_WRF",sep="")
  filename.wrf.grid      = paste(filename.wrf,"_grid.Rdata",sep="")
  
  # Read the ncdf data
  #
  WRF                = nc_open(filename.wrf,write=F)
  
  # Get the elevation and the grid
  #
  WRF.data           = ncvar_get(WRF,WRF.layer)
  
  # This is in case the WRF.data has more than 2 dimensions
  #
  if ( length(dim(WRF.data))>2) {
    WRF.data = WRF.data[,,layer.id]
    filename.wrf.esri      = paste(filename.wrf.esri,layer.id,sep="_")
    filename.raster.esri   = paste(filename.raster.esri,layer.id,sep="_")
  }
  
  
  # This i
  if ( !file.exists(filename.wrf.grid) ) {
    WRF.grid.polygons  = WRF.grid2polygons(WRF)
    WRF.grid.sp        = SpatialPolygons(WRF.grid.polygons, proj4string = CRS.WGS84)
    save(WRF.grid.sp, file=filename.wrf.grid)
  } else {
    load(filename.wrf.grid)
  }
  
  if ( WRF.layer == "F" || WRF.layer =="E") {
    data.raster   = ncvar_get(WRF,"CLAT")
  } else {
    # Read the raster data and make sure it is in WGS84
    #
    data.raster  =  raster(filename.raster)
    data.raster  =  RGIS.raster2WGS84( data.raster )  
    
    if ( filename.raster2 != "") {
      data.raster2  =  raster(filename.raster2)
      data.raster2  =  RGIS.raster2WGS84( data.raster2 )  
      
    }
    
  }
  
  
  
  # Check if we need to reproject the data in WRF spherical coordinates
  #
  if ( shift.to.sphere == TRUE ) {
    print("Converting raster to WRF coordinates")
    
    if ( WRF.layer == "F" || WRF.layer =="E") {
      data.raster = transform2sphere( data.raster )
    } else {
      extent(data.raster)[3:4] = transform2sphere( extent(data.raster)[3:4] )
      
      if ( filename.raster2 != "") {
        extent(data.raster2)[3:4] = transform2sphere( extent(data.raster2)[3:4] )
      }
      
    }
    filename.raster.esri  = paste(filename.raster.esri,"-SHIFT",sep="")
  }
  
  
  
  # Compute the average for each grid cell using the higher resolution DEM using parallel computation
  # Remember that in all computation we must follow the order in which the WRF.grid is created.   Because
  # We assume the matrix by row (we loop first the rows and then the cols when we call grid2polygons)
  # all computations must be done by row (hence the argument passed to matrix)
  #
  print("Computing the masking of the raster data to the WRF grid")
  
  if ( WRF.layer == "LU_INDEX" ) {
    res       = mclapply(1:length(WRF.grid.sp), RGIS.rasterSpatialPolygonsNLCDIndex2, polygons=WRF.grid.sp, data.raster=data.raster, data.raster2 = data.raster2, mc.cores=cores)
    res.mat   = matrix(as.numeric(unlist(res)),ncol=2,byrow=T)  
    
    # Convert to WRF cover classes
    #
    conv = matrix(c(95,90,82,81,71,52,43,42,41,31,24,23,22,21,11,
                    40,39,38,37,33,32,30,29,28,27,26,25,24,23,17),ncol=2)
    
    for (row in 1:nrow(conv) ) {
      res.mat[ res.mat[,1]==conv[row,1] ,1] = conv[row,2]
    }    
    
  } else if ( WRF.layer == "F" || WRF.layer == "E" ) {
    time = 23.9344699*60*60
    
    #Recalculate the Coriolis F, E
    #
    if ( WRF.layer == "F" ) {
      coriolis         = 4*pi / (time) * sin( deg2rad(data.raster) ) # always 1. use radius not degree 2. check the unit
    } else {
      coriolis         = 4*pi / (time) * cos( deg2rad(data.raster) )
    }
    res.mat = cbind( as.vector(t(coriolis)),1:length(coriolis) )
    
  } else {
    res       = mclapply(1:length(WRF.grid.sp), 
                         RGIS.rasterSpatialPolygonsAverageIndex, 
                         polygons=WRF.grid.sp, data.raster=data.raster, na.rm=F, 
                         mc.cores=cores)
    res.mat   = matrix(as.numeric(unlist(res)),ncol=2,byrow=T)  
  }
  
  print("Saving files")  
  
  # To write the shapefiles
  #
  if ( write.shapefile ) {
    
    res.df                 = data.frame(CellID = res.mat[,2], Data = res.mat[,1])
    WRF.grid.spdf          = SpatialPolygonsDataFrame(WRF.grid.sp[1:nrow(res.df)], data=res.df) 
    writeOGR(WRF.grid.spdf, dirname(filename.wrf), filename.raster.esri, driver="ESRI Shapefile",overwrite_layer=T)
    
    res.df                 = data.frame(CellID = 1:length(WRF.data), Data=as.vector(t(WRF.data)))   
    WRF.grid.spdf          = SpatialPolygonsDataFrame(WRF.grid.sp, data=res.df)
    writeOGR(WRF.grid.spdf, dirname(filename.wrf), filename.wrf.esri, driver="ESRI Shapefile",overwrite_layer=T)
  }
  
  
  WRF.data.HR            = t(matrix(res.mat[,1],nrow=ncol(WRF.data)))
  filename.Rdata         = paste(dirname(filename.wrf),"/",filename.raster.esri,".Rdata",sep="")
  save(WRF.data.HR, file = filename.Rdata)
  
  close.ncdf(WRF)
}





WRF.plotLayer = function(  filename.wrf,  WRF.layer, layer.id=1 ) {
    
  filename.wrf.esri      = paste(basename(file_path_sans_ext(filename.wrf)),"_",WRF.layer,"_WRF",sep="")

  WRF.grid.spdf = getDataWRF( filename.wrf, WRF.layer="Temp", layer=2 )[[1]]
  writeOGR(WRF.grid.spdf, dirname(filename.wrf), filename.wrf.esri, driver="ESRI Shapefile",overwrite_layer=T)
}


UVtoSpd = function(U, V) {
  valid           = !is.na(U) & !is.na(V)
  ret             = rep(NA, length(U))
  ret[valid] = sqrt(U[valid]^2+V[valid]^2)
  return(ret)
}


UVtoDir = function ( U, V ) {
  
  
  valid           = !is.na(U) & !is.na(V)
  angle           = rep(NA, length(U))
  
  angle[valid]    = (180/pi * atan2(U[valid],V[valid])) %% 360
  angle[valid][angle[valid]<0]  = angle[valid][angle[valid]<0]+360
  
  return( angle )
}

getDataWRF = function( filename.wrf, WRF.layer="Temp", layer=2, type="input" ) {
  

  print(paste("Processing:",filename.wrf,WRF.layer,layer))
  
  WRF                = nc_open(filename.wrf,write=F)
  
  temp               = ncvar_get(WRF,"CLAT")
  
  if ( WRF.layer == "Temp") {
    #Calculate temprature, combing:
    #http://gradsusr.org/pipermail/gradsusr/2011-December/031698.html 
    ###this one is wrong since WRF Pressure is in pa
    ###it should be converted to mbar. so I used (PB+B)/100
    #http://mailman.ucar.edu/pipermail/wrf-users/2010/001896.html
    P               = ncvar_get(WRF,"P")[,,layer]
    PB              = ncvar_get(WRF,"PB")[,,layer]
    TotalPotentialTemperature = T+300
    WRF.data     = TotalPotentialTemperature*((((PB+P)/100)/1000)^(2/7))-273.15
  } else if ( WRF.layer == "WindSpeed"){
    U               = ncvar_get(WRF,"U")[,,layer]
    U               = U[1:nrow(temp),1:ncol(temp)]
    V               = ncvar_get(WRF,"V")[,,layer]
    V               = V[1:nrow(temp),1:ncol(temp)]
    spd             = UVtoSpd(as.vector(U),as.vector(V))
    WRF.data        = matrix(spd,nrow = nrow(temp),ncol =ncol(temp))
  } else if ( WRF.layer == "WindDir"){
    U               = ncvar_get(WRF,"U")[,,layer]
    U               = U[1:nrow(temp),1:ncol(temp)]
    V               = ncvar_get(WRF,"V")[,,layer]
    V               = V[1:nrow(temp),1:ncol(temp)]
    dir             = UVtoDir(as.vector(U),as.vector(V))
    WRF.data        = matrix(dir,nrow = nrow(temp),ncol =ncol(temp))
  } else {
    WRF.data     = ncvar_get(WRF,WRF.layer)
    if ( length(dim(WRF.data)) > 2) { WRF.data = WRF.data[,,layer]}
    
    WRF.data     = WRF.data[1:nrow(temp),1:ncol(temp)]
  }
  
  
  filename.wrf.grid      = paste(filename.wrf,"_grid.Rdata",sep="")
  
  if ( !file.exists(filename.wrf.grid) ) {
    WRF.grid.polygons  = WRF.grid2polygons(WRF, type=type)
    WRF.grid.sp        = SpatialPolygons(WRF.grid.polygons, proj4string = CRS.WGS84)
    save(WRF.grid.sp, file=filename.wrf.grid)
  } else {
    load(filename.wrf.grid)
  }
  
  res.df                 = data.frame(CellID = 1:length(WRF.data), Data=as.vector(t(WRF.data)))   
  WRF.grid.spdf          = SpatialPolygonsDataFrame(WRF.grid.sp, data=res.df)
  
  return( list(WRF.grid.spdf, WRF.data ) )
}


getRasterWRF = function( filename.wrf, WRF.layer="Temp", layer=2, type="input" ) {
  
  temp           = getDataWRF(filename.wrf, WRF.layer, layer, type)
  WRF.grid.spdf  = temp[[1]]
  WRF.data       = temp[[2]]
  
  raster.temp            = flip(raster(t(WRF.data)),direction="y")
  extent(raster.temp)    = extent(WRF.grid.spdf)
  ncol(raster.temp)      = nrow(WRF.data)  # thins in R are rotated and flipped compared to raster
  nrow(raster.temp)      = ncol(WRF.data)
  
  #return(raster.temp)
  
  data.raster            = rasterize( WRF.grid.spdf, raster.temp, field='Data', fun=mean)
  
  #writeRaster(data.raster, filename = "temp.tif", format="GTiff")
  #CRS.WRF    = CRS("+proj=lcc +lat_1=30 +lat_2=60 +lat_0=41.8389129639 +lon_0=-77 +x_0=0 +y_0=0 +a=6378137 +b=6356752.314 +units=m +no_defs")
  #gdalwarp("temp.tif", "tempout.tif",t_srs=(CRS.WGS84), t_srs=(CRS.WRF), verbose=TRUE, overwrite = FALSE) 
  #tempout    = raster("tempout.tif")
  #extent(tempout) = extent(WRF.grid.spdf)
  
  return (data.raster)
}



diff.circular = function(dir1, dir2 ) {
  #  if ( is.na(dir1) || is.na(dir2) ) { return( NA )}
  
  res = matrix( NA, nrow=2, ncol=length(dir1) )
  
  res[1,] = abs(dir1 - dir2)
  res[2,] = abs(res[1,] - 360)
  
  diff = apply(res,2,min)
  
  return ( diff )
  
}



sortMat <- function (Mat, Sort, ascending=TRUE)
  #  by   Renaud Lancelot <lancelot at telecomplus.sn>a
  #  Sort matrix or dataframe 'Mat', by column(s) 'Sort'.
  #  e.g. sortmat.R(datafr, c(2,4,1))
{
  m <- do.call("order", as.data.frame(Mat[, Sort]))
  
  if (!ascending)
    m=rev(m)
  
  Mat[m, ]
}

WRF.compare.obs.wrf = function(wrf.raster.stub) {# function of compare weather observation station data and wrf output
  
  wrf.out.path    = "~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Feb2016/d03_HR/"
  obs.path        = "~/geolab_storage/data/Marcellus/Reprojection/WeatherStation/Marcellus2015/"
  wrf.raster.path = "~/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Comparisons-Feb2016/"
  obs.file.stub   = paste(obs.path,"qc_obs_used.d01.2015-05-DAY_TIME.0000HR",sep="")
  wrf.out.files   = dir(path=wrf.out.path,pattern="\\.nc$") #files names in WRF data dir
  
  res         = NULL
  for (f in 1:25) {
    
    
    day = 14
    
    if (f > 12) {
      day = 15
    }  
    
    date           = gsub("-",":",unlist(strsplit( unlist(strsplit(wrf.out.files[ f ],"_"))[4], "\\."))[1]) #the hour in the file name eg. when f=9, date 20:00:00
    obs.file       = gsub("TIME", date, obs.file.stub)
    obs.file       = gsub("DAY", day, obs.file)
    
    wrf.file       = paste(wrf.raster.path, wrf.raster.stub, f, "_1.tif",sep="")
    
    if ( file.exists(obs.file) && file.exists(wrf.file) ) {
      print(obs.file)
      print(wrf.file)
      obs        = cbind( WRF.qc.obs(obs.file), day = day, time = date)
      wrf.raster = raster(wrf.file)
      wrf        = extract(wrf.raster, cbind(obs$long, obs$lat) )
      
      res        = rbind(res, cbind(obs, WRF=wrf) )
    }
  }
  res$day_time   =with(res,paste0(day,time))
  res$sq_obs_wrf =(res$Temp_C-res$WRF) ** 2   #difference between obs and wrf in square
  return (res)
}



