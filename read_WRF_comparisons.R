
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


library(raster)
library(stringi)

dir="/Volumes/geolab_storage/data/Marcellus/Reprojection/WRF-Output-Comparisons-Feb2016"
setwd(dir)

str_name <- list.files(pattern = "diff.HR-HR.RESHIFT.WindSpeed_.*\\.tif$", recursive=TRUE)

raster <- raster(str_name[3])
n<- na.omit(getValues(raster))
col.name = stri_sub(str_name[3],25,-5)
dat = data.frame(data=n, names = c(rep(col.name, length(n))))

order = c(8,13,5,10,15,1,6,11,2,7,12,4,9,14)
for (i in order) {
  str_name[i]
  raster <- raster(str_name[i])
  n<- na.omit(getValues(raster))
  col.name = stri_sub(str_name[i],25,-5)
  dat = rbind(dat,data.frame(data=n,names=col.name))
  #hist(n,main = i,breaks="FD")
}

boxplot(data ~ names, data =dat, 
        outline = FALSE,
        ylab ="Wind Speed. Diff. (m/s)", 
        xlab ="",
        las =2, 
        col = c("red","sienna","palevioletred1"),
        at = c(1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19),
        names = c("5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM", "5.14 12PM", "5.14 23PM", "5.15 12PM"),
        par(mar = c(12, 5, 4, 2)+ 0.1)
        )

rect(0,-2,4,2,col="lightgrey",border=NA)
rect(8,-2,12,2,col="lightgrey",border=NA)
rect(16,-2,20,2,col="lightgrey",border=NA)
lines(c(0,100),c(0,0),lty=2)

title("HR-HR_Reshift")
mtext("Layer 1        Layer 7       Layer 13      Layer 19       Layer 25",3)

boxplot(data ~ names, data =dat, 
        outline = FALSE,
        ylab ="Wind Speed. Diff. (m/s)", 
        xlab ="",
        las =2, 
        col = c("red","sienna","palevioletred1"),
        at = c(1,2,3, 5,6,7, 9,10,11, 13,14,15, 17,18,19),
        names = c("5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM","5.14 12PM", "5.14 23PM", "5.15 12PM", "5.14 12PM", "5.14 23PM", "5.15 12PM"),
        par(mar = c(12, 5, 4, 2)+ 0.1),add=T)


#mtext("Difference",1,line=9)

