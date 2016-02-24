# Calculate t statistic
# for modified t-test for small samples (ne<=30)
# Zwiers and von Storch, 1995, J. Clim, Eq. 16

calc_Tstat_small <- function(x,y,xAve,yAve){
source(paste("/home/z3441306/scripts/stat_tests_correlated_climdata/modTtest/calc_poolSampVar_3d.R",sep="/"))

dimsX <- dim(x)
xnyrs<-dimsX[3]

dimsY<-dim(y)
ynyrs<-dimsY[3]

if ( dimsX[2] != dimsY[2]){
   stop("Error: X and Y have different number of latitudes")
} else if ( dimsX[1] != dimsY[1]){
   stop("Error: X and Y have different number of longitudes")
} else{

	nlat<-dimsX[2]
	nlon<-dimsX[1]
}


#Pooled sample variance, Eq. 13
sp2<-calc_poolSampVar(x,y,xAve,yAve)

#t-statistic, Eq. 16
tval<-(yAve-xAve)/(sqrt(sp2)*(sqrt(1/xnyrs+1/ynyrs)))


return(tval)
}
