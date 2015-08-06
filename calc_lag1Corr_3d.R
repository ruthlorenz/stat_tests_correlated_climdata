# Calculate lag 1 correlation coefficient
# for modified t-test (large and small samples)
# Zwiers and von Storch, 1995, J.Clim, Eq. 14 (two sample case)

calc_lag1Corr <- function(x,y = NULL,xAve,yAve = NULL){

dimsX<-dim(x)
xnyrs<-dimsX[3]

if (!is.null(y)){
   dimsY<-dim(y)
   ynyrs<-dimsY[3]

    if (  dimsX[2] != dimsY[2]){
       stop("Error: X and Y have different number of latitudes")
    } else if ( dimsX[1] != dimsY[1]){
       stop("Error: X and Y have different number of longitudes")
    }
}
nlat<-dimsX[2]
nlon<-dimsX[1]

xsumsq<-array(0,dim=c(nlon,nlat,1))
xsumlag<-array(0,dim=c(nlon,nlat,1))
if (!is.null(y)){
   ysumsq<-array(0,dim=c(nlon,nlat,1))
   ysumlag<-array(0,dim=c(nlon,nlat,1))
}
r1<-array(NA,dim=c(nlon,nlat,1))

	    #Sum of Squares
	    for (yr in 1:xnyrs){
	    xsumsq<-xsumsq+
		((x[,,yr, drop = FALSE]-xAve)^2)
	    }
	    if (!is.null(y)){
		for (yr in 1:ynyrs){
		ysumsq<-ysumsq+
		    ((y[,,yr, drop = FALSE]-yAve)^2)
		}
	    xysumsq<-xsumsq+ysumsq
	    }

	    for (yr in 2:xnyrs){
	    xsumlag<-xsumlag+
		((x[,,yr, drop = FALSE]-xAve)*
		(x[,,yr-1, drop = FALSE]-xAve))
	    }
	    if (!is.null(y)){
		for (yr in 2:ynyrs){
		ysumlag<-ysumlag+
		    ((y[,,yr, drop = FALSE]-yAve)*
		    (y[,,yr-1, drop = FALSE]-yAve))
		}
	    xysumlag<-xsumlag+ysumlag
	    r1<-xysumlag/xysumsq
	    } else {
	    r1<-xsumlag/xsumsq
	    }
return(r1)

}
