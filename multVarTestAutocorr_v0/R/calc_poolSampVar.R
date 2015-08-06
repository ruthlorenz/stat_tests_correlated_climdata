calc_poolSampVar <-
function(x,y,xAve,yAve){

dimsX<-dim(x)
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

xsumsq<-array(0,dim=c(nlon,nlat,1))
ysumsq<-array(0,dim=c(nlon,nlat,1))

for (yr in 1:xnyrs){
    xsumsq<-xsumsq+((x[,,yr, drop = FALSE]-xAve)^2)
}
for (yr in 1:ynyrs){
    ysumsq<-ysumsq+((y[,,yr, drop = FALSE]-yAve)^2)
}	    
sp2<-((xnyrs-1)/xnyrs*xsumsq+(ynyrs-1)/ynyrs*ysumsq)/(xnyrs+ynyrs-2)

return(sp2)
}
