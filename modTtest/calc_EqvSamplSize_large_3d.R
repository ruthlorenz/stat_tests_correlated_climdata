# Calculate equivalent sample size 
# for Modified t-test for large samples (ne>=30)
# Zwiers and von Storch, 1995, J.Clim, Eq. 10

calc_EqvSamplSize_large <- function(x,y=NULL,r1){

dimsX<-dim(x)
xnyrs<-dimsX[3]

if (!is.null(y)){
   dimsY<-dim(y)
   ynyrs<-dimsY[3]

    if ( dimsX[2] != dimsY[2]){
       stop("Error: X and Y have different number of latitudes")
    } else if ( dimsX[1] != dimsY[1]){
       stop("Error: X and Y have different number of longitudes")
    }
}

nlat<-dimsX[2]
nlon<-dimsX[1]

xne_p<-array(NA,dim=c(nlon,nlat,1))
ne<-array(NA,dim=c(nlon,nlat,1))

if ((xnyrs>50) & (any(abs(r1)!=1.0))){
	xne<-xnyrs*((1-r1)/(1+r1))
} else {
	den<-array(0,dim=c(nlon,nlat,1))
	for (n in 1:xnyrs-1){
		if (!(all(is.na(x[,,n])))){
		den<-den+(1-n/xnyrs)*r1^n
		} 
	}
	xne<-xnyrs/(1+2*den)
}

if (!is.null(y)){
    yne_p<-array(NA,dim=c(nlon,nlat,1))
    if ((ynyrs>50) && (abs(r1)!=1.0)){
	    yne<-ynyrs*((1-r1)/(1+r1))
    } else {
	    den<-array(0,dim=c(nlon,nlat,1))
	    for (n in 1:ynyrs-1){
		    if (!(all(is.na(y[,,n])))){
		    den<-den+(1-n/ynyrs)*r1^n
		    }
	    }
	    yne<-ynyrs/(1+2*den)
    }
}
    for (lat in 1:nlat){
    	for (lon in 1:nlon){
	    if (all(is.na(x[lon,lat,]))) {
	    xne_p[lon,lat,1]<-NA
	    if (!is.null(y)){
	       yne_p[lon,lat,1]<-NA
	    }
	    ne[lon,lat,1]<-NA
	    }
	    else {
		if (r1[lon,lat,1]==0){
			xne_p[lon,lat,1]<-xnyrs
			if (!is.null(y)){
			   yne_p[lon,lat,1]<-ynyrs
			}
		} else {
	    		if (xne[lon,lat,1] <= 2){
	       			xne_p[lon,lat,1]<-2
	    		} else if ((xne[lon,lat,1]>2) && (xne[lon,lat,1]<=xnyrs)){
	       			xne_p[lon,lat,1]<-xne[lon,lat,1]
	   		} else {
	       			xne_p[lon,lat,1]<-xnyrs
	    		}
			if (!is.null(y)){
			    if (yne[lon,lat,1] <= 2){
				    yne_p[lon,lat,1]<-2
			    } else if ((yne[lon,lat,1]>2) && (yne[lon,lat,1]<=ynyrs)){
				    yne_p[lon,lat,1]<-yne[lon,lat,1]
			    } else {
				    yne_p[lon,lat,1]<-ynyrs
			    }
			}
		}
	}
	}# end lon
    }#end lat
if (!is.null(y)){
   ne<-(sqrt(1/xne_p+1/yne_p))
   #ne<-(1/sqrt(xne_p)+1/sqrt(yne_p))
   return(list(Xne_p=xne_p,Yne_p=yne_p,NE=ne))	
} else {
       ne<-xne_p
       return(ne)
}

}
