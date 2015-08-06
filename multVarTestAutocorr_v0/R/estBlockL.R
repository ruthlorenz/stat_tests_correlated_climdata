estBlockL <-
function(x,y){

dimsX <- dim(x)
dimsY<-dim(y)

if ( dimsX[3] != dimsY[3]){
   stop("Error: X and Y have different number of time")
} else if ( dimsX[2] != dimsY[2]){
   stop("Error: X and Y have different number of latitudes")
} else if ( dimsX[1] != dimsY[1]){
   stop("Error: X and Y have different number of longitudes")
} else{
        ntim<-dimsX[3]
        nlat<-dimsX[2]
        nlon<-dimsX[1]
}

findL <- function(l,length,b,ntim){
       	  while ( length < (ntim-l +1)^b ) {
      	  	l=l+1
      		length=l
	  }
	  return(l)
}

ACF_x=apply(x,c(1:2),function(z) acf(z,lag.max=1,type="correlation",plot=FALSE,na.action = na.pass))
ACF_y=apply(y,c(1:2),function(w) acf(w,lag.max=1,type="correlation",plot=FALSE,na.action = na.pass))

r1_x=apply(ACF_x,c(1,2),function(v) v[[1]]$acf[2])
r1_y=apply(ACF_y,c(1,2),function(u) u[[1]]$acf[2])

k <-1:(ntim-1)
V_x=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    V_x[lon,lat]=1+2*sum((1-k/ntim)*r1_x[lon,lat]^k,na.rm=T)
	}#end lon
} #end lat

V_prime_x=V_x*exp((2*V_x)/ntim)

b_x<-(2/3)*(1-1/V_prime_x)

l_x_2d=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    l_x <-1
	    length_x<-l_x
	    l_x_2d[lon,lat]=findL(l_x,length_x,b_x[lon,lat],ntim)

	}#end lon
} #end lat

L_x<-mean(l_x_2d,na.rm=T)

V_y=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    V_y[lon,lat]=1+2*sum((1-k/ntim)*r1_y[lon,lat]^k,na.rm=T)
	}#end lon
} #end lat

V_prime_y=V_y*exp((2*V_y)/ntim)

b_y<-(2/3)*(1-1/V_prime_y)

l_y_2d=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    l_y <-1
	    length_y<-l_y
	    l_y_2d[lon,lat]=findL(l_y,length_y,b_y[lon,lat],ntim)

	}#end lon
} #end lat

L_y<-mean(l_y_2d,na.rm=T)

L=floor((L_x+L_y)/2)

return(list(L=L,V_prime_x=V_prime_x,V_prime_y=V_prime_y))

}
