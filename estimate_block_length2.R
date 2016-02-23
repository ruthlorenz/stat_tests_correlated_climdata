# calculate block length for moving blocks bootstrapping AR(2) processs
# Wilks, 1997, J.Climate
# L=(n-L+1)^((2/3)(1-1sqrt(4V')))
# "largest integer no greater than L", evaluated iteratively
# r1=lag-1 autocorrelation coefficient, Eq. 6
# r2=lag-2 autocorrelation coefficient. Eq. 21

estBlockL2 <- function(x, y = NULL){

dimsX <- dim(x)
if( !is.null(y) ) {
    dimsY<-dim(y)
}
ntim<-dimsX[3]
nlat<-dimsX[2]
nlon<-dimsX[1]

findL <- function(l,length,b,ntim){
       	  while ( length < (ntim-l +1)^b ) {
      	  	l=l+1
      		length=l
	  }
	  return(l)
}
calcR1 <- function(x){
       	  ntim <- length(x)
	  xAve <- mean(x,na.rm=TRUE)
	  xsumsq=sum((x[1:ntim]-xAve)^2,na.rm=TRUE)
	  xsumlag=0
          for (yr in 1:(ntim-1)){
              xsumlag<-xsumlag+((x[yr]-xAve)*(x[yr+1]-xAve))
          }
       	  R1 <- xsumlag/xsumsq
	  return(R1)
}
calcR2 <- function(x){
       	  ntim <- length(x)
	  xAve <- mean(x,na.rm=TRUE)
	  xAve_m <- mean(x[1:(ntim-2)],na.rm=TRUE)
	  xAve_p <- mean(x[3:ntim],na.rm=TRUE)
	  xsumsq1=sum((x[1:(ntim-2)]-xAve_m)^2,na.rm=TRUE)
	  xsumsq2=sum((x[3:ntim]-xAve_p)^2,na.rm=TRUE)
	  xsumlag=0
          for (yr in 1:(ntim-2)){
	      xsumlag<-xsumlag+((x[yr]-xAve_m)*(x[yr+2]-xAve_p))
          }
       	  R2 <- xsumlag/(xsumsq1*xsumsq2)^(1/2)
	  return(R2)
}
rk_x=array(NA,dim=c(nlon,nlat,ntim))

rk_x[,,1]<- apply(x,c(1:2),calcR1)
rk_x[,,2]<- apply(x,c(1:2),calcR2)

ph1_x <- rk_x[,,1]*(1-rk_x[,,2])/(1-rk_x[,,1]^2)
ph2_x <- (rk_x[,,2]-rk_x[,,1]^2)/(1-rk_x[,,1]^2)

V_x=array(NA,dim=c(nlon,nlat))

k <-1:(ntim-1)
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    for (ka in 3:(ntim-1)){
	    	rk_x[lon,lat,ka]=ph1_x[lon,lat]*rk_x[lon,lat,ka-1]+ph2_x[lon,lat]*rk_x[lon,lat,ka-2]
	    }
	    V_x[lon,lat]=1+2*sum((1-k/ntim)*rk_x[lon,lat,k],na.rm=T)
	}#end lon
} #end lat

V_prime_x=V_x*exp((3*V_x)/ntim)

b_x<-(2/3)*(1-1*sqrt(4*V_prime_x))

l_x_2d=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    l_x <-1
	    length_x<-l_x
	    l_x_2d[lon,lat]=findL(l_x,length_x,b_x[lon,lat],ntim)
	}#end lon
} #end lat

L_x<-mean(l_x_2d,na.rm=TRUE)

if( !is.null(y) ) {
    rk_y=array(NA,dim=c(nlon,nlat,ntim))
    rk_y[,,1]=apply(y,c(1:2),calcR1)
    rk_y[,,2]=apply(y,c(1:2),calcR2)
    ph1_y=rk_y[,,1]*(1-rk_y[,,2])/(1-rk_y[,,1]^2)
    ph2_y=(rk_y[,,2]-rk_y[,,1]^2)/(1-rk_y[,,1]^2)

    V_y=array(NA,dim=c(nlon,nlat))
    
    for (lat in 1:nlat){
	    for (lon in 1:nlon){
	    	for (ka in 3:(ntim-1)){
		    rk_y[lon,lat,ka]=ph1_y[lon,lat]*rk_y[lon,lat,ka-1]+
					ph2_y[lon,lat]*rk_y[lon,lat,ka-2]
		}
		V_y[lon,lat]=1+2*sum((1-k/ntim)*rk_y[lon,lat,k],na.rm=T)
	    }#end lon
    } #end lat

    V_prime_y=V_y*exp((3*V_y)/ntim)

    b_y<-(2/3)*(1-1*sqrt(4*V_prime_y))

    l_y_2d=array(NA,dim=c(nlon,nlat))
    for (lat in 1:nlat){
	    for (lon in 1:nlon){
		l_y <-1
		length_y<-l_y
		if (is.na(b_y[lon,lat])){
		   l_y_2d[lon,lat]=NA
		} else {
		  l_y_2d[lon,lat]=findL(l_y,length_y,b_y[lon,lat],ntim)
		}
	    }#end lon
    } #end lat

    L_y<-mean(l_y_2d,na.rm=TRUE)

    L=floor((L_x+L_y)/2)
} else {
  V_prime_y = NULL
  L=floor(L_x)
}

return(list(L=L,V_prime_x=V_prime_x,V_prime_y=V_prime_y))

}

