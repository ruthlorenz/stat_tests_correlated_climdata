# calculate block length for moving blocks bootstrapping
# Wilks, 1997, J.Climate
# L=(n-L+1)^((2/3)(1-n'/n))
# "largest integer no greater than L", evaluated iteratively
# n'~=n((1-r1)/(1+r1))
# r1=lag-1 autocorrelation coefficient

estBlockL <- function(x,y=NULL){

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
	  #for (yr in 1:ntim){
          #    xsumsq<-xsumsq+((x[yr]-xAve)^2)
          #}
	  xsumlag=0
          for (yr in 1:(ntim-1)){
              xsumlag<-xsumlag+((x[yr]-xAve)*(x[yr+1]-xAve))
          }
       	  R1 <- xsumlag/xsumsq
	  return(R1)
}
ACF_x=apply(x,c(1:2),function(z) acf(z,lag.max=1,type="correlation",plot=FALSE,na.action = na.pass))
r1_x=apply(ACF_x,c(1,2),function(v) v[[1]]$acf[2])
#r1_x=apply(x,c(1:2),calcR1)

k <-1:(ntim-1)
V_x=array(NA,dim=c(nlon,nlat))
for (lat in 1:nlat){
        for (lon in 1:nlon){
	    V_x[lon,lat]=1+2*sum((1-k/ntim)*r1_x[lon,lat]^k,na.rm=T)
	    #if (lat==1 & lon==1) print(r1_x[lon,lat]^k)
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

if( !is.null(y) ) {
    ACF_y=apply(y,c(1:2),function(w) acf(w,lag.max=1,type="correlation",plot=FALSE,na.action = na.pass))
    r1_y=apply(ACF_y,c(1,2),function(u) u[[1]]$acf[2])
    #r1_y=apply(y,c(1:2),calcR1)

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
} else {
  V_prime_y = NULL
  L=floor(L_x)
}

return(list(L=L,V_prime_x=V_prime_x,V_prime_y=V_prime_y))

}

