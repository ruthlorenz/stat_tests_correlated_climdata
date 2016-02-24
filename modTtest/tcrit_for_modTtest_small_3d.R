# evaluate tcrit for modified t-test for small samples based on
# Zwiers and von Storch, 1995, J. Clim 

lookupTcrit <- function(r1,siglev,nobs,y=NULL){

source(paste("/home/z3441306/scripts/stat_tests_correlated_climdata/modTtest/create_LookupTabs_modTtest.R",sep="/"))

a<-seq(-0.35,0.95,length=27)
b<-a-0.025
c<-a+0.025

alpha<-array(NA,dim=c(27,3))
alpha[,1]<-a
alpha[,2]<-b
alpha[,3]<-c

stat_val<-lookupTabs(nobs)
if (is.null(y)) {
    if (isTRUE(all.equal(siglev,0.10))){
      tab<-stat_val[,1]
    } else if (isTRUE(all.equal(siglev,0.05))){
      tab<-stat_val[,2]
    } else if (isTRUE(all.equal(siglev, 0.025))){
      tab<-stat_val[,3]
    } else if (isTRUE(all.equal(siglev, 0.01))){
      tab<-stat_val[,4]
    } else if (isTRUE(all.equal(siglev, 0.005))){
      tab<-stat_val[,5]
    } else stop("significance level not found in lookup table, exiting")

} else {

    if (isTRUE(all.equal(siglev, 0.20))){
      tab<-stat_val[,1]
    } else if (isTRUE(all.equal(siglev, 0.10))){
      tab<-stat_val[,2]
    } else if (isTRUE(all.equal(siglev, 0.05))){
      tab<-stat_val[,3]
    } else if (isTRUE(all.equal(siglev, 0.02))){
      tab<-stat_val[,4]
    } else if (isTRUE(all.equal(siglev, 0.01))){
      tab<-stat_val[,5]
    }
}

dims<-dim(r1)
nlat<-dims[2]
mlon<-dims[1]

tcrit=array(NA,dim=c(mlon,nlat,1))

    for (lat in 1:nlat){
    	for (lon in 1:mlon){
	    r=r1[lon,lat,1]
	    if (r<min(alpha[,2]) || r>max(alpha[,3])){
    	       tcrit[lon,lat,1]<-NA
    	       print('Error in tcrit_for_modTtest_small_3d.R, no tcrit found')
	    } else if (min(which(alpha[,3]>=r),na.rm=T) && max(which(alpha[,2]<r),na.rm=T)){
	    	 tindx<-min(which(alpha[,3]>=r),na.rm=T)
		 tcrit[lon,lat,1]<-tab[tindx]
	    } else {
		tcrit[lon,lat,1]<-NA
		print('Error in tcrit_for_modTtest_small_3d.R, no tcrit found')
	    }
	}
    }

return(tcrit)
}
