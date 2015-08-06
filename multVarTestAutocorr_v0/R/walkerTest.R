walkerTest <-
function(p_val, siglev=0.05, ...){

        dims_p <- dim(p_val)
        nlon<-dims_p[1]
        nlat<-dims_p[2]
	if (!is.na(dims_p[3])){
           ntim<-dims_p[3]
	} else {
	  ntim<-1
	  tmp<-p_val
	  p_val<-array(NA,dim=c(nlon,nlat,ntim))
	  p_val[,,1]<-tmp
	}

	h_val<-array(NA,dim=c(nlon,nlat,ntim))

	for (t in 1:ntim){
	    for (lat in 1:nlat){
		for (lon in 1:nlon){
			if (is.na(p_val[lon,lat,t])){
				h_val[lon,lat,t]<-NA
			} else if (p_val[lon,lat,t] < siglev){
				h_val[lon,lat,t]<-1
			} else {
				h_val[lon,lat,t]<-0
			}
		}
	    }
	}

	   K<-sum(!is.na(p_val[,,1]))


                sig_walk<-array(NA,dim=c(nlon,nlat,ntim))
                p_min<-array(NA,dim=c(ntim))

                #walker criterium (Wilks, 2006)
                p_walk<-1-(1-siglev)^(1/K)

                # test if minimum p-value per season is larger than p_walk -> not significant
                for (j in 1:ntim){
                    p_min[j]<-min(p_val[,,j],na.rm=T)

                    if (p_min[j]>p_walk){
                       sig_walk[,,j]<-0
                    }  else {
                       sig_walk[,,j]<-h_val[,,j]
                    }
                }

		sig_pts<-array(NA,dim=c(ntim))
		for (j in 1:ntim){
		       sig_pts[j]<-(sum(sig_walk[,,j],na.rm=T))
		}

	method <- paste("Walkers test for field significance")
	rval <- list(h.value=sig_walk, p.value=p_val, field.sig = siglev,
	     	nr.sigpt=sig_pts, total.test=K, method=method, call=match.call())
	class(rval) <- "walkFS"
	return(rval)
}
