fdrTest <-
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

	fdr<-array(0,dim=c(nlon,nlat,ntim))
	sig_FDR<-array(0,dim=c(nlon,nlat,ntim))
	p<-array(NA,dim=c(nlon*nlat))

	#put all p-values in 1D vector
	prob_1D<-(c(p_val))

	# sort vector increasing
	p_sort<-sort(prob_1D, decreasing = FALSE)

	# reject those local tests for which max[p(k)<=(siglev^(k/K)]
	for (k in 1:K){
	    if (p_sort[k]<=(siglev*(k/K))){
	       p[k]<-p_sort[k]
	    } else {
	      p[k]<-0.0
	    }
	}
	p_fdr<-max(p,na.rm=T)

	fdr[which(p_val<=p_fdr)] <- 1
	sig_FDR[which(fdr==1 & h_val==1)] <- 1

	sig_pts<-array(NA,dim=c(ntim))
	for (j in 1:ntim){
	    sig_pts[j]<-(sum(sig_FDR[,,j],na.rm=T))
	}

	method <- paste("False Discovery Rate for field significance")
	rval <- list(h.value=sig_FDR, p.value=p_val, field.sig = siglev,
	     	nr.sigpt=sig_pts, total.test=K, method=method, call=match.call())
	class(rval) <- "fdrFS"
	return(rval)
}
