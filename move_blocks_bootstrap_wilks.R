# function for calculating moving blocks bootstrapping
# after Wilks, 1997, J. Clim
# bootstrap two samples separately
# take into account multiplicity, spatial correlation of data
# x,y: timesseries of exp and ctl
# L: block length
# nb: number of bootstraps (iboot)
# b: number of blocks

moveBlocksBoot <- function (x, y = NULL, mu, L, nb, verbose) {

require(abind)
if(require(abind)){
    #print("abind is loaded correctly")
} else {
    #print("trying to install abind")
    install.packages("abind",lib="/home/z3441306/scripts/plot_scripts/R_scripts/libraries_zip/abind_1.4-0.tar.gz", repos = NULL, type="source")
    if(require(abind)){
        #print("abind installed and loaded")
    } else {
        stop("could not install abind")
    }
}

dimsX <- dim(x)
if( !is.null(y) ) {
    dimsY<-dim(y)
}

nyears<-dimsX[3]
nlat<-dimsX[2]
nlon<-dimsX[1]

#------ number of blocks: ------#
b=ceiling(nyears/L)

x_mean = apply(x,c(1:2),mean,na.rm=T)
if( !is.null(y) ) {
    y_mean = apply(y,c(1:2),mean,na.rm=T)
    anom = x_mean-y_mean
} else {
  anom = x_mean-mu
}

x_bs <- array(NA,dim=c(nlon,nlat,nyears,nb))
x_mean_bs <- array(NA,dim=c(nlon,nlat,nb))
boot_mean <- array(NA,dim=c(nlon,nlat,nb))
x_star_i_jack <- array(NA,dim=c(nlon,nlat,b))
d_star <- array(NA,dim=c(nlon,nlat,nb))
D_star <- array(NA,dim=c(nb))
sigma_s_star<- array(NA,dim=c(nlon,nlat,nb))
if( !is.null(y) ) {
    y_bs <- array(NA,dim=c(nlon,nlat,nyears,nb))
    y_mean_bs <- array(NA,dim=c(nlon,nlat,nb))
    y_star_i_jack <- array(NA,dim=c(nlon,nlat,b))
}

for (iboot in 1:nb) {
    #------ First bootstrap exp ------#
    endpoints <- sample(L:(nyears), size=b, replace=T)
    #------ Find vector of indexes of the block samples ------#
    block_indices <- lapply(endpoints, function(i) seq(max(0, i-L+1,na.rm=T), i))
    blocks <- lapply(block_indices, function(i) x[,,i])  # select the blocks in the order of the indices
    tmp1_bs <- abind(blocks, along=3)
    x_bs[,,,iboot] <- tmp1_bs[,,1:nyears]

    #------ x_star_i for jacknife-after-bootstrap variance ------#
    x_star_i_list=lapply(blocks,apply,c(1,2),sum,na.rm=T)
    x_star_i=abind(x_star_i_list, along=3)

    if( !is.null(y) ) {
	#------ Second bootstrap ctl------#
	endpoints <- sample(L:(nyears), size=b, replace=T)
	#------ Find vector of indexes of the block samples ------#
	block_indices <- lapply(endpoints, function(i) seq(max(0, i-L+1,na.rm=T), i))
	blocks <- lapply(block_indices, function(i) y[,,i])  # select the blocks in the order of the indices
	tmp2_bs <- abind(blocks, along=3)
	y_bs[,,,iboot] <- tmp2_bs[,,1:nyears]

	#------ y_star_i for jacknife-after-bootstrap variance ------#
	y_star_i_list=lapply(blocks,apply,c(1,2),sum,na.rm=T)
	y_star_i=abind(y_star_i_list, along=3)
    }
    #------ Estimate sampling variability through jacknife variance estimates for x_bs and y_bs ------#
    for (i in 1:b){
    	x_star_i_jack[,,i]=apply(x_star_i,c(1,2),function(j) (sum(j,na.rm=T)-j[i])/(b*L-L))
	if( !is.null(y) ) {
	    y_star_i_jack[,,i]=apply(y_star_i,c(1,2),function(k) (sum(k,na.rm=T)-k[i])/(b*L-L))
	}
    }

    x_star_dot<-array(NA,dim=c(nlon,nlat))
    x_star_dot<-apply(x_star_i_jack,c(1,2),mean,na.rm=T)
    sigma_x_2<-array(NA,dim=c(nlon,nlat))
    if( !is.null(y) ) {
    	y_star_dot<-array(NA,dim=c(nlon,nlat))
    	y_star_dot<-apply(y_star_i_jack,c(1,2),mean,na.rm=T)
	sigma_y_2<-array(NA,dim=c(nlon,nlat))
    }


    for (lat in 1:nlat){
    	for (lon in 1:nlon){
	    sigma_x_2[lon,lat]=(b*L/nyears)*((b-1)/b)*
				sum((x_star_i_jack[lon,lat,1:b]-x_star_dot[lon,lat])^2,na.rm=T)
	    if( !is.null(y) ) {
	    	sigma_y_2[lon,lat]=(b*L/nyears)*((b-1)/b)*
				    sum((y_star_i_jack[lon,lat,1:b]-y_star_dot[lon,lat])^2,na.rm=T)
	    }
	}
    }
    if( !is.null(y) ) {
    	sigma_s_star[,,iboot]=(sigma_x_2+sigma_y_2)^(1/2)
    } else {
      	   sigma_s_star[,,iboot]=(sigma_x_2)^(1/2)
    }
    x_mean_bs[,,iboot]=apply(x_bs[,,,iboot],c(1:2),mean,na.rm=T)
    if( !is.null(y) ) {
    	y_mean_bs[,,iboot]=apply(y_bs[,,,iboot],c(1:2),mean,na.rm=T)
    	boot_mean[,,iboot]=x_mean_bs[,,iboot]-y_mean_bs[,,iboot]
    } else {
      	   boot_mean[,,iboot]=x_mean_bs[,,iboot]-mu
    }
    #------ Calculate d_star for each bootstrap ------#
    d_star[,,iboot]=(boot_mean[,,iboot]-anom)/sigma_s_star[,,iboot]
    D_star[iboot]=max(abs(d_star[,,iboot]),na.rm=T)

#------ If verbose report every 100 bootstrap passes ------#
if(verbose & iboot%%100==0) print(iboot)
}

return(list(d_star=d_star,D_star=D_star))
}
