moveBlocksBoot <-
function (x, y, L=1, nb=1000, verbose=FALSE) {

dimsX <- dim(x)
dimsY<-dim(y)

if ( dimsX[3] != dimsY[3]){
   stop("Error: X and Y have different number of months")
} else if ( dimsX[2] != dimsY[2]){
   stop("Error: X and Y have different number of latitudes")
} else if ( dimsX[1] != dimsY[1]){
   stop("Error: X and Y have different number of longitudes")
} else{
        nyears<-dimsX[3]
        nlat<-dimsX[2]
        nlon<-dimsX[1]
}

#------ number of blocks: ------#
b=ceiling(nyears/L)

x_mean = apply(x,c(1:2),mean,na.rm=T)
y_mean = apply(y,c(1:2),mean,na.rm=T)
anom = x_mean-y_mean

exp_bs <- array(NA,dim=c(nlon,nlat,nyears,nb))
ctl_bs <- array(NA,dim=c(nlon,nlat,nyears,nb))
exp_mean_bs <- array(NA,dim=c(nlon,nlat,nb))
ctl_mean_bs <- array(NA,dim=c(nlon,nlat,nb))
boot_mean <- array(NA,dim=c(nlon,nlat,nb))
x_exp_star_i_jack <- array(NA,dim=c(nlon,nlat,b))
x_ctl_star_i_jack <- array(NA,dim=c(nlon,nlat,b))
d_star <- array(NA,dim=c(nlon,nlat,nb))
D_star <- array(NA,dim=c(nb))
sigma_s_star<- array(NA,dim=c(nlon,nlat,nb))

for (iboot in 1:nb) {
    #------ First bootstrap exp ------#
    endpoints <- sample(L:(nyears), size=b, replace=T)
    #------ Find vector of indexes of the block samples ------#
    block_indices <- lapply(endpoints, function(i) seq(max(0, i-L+1,na.rm=T), i))
    blocks <- lapply(block_indices, function(i) x[,,i])  # select the blocks in the order of the indices
    tmp_bs <- abind::abind(blocks, along=3)
    exp_bs[,,,iboot] <- tmp_bs[,,1:nyears]

    #------ x_exp_star_i for jacknife-after-bootstrap variance ------#
    x_exp_star_i_list=lapply(blocks,apply,c(1,2),sum,na.rm=T)
    x_exp_star_i=abind::abind(x_exp_star_i_list, along=3)

    #------ Second bootstrap ctl------#
    endpoints <- sample(L:(nyears), size=b, replace=T)
    #------ Find vector of indexes of the block samples ------#
    block_indices <- lapply(endpoints, function(i) seq(max(0, i-L+1,na.rm=T), i))
    blocks <- lapply(block_indices, function(i) y[,,i])  # select the blocks in the order of the indices
    tmp_bs <- abind::abind(blocks, along=3)
    ctl_bs[,,,iboot] <- tmp_bs[,,1:nyears]

    #------ x_ctl_star_i for jacknife-after-bootstrap variance ------#
    x_ctl_star_i_list=lapply(blocks,apply,c(1,2),sum,na.rm=T)
    x_ctl_star_i=abind::abind(x_ctl_star_i_list, along=3)

    #------ Estimate sampling variability through jacknife variance estimates for exp_bs and ctl_bs ------#
    for (i in 1:b){
    	x_exp_star_i_jack[,,i]=apply(x_exp_star_i,c(1,2),function(j) (sum(j,na.rm=T)-j[i])/(b*L-L))
	x_ctl_star_i_jack[,,i]=apply(x_ctl_star_i,c(1,2),function(k) (sum(k,na.rm=T)-k[i])/(b*L-L))
    }

    x_exp_star_dot<-array(NA,dim=c(nlon,nlat))
    x_ctl_star_dot<-array(NA,dim=c(nlon,nlat))
    x_exp_star_dot<-apply(x_exp_star_i_jack,c(1,2),mean,na.rm=T)
    x_ctl_star_dot<-apply(x_ctl_star_i_jack,c(1,2),mean,na.rm=T)

    sigma_exp_2<-array(NA,dim=c(nlon,nlat))
    sigma_ctl_2<-array(NA,dim=c(nlon,nlat))
    for (lat in 1:nlat){
    	for (lon in 1:nlon){
	    sigma_exp_2[lon,lat]=(b*L/nyears)*((b-1)/b)*sum((x_exp_star_i_jack[lon,lat,1:b]-x_exp_star_dot[lon,lat])^2,na.rm=T)
	    sigma_ctl_2[lon,lat]=(b*L/nyears)*((b-1)/b)*sum((x_ctl_star_i_jack[lon,lat,1:b]-x_ctl_star_dot[lon,lat])^2,na.rm=T)
	}
    }
    sigma_s_star[,,iboot]=(sigma_exp_2+sigma_ctl_2)^(1/2)
    
    exp_mean_bs[,,iboot]=apply(exp_bs[,,,iboot],c(1:2),mean,na.rm=T)
    ctl_mean_bs[,,iboot]=apply(ctl_bs[,,,iboot],c(1:2),mean,na.rm=T)
    boot_mean[,,iboot]=exp_mean_bs[,,iboot]-ctl_mean_bs[,,iboot]

    #------ Calculate d_star for each bootstrap ------#
    d_star[,,iboot]=(boot_mean[,,iboot]-anom)/sigma_s_star[,,iboot]
    D_star[iboot]=max(abs(d_star[,,iboot]),na.rm=T)

#------ If verbose report every 100 bootstrap passes ------#
if(verbose & iboot%%100==0) print(iboot)
}

return(list(d_star=d_star,D_star=D_star))
}
