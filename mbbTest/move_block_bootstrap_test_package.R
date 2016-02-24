# Moving blocks bootstrapping after Wilks 1997, J.Clim
# test for 3d data (lon,lat,time), taking into account
# autocorrelation in time and multiplicity, correlation in space

mbbTest <- function(x, ...) UseMethod("mbbTest")

mbbTest.default <- function(x, y = NULL, mu = 0, model = "AR1",
		L = NULL, siglev=0.05, alpha=0.05, nb=1000, verbose = FALSE,
		timevals=NULL, na.rm = TRUE, ...){

source("/home/z3441306/scripts/plot_scripts/R_scripts/tests_for_autocorr_fields/estimate_block_length.R")
source("/home/z3441306/scripts/plot_scripts/R_scripts/tests_for_autocorr_fields/estimate_block_length2.R")
source("/home/z3441306/scripts/plot_scripts/R_scripts/tests_for_autocorr_fields/move_blocks_bootstrap_wilks.R")

#------ Find and check dimensions ------#
dimsX <- dim(x)
if( !is.null(y) ) {
    dimsY<-dim(y)

    if (all(dimsX!=dimsY)){
       stop("wrong dimensions, dim of X must be the same as dim of Y")
       if ( dimsX[3] != dimsY[3]){
	  stop("Error: X and Y have different number of time")
       } else if ( dimsX[2] != dimsY[2]){
	 stop("Error: X and Y have different number of latitudes")
       } else if ( dimsX[1] != dimsY[1]){
	 stop("Error: X and Y have different number of longitudes")
       }
    }
}
nyears<-dimsX[3]
nlat<-dimsX[2]
nlon<-dimsX[1]  

if(!is.null(timevals)){
	if(nyears!=timevals){
		stop("Error: wrong dimension, the time in X and Y does not
			     correspond to timevals")
	}
}

#------ Initialize output list ------#
out <- list()
if( !is.null(y) ) {
    data_name <- c(as.character(substitute(x)),as.character(substitute(y)))
    names(data_name) <- c("Control","Experiment")
    out$data.name <- data_name
} else {
    data_name <- c(as.character(substitute(x)))
    names(data_name) <- c("Mean ")
    out$data.name <- data_name
}
#------ Calculate time average ------#
if (any(is.na(x))) {
    if (na.rm) x_avg <- apply(x,c(1:2),mean,na.rm=TRUE)
    else stop("'x' contains missing values")
} else {
  x_avg <- apply(x,c(1:2),mean)
}

if( !is.null(y) ) {
    if (any(is.na(y))) {
       if (na.rm) y_avg <- apply(y,c(1:2),mean,na.rm=TRUE)
       else stop("'y' contains missing values")
    } else {
           y_avg <- apply(y,c(1:2),mean)
    }
    anom_avg <- x_avg-y_avg
} else {
    anom_avg <- x_avg - mu
}

#------ Estimate block length ------#
if ( model == "AR1" ){
   est_L <- estBlockL(x,y)
} else if ( model == "AR2" ){
   est_L <- estBlockL2(x,y)
} else {
   stop("'model' can only be 'AR1' or 'AR2'")
}
if (is.null(L)){
   L <- est_L$L
}
if(verbose) print(L)

V_prime_x <-est_L$V_prime_x

if( !is.null(y) ) {
    V_prime_y <- est_L$V_prime_y
}

#------ Calculate local test statistic d=anom_seas/sigma_s ------#
x_var <- apply(x,c(1:2),var)
if( !is.null(y) ) {
    y_var <- apply(y,c(1:2),var)
    sigma_s <- (V_prime_x*x_var/nyears+V_prime_y*y_var/nyears)^(1/2)
} else {
    sigma_s <- (V_prime_x*x_var/nyears)^(1/2)
}

d <- anom_avg/sigma_s
D <- max(abs(d))

#------ Perform moving blocks bootstrapping ------#
if( !is.null(y) ) {
    move_boot <- moveBlocksBoot(x=x,y=y,L=L,nb=nb,verbose=verbose)
    method <- paste("Two-sample moving blocks bootstrapping ",model,sep="")
} else {
       move_boot <- moveBlocksBoot(x=x,mu=mu,L=L,nb=nb,verbose=verbose)
       method <- paste("One sample moving blocks bootstrapping ",model,sep="")
}

d_star <- move_boot$d_star
D_star <- move_boot$D_star

#------ Compare d-statistic to bootstrapped distribution of d*-statistic ------#
# for local significance
# test rejects if H0  #(d*>=d)/(nb+1) <=alpha/2 for upper tail and
#                     #(d*<=d)/(nb+1) <=alpha/2 for lower tail
 signif<-array(NA,dim=c(nlon,nlat))
for (lon in 1:nlon){
	for (lat in 1:nlat){
        	I_u <- sum(d_star[lon,lat,1:nb]>=d[lon,lat],na.rm=TRUE)
                I_l <- sum(d_star[lon,lat,1:nb]<=d[lon,lat],na.rm=TRUE)

		if (all(is.na(x[lon,lat,]))){
		   signif[lon,lat] <- NA
		} else {
		    if (I_u/(nb+1)<= (as.numeric(siglev)/2)) {
			    signif[lon,lat] <- 1
		    } else if (I_l/(nb+1)<= (as.numeric(siglev)/2)) {
			    signif[lon,lat] <- 1
		    } else {
			    signif[lon,lat] <- 0}
		}
	}
}
print(sum(signif,na.rm=T))

loc_sig_coverage <- sum(signif,na.rm=TRUE)/sum(!is.na(x[,,1])) 

#------ Compare D-statistic to bootstrapped distribution of D*-statistic ------#
#------ for field significance

I_u=sum(D_star[1:nb]>=D,na.rm=TRUE)
I_l=sum(D_star[1:nb]<=D,na.rm=TRUE)

if (I_u/(nb+1)<= (as.numeric(siglev)/2)) {
	fs = 1
} else if (I_l/(nb+1)<= (as.numeric(siglev)/2)) {
	fs = 1
} else {
	fs = 0}

out$estimate <- anom_avg
out$loc.sig <- siglev
out$field.sig <- alpha
out$h.value <- signif
out$h.valueFS <- fs
out$nr.sigpt <- sum(signif,na.rm=T)
out$total.test <- sum(!is.na(x[,,1]))
out$loc.sig.coverage <- loc_sig_coverage
out$FS <- if(fs==1) TRUE else FALSE
out$nboot <- nb
out$method <- method
out$call <- match.call()

class(out) <- "mbbtest"
return(out)
} #end of 'mbbTest'

print.mbbtest <- function(x, ...){
	    cat("Call:\n")
	    print(x$call)
	    cat("\nNumber of Significant points:", x$nr.sigpt,"out of ",
	    		   x$total.test," tests.","\n")
} 

summary.mbbtest <- function(object, ...){
   cat("\n")
   msg <- paste("Results for ",object$method," ",object$data.name[2],
        " compared against ", object$data.name[1], sep="")
   print(msg)
   cat("\n")
   cat("Local significance level: ", object$loc.sig, "\n")
   cat("Field significance level: ", object$field.sig, "\n")
   cat("Number of locally significant tests: ",object$nr.sigpt, "\n")
   cat("Total number of tests: ", object$total.test, "\n")
   cat("Local significance coverage: ", object$loc.sig.coverage, "\n")
   cat("Result of field significance test: ", object$FS, "\n")
   invisible()
}
