# stat_tests_correlated_climdata
Statistical tests for temporally autocorrelated and/or spatially correlated data. Implemented for climate model output (3 dimesions: time, lat, lon).

Functions for modified t-test (effective sample size t-test), two field significance tests (Walkers test and False Discovery Rate, can be used with modified t-test) and moving blocks bootstrap. The moving blocks bootstrap also includes a test for field significance.

This package contains four main functions, modTtest3d, mbbTest, walkerTest and fdrTest. All functions are designed for spatial data (longitude x latitude). modTtest3d and mbbTest need multiple timesteps while walkerTest and fdrTest can have multiple time steps but do not need to also work with longitude x latitide data only.

References:
Zwiers and von Storch, 1995, Taking serial correlation into account in tests of the mean, J. Clim, 8, p. 336--351.
Wilks, D.S., 1997, Resampling Hypothesis Tests for Autocorrelated Fields, J. Clim., 10, p.65--82.
Wilks, D.S., 2006, On "Field Significance" and the False Discovery Rate, J. Appl. Meteorol. Climatol., 45, p. 1181--1189.

Example:
x<-array(NA,dim=c(20,50,31))
y<-array(NA,dim=c(20,50,31))

for (lon in 1:20){
        for (lat in 1:50){
                #create timeseries with AR(1) correlation
                x[lon,lat,]<-arima.sim(list(ar = 0.3),n=31,rand.gen=rnorm,sd=0.1,mean=0)
                y[lon,lat,]<-arima.sim(list(ar = 0.3),n=31,rand.gen=rnorm,sd=0.1,mean=0)
        }
}

test3d<-modTtest3d(x,y,alternative = c("two.sided"),conf.level=0.95)
print(test3d)

walk<-walker.test(test3d$p.value, siglev=0.05)
print(walk)

fdr<-fdr.test(test3d$p.value, siglev=0.05)
print(fdr)

mbb <- mbbTest(x,y,siglev=0.05,nb=10,verbose=FALSE)
print(mbb)
