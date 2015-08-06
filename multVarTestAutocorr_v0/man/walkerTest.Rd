\name{walkerTest}
\alias{walkerTest}
\alias{print.walkFS}
\alias{summary.walkFS}
\title{Walkers test for field significance
}
\description{
Method to account for multiplicity and spatial correlation in the data. The Walker's test says that a global null hypotheses may be rejected at the \eqn{\alpha}{alpha} level if the smallest of the local p-values is smaller or equal than \eqn{p_{Walker}}{p_Walker}. 
}
\usage{
walkerTest(p_val, siglev = 0.05, ...)

\method{print}{walkFS}(x, ...)
\method{summary}{walkFS}(x, ...)

}
\arguments{
  \item{p_val}{ p-values as returned by \code{modTtest3d} or \code{t.test}. Need p-values for grid longitude x latitude (dimension time is optional).

}
  \item{siglev}{ field significance level alpha used for the test

}
  \item{x}{ data of class \dQuote{walkFS} as returned by \code{walkerTest}
}
  \item{\dots}{ further arguments to be passed to or from methods, not used at the moment.
}
}
\details{
The Walker's test is based on the question \dQuote{How small must the smallest local p-value be in order to reject the global null hypothese that all of the local null hypotheses are true?}
\deqn{p_{Walker}=1-(1-\alpha)^\frac{1}{K}}{p_walker=1-alpha^(1/K)}
with K=number of local null hypotheses.
}
\value{
   A list with class \code{walkFS} containing the following components:
  \item{h.value }{ numeric 0 or 1, 1=significant if minimum p-value smaller than p_walker}
  \item{p.value }{ p-values of local tests}
  \item{field.sig }{ numeric value of field significance level alpha}
  \item{nr.sigpt }{ number of singnificant tests}
  \item{total.test }{ total number of tests K}
  \item{method }{ string with name of the method}
  \item{call }{ how the function was called}

}
\references{
Katz R.W., Brown B.G., 1991, The problem of multiplicity in research on teleconnections, Int. J. Climatol., \bold{11}, p.505--513

Wilks, D.S., 2006, On \dQuote{Field Significance} and the False Discovery Rate, J. Appl. Meteorol. Climatol., \bold{45}, p. 1181--1189.
}
\author{
Ruth Lorenz
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{fdrTest}}
}
\examples{
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

walk<-walkerTest(test3d$p.value, siglev=0.05)
print(walk)

## The function is currently defined as
function (p_val, siglev = 0.05, ...) 
{
    dims_p <- dim(p_val)
    nlon <- dims_p[1]
    nlat <- dims_p[2]
    if (!is.na(dims_p[3])) {
        ntim <- dims_p[3]
    }
    else {
        ntim <- 1
        tmp <- p_val
        p_val <- array(NA, dim = c(nlon, nlat, ntim))
        p_val[, , 1] <- tmp
    }
    h_val <- array(NA, dim = c(nlon, nlat, ntim))
    for (t in 1:ntim) {
        for (lat in 1:nlat) {
            for (lon in 1:nlon) {
                if (is.na(p_val[lon, lat, t])) {
                  h_val[lon, lat, t] <- NA
                }
                else if (p_val[lon, lat, t] < siglev) {
                  h_val[lon, lat, t] <- 1
                }
                else {
                  h_val[lon, lat, t] <- 0
                }
            }
        }
    }
    K <- sum(!is.na(p_val[, , 1]))
    sig_walk <- array(NA, dim = c(nlon, nlat, ntim))
    p_min <- array(NA, dim = c(ntim))
    p_walk <- 1 - (1 - siglev)^(1/K)
    for (j in 1:ntim) {
        p_min[j] <- min(p_val[, , j], na.rm = T)
        if (p_min[j] > p_walk) {
            sig_walk[, , j] <- 0
        }
        else {
            sig_walk[, , j] <- h_val[, , j]
        }
    }
    sig_pts <- array(NA, dim = c(ntim))
    for (j in 1:ntim) {
        sig_pts[j] <- (sum(sig_walk[, , j], na.rm = T))
    }
    method <- paste("Walkers test for field significance")
    rval <- list(h.value = sig_walk, p.value = p_val, field.sig = siglev, 
        nr.sigpt = sig_pts, total.test = K, method = method, 
        call = match.call())
    class(rval) <- "walkFS"
    return(rval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
