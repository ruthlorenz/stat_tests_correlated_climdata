\name{fdrTest}
\alias{fdrTest}
\alias{print.fdrFS}
\alias{summary.fdrFS}
\title{
False discovery rate for field significance testing
}
\description{
Method to account for multiplicity and spatial correlation in the data.The False Discovery Rate (FDR) identifies locally significant tests by controlling the expected proportion of rejected local null hypotheses that are actually true. Hence, field significance is declared if at least one local null hypotheses is rejected at the level alpha (e.g. alpha=0.05)
}
\usage{
fdrTest(p_val, siglev = 0.05, ...)

\method{print}{fdrFS}(x, ...)
\method{summary}{fdrFS}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p_val}{ p-values as returned by \code{modTtest3d} or \code{t.test}. Need p-values for grid longitude x latitude.
}
  \item{siglev}{ field significance level alpha used for the test
}
  \item{x}{ data of class \dQuote{fdrFS} as returned by \code{fdrTest}
}
  \item{\dots}{ further arguments to be passed to or from methods, not used at the moment.
}
}
\details{
We reject the local test where p(i) is no greater than
\deqn{p_{FDR}&=&max_{j=1.....,k} [p(j):p(j)\leq q(\frac{j}{K})] {}\nonumber\\
&& {}	=max_{j=1.....,k} [p(j):p(j)\leq \alpha(\frac{j}{K})]}{p_FDR=max_{j=1.....,k} [p(j):p(j)<= q(j/K)]
=max_{j=1.....,k} [p(j):p(j)<= alpha(j/K)]}
with p(j) being the ith smallest p-value.
Note that the local significance level q is numerically equal to the global test level \eqn{\alpha}{alpha}.
Despite being more permissive than the Monte Carlo and Walker's test, the FDR has the advantage of not only giving a yes or no answer for a field but checking at each grid point if the p-value is small enough for field significance.
}
\value{
   A list with class \code{fdrFS} containing the following components:
  \item{h.value }{ numeric 0 or 1, 1=significant if local p-value smaller than p_fdr}
  \item{p.value }{ p-values of local tests}
  \item{field.sig }{ numeric value of field significance level alpha}
  \item{nr.sigpt }{ number of singnificant tests}
  \item{total.test }{ total number of tests K}
  \item{method }{ string with name of the method}
  \item{call }{ how the function was called}
}
\references{
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
\code{\link{walkerTest}}
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

dfr<-fdrTest(test3d$p.value, siglev=0.05)
print(fdr)

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
    fdr <- array(0, dim = c(nlon, nlat, ntim))
    sig_FDR <- array(0, dim = c(nlon, nlat, ntim))
    p <- array(NA, dim = c(nlon * nlat))
    prob_1D <- (c(p_val))
    p_sort <- sort(prob_1D, decreasing = FALSE)
    for (k in 1:K) {
        if (p_sort[k] <= (siglev * (k/K))) {
            p[k] <- p_sort[k]
        }
        else {
            p[k] <- 0
        }
    }
    p_fdr <- max(p, na.rm = T)
    fdr[which(p_val <= p_fdr)] <- 1
    sig_FDR[which(fdr == 1 & h_val == 1)] <- 1
    sig_pts <- array(NA, dim = c(ntim))
    for (j in 1:ntim) {
        sig_pts[j] <- (sum(sig_FDR[, , j], na.rm = T))
    }
    method <- paste("False Discovery Rate for field significance")
    rval <- list(h.value = sig_FDR, p.value = p_val, field.sig = siglev, 
        nr.sigpt = sig_pts, total.test = K, method = method, 
        call = match.call())
    class(rval) <- "fdrFS"
    return(rval)
  }
}
\keyword{ multivariate }