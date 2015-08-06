\name{modTtest3d}
\alias{modTtest3d}
\alias{modTtest3d.default}
\alias{modTtest3d.formula}
\title{
Modified t-test for [lon,lat,time] data
}
\description{
Modified student t-test for autocorrelated data
using equivalent sample size. The method is based on
Zwiers and von Storch, J. Clim 1995, p. 336-351.
For small sample sizes (\eqn{\leq 30}{<=30}) a lookup table test is
used (also described in above publication). This function only works for 3d [lon,lat,time] data.
}
\usage{
modTtest3d(x, ...)

\method{modTtest3d}{default}(x, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95, timevals=NULL, \dots)

\method{modTtest3d}{formula}(formula, data, subset, na.action, \dots)
}

}

\arguments{
  \item{x}{a (non-empty) numeric array of data values with dimensions longitude, latitude, time where time is autocorrelated AR(1)
}
  \item{y}{an optional (non-empty) numeric array of data values with dimensions longitude, latitude, time where time is autocorrelated AR(1)
}
  \item{alternative}{a character string specifying the alternative
    hypothesis, must be one of \code{"two.sided"} (default),
    \code{"greater"} or \code{"less"}.  You can specify just the initial
    letter.
}
  \item{mu}{a number indicating the true value of the mean (or difference in means if you are performing a two sample test).
}
  \item{paired}{a logical indicating whether you want a paired test.
}
  \item{var.equal}{a logical variable indicating whether to treat the two variances as being equal. If \code{TRUE} then the pooled variance is used to estimate the variance, if \code{FALSE} a warning is issued because the method assumes equal variances.
}
  \item{conf.level}{confidence level of the interval.
}
  \item{timevals}{an optional number of time values in the x (y) timeseries. If provided the function checks for correct dimensions.
}
  \item{formula}{a formula of the form \code{lhs ~ rhs} where \code{lhs} is a numeric variable giving the data values and \code{rhs} a factor with two levels giving the corresponding groups.
}
  \item{data}{an optional matrix or data frame (or similar: see \code{\link{model.frame}}) containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
}
  \item{subset}{an optional vector specifying a subset of observations to be used.
}
  \item{na.action}{a function which indicates what should happen when the data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
}
  \item{\dots}{further arguments to be passed to or from methods, not used at the moment.

}
}
\details{
Performs a modified t-test using the equivalent sample size in case of autocorrelation in the data (Zwiers and von Storch, 1995).  When the sample size is too small (\eqn{\leq 30}{<= 30}) to obtain robust results a look-up table test is performed.

Either one sample or two sample tests can be performed, the function needs timeseries of at least 12, otherwise the lookup table test fails. The lookup table test also is not suited for negative correlations \eqn{\leq -0.3}{<= -0.3}. For negative autocorrelations and sample sizes larger than 30 a standard student t-test is performed.

The formula interface is only applicable for the 2-sample tests.

alternative = "greater" is the alternative that \code{x} has a larger mean than \code{y}.

If paired is \code{TRUE} then both x and y must be specified and they must be the same length. Missing values are silently removed (in pairs if paired is \code{TRUE}).

If the input data are effectively constant (compared to the larger of the two means) an error is generated.
}
\value{
A list with class \code{htest} containing the following components:

\item{statistic }{the value of the t-statistic.}
\item{parameter }{the degrees of freedom for the t-statistic. Is NA for the lookup table test}
\item{p.value }{the p-value for the test. Is NA for the lookup table test}
\item{conf.int }{a confidence interval for the mean appropriate to the specified alternative hypothesis.}
\item{estimate }{the estimated mean or difference in means depending on whether it was a one-sample test or a two-sample test.}
\item{null.value }{the specified hypothesized value of the mean or mean difference depending on whether it was a one-sample test or a two-sample test.}
\item{alternative }{a character string describing the alternative hypothesis.}
\item{method }{a character string indicating what type of t-test was performed.}
\item{data.name }{a character string giving the name(s) of the data.}

}
\references{
Zwiers and von Storch, 1995, Taking serial correlation into account in tests of the mean, J. Clim, \bold{8}, p. 336--351.
}
\author{
Ruth Lorenz
}

\seealso{
 \code{\link{t.test}},  \code{\link[modTtest]{modTtest}}
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

## The function is currently defined as
function (x, ...) 
UseMethod("modTtest3d")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ univar }
\keyword{ ts }
