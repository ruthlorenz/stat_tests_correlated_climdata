\name{modTtest3d.default}
\alias{modTtest3d.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
modTtest3d.default(x, y = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95, time.vals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{alternative}{
%%     ~~Describe \code{alternative} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{paired}{
%%     ~~Describe \code{paired} here~~
}
  \item{var.equal}{
%%     ~~Describe \code{var.equal} here~~
}
  \item{conf.level}{
%%     ~~Describe \code{conf.level} here~~
}
  \item{time.vals}{
%%     ~~Describe \code{time.vals} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x, y = NULL, alternative = c("two.sided", "less", "greater"), 
    mu = 0, paired = FALSE, var.equal = TRUE, conf.level = 0.95, 
    time.vals = NULL, ...) 
{
    library(abind)
    source("calc_lag1Corr_3d.R")
    source("calc_EqvSamplSize_large_3d.R")
    source("tcrit_for_modTtest_small_3d.R")
    source("calc_poolSampVar_3d.R")
    alternative <- match.arg(alternative)
    if (!missing(mu) && (length(mu) != 1 || is.na(mu))) 
        stop("'mu' must be a single number")
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) || 
        conf.level < 0 || conf.level > 1)) 
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        if (paired) 
            xok <- yok <- complete.cases(x[1, 1, ], y[1, 1, ])
        else {
            yok <- !is.na(y[1, 1, ])
            xok <- !is.na(x[1, 1, ])
        }
        y <- y[, , yok, drop = FALSE]
    }
    else {
        dname <- deparse(substitute(x))
        if (paired) 
            stop("'y' is missing for paired test")
        xok <- !is.na(x[1, 1, ])
        yok <- NULL
    }
    x <- x[, , xok, drop = FALSE]
    if (paired) {
        x <- x - y
        y <- NULL
    }
    dimsX <- dim(x)
    nx <- dimsX[3]
    if (!is.null(time.vals)) {
        if (nx != time.vals) {
            stop("Error: wrong dimension, the time in X does not correspond to time.vals")
        }
    }
    if (!is.null(y)) {
        dimsY <- dim(y)
        ny <- dimsY[3]
        if (dimsX[2] != dimsY[2]) {
            stop("Error: X and Y have different number of latitudes")
        }
        else if (dimsX[1] != dimsY[1]) {
            stop("Error: X and Y have different number of longitudes")
        }
    }
    nlat <- dimsX[2]
    nlon <- dimsX[1]
    mx <- array(NA, dim = c(nlon, nlat, 1))
    vx <- array(NA, dim = c(nlon, nlat, 1))
    mx[, , 1] <- apply(x, c(1:2), mean, na.rm = TRUE)
    vx[, , 1] <- apply(x, c(1:2), var, na.rm = TRUE)
    if (is.null(y)) {
        if (nx < 2) 
            stop("not enough 'x' observations")
        r1 <- calc_lag1Corr(x, mx)
        if (nx >= 30) {
            ne <- calc_EqvSamplSize_large(x, r1)
            df <- ne - 1
            stderr <- sqrt(vx/ne)
            if (all(stderr < 10 * .Machine$double.eps * abs(mx))) 
                stop("data are essentially constant")
            tstat <- (mx - mu)/stderr
            method <- "One Sample modified t-test for spatial data"
            estimate <- setNames(mx, "mean of x")
        }
        else {
            df <- NA
            stderr <- sqrt(vx/nx)
            tstat <- (mx - mu)/stderr
            tcrit <- lookupTcrit(r1, conf.level - 1, nx)
            method <- if (paired) 
                "Paired modified t-test for spatial data"
            else "One Sample modified lookup table t-test for spatial data"
            estimate <- setNames(mx, if (paired) 
                "mean of the differences"
            else "mean of x")
        }
    }
    else {
        dimsY <- dim(y)
        ny <- dimsY[3]
        my <- array(NA, dim = c(nlon, nlat, 1))
        vy <- array(NA, dim = c(nlon, nlat, 1))
        my[, , 1] <- apply(y, c(1:2), mean, na.rm = TRUE)
        vy[, , 1] <- apply(y, c(1:2), var, na.rm = TRUE)
        if (nx < 2) 
            stop("not enough 'x' observations")
        if (ny < 2) 
            stop("not enough 'y' observations")
        if (!var.equal) 
            stop("the two sample modified t-test assumes equal variance")
        if (nx + ny < 3) 
            stop("not enough observations")
        sp2 <- calc_poolSampVar(x, y, mx, my)
        estimate <- my - mx
        if (nx + ny >= 30) {
            r1 <- calc_lag1Corr(x, y, mx, my)
            EqvSize <- calc_EqvSamplSize_large(x, y, r1)
            ne <- EqvSize$NE
            xne_p <- EqvSize$Xne_p
            yne_p <- EqvSize$Yne_p
            df <- xne_p + yne_p - 2
            stderr <- (sqrt(sp2) * ne)
            tstat <- (my - mx)/stderr
            method <- paste("Two Sample modified t-test for spatial data")
        }
        else {
            df <- NA
            stderr <- (sqrt(sp2) * (sqrt(1/nx + 1/ny)))
            tstat <- (my - mx)/stderr
            r1 <- calc_lag1Corr(x, y, mx, my)
            tcrit <- lookupTcrit(r1, 1 - conf.level, nx + ny)
            method <- paste("Two Sample modified lookup table t-test for spatial data")
        }
        if (all(stderr < 10 * .Machine$double.eps * max(abs(mx), 
            abs(my)))) 
            stop("data are essentially constant")
    }
    cint <- array(NA, dim = c(2, nlon, nlat, 1))
    if (all(!is.na(df))) {
        if (alternative == "less") {
            pval <- pt(tstat, df)
            cint[1, , , ] <- -Inf
            cint[2, , , ] <- tstat + qt(conf.level, df)
        }
        else if (alternative == "greater") {
            pval <- pt(tstat, df, lower.tail = FALSE)
            cint[1, , , ] <- tstat - qt(conf.level, df)
            cint[2, , , ] <- Inf
        }
        else {
            pval <- 2 * pt(-abs(tstat), df)
            alpha <- 1 - conf.level
            ci <- qt(1 - alpha/2, df)
            cint[1, , , ] <- tstat - ci
            cint[2, , , ] <- tstat + ci
        }
    }
    else {
        pval <- NA
        if (alternative == "less") {
            cint[1, , , ] <- -Inf
            cint[2, , , ] <- tcrit
        }
        else if (alternative == "greater") {
            cint[1, , , ] <- tcrit
            cint[2, , , ] <- Inf
        }
        else {
            ci <- tcrit
            cint[1, , , ] <- tstat - ci
            cint[2, , , ] <- tstat + ci
        }
    }
    cint[1, , , ] <- mu + adrop(cint[1, , , , drop = FALSE], 
        drop = 1) * stderr
    cint[2, , , ] <- mu + adrop(cint[2, , , , drop = FALSE], 
        drop = 1) * stderr
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y)) 
        "difference in means"
    else "mean"
    attr(cint, "conf.level") <- conf.level
    rval <- list(statistic = drop(tstat), parameter = drop(df), 
        p.value = drop(pval), conf.int = drop(cint), estimate = drop(estimate), 
        null.value = mu, alternative = alternative, method = method, 
        data.name = dname)
    class(rval) <- "htest"
    return(rval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
