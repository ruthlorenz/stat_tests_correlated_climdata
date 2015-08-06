\name{mbbTest.default}
\alias{mbbTest.default}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mbbTest.default(x, y, L = 1, siglev = 0.05, alpha = 0.05, nb = 1000, verbose = FALSE, timevals = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{L}{
%%     ~~Describe \code{L} here~~
}
  \item{siglev}{
%%     ~~Describe \code{siglev} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{nb}{
%%     ~~Describe \code{nb} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
}
  \item{timevals}{
%%     ~~Describe \code{timevals} here~~
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
function (x, y, L = 1, siglev = 0.05, alpha = 0.05, nb = 1000, 
    verbose = FALSE, timevals = NULL, ...) 
{
    source("/home/z3441306/scripts/plot_scripts/R_scripts/tests_for_autocorr_fields/estimate_block_length.R")
    source("/home/z3441306/scripts/plot_scripts/R_scripts/tests_for_autocorr_fields/move_blocks_bootstrap_wilks.R")
    dimsX <- dim(x)
    dimsY <- dim(y)
    if (all(dimsX != dimsY)) {
        stop("wrong dimensions, dim of X must be the same as dim of Y")
        if (dimsX[3] != dimsY[3]) {
            stop("Error: X and Y have different number of time")
        }
        else if (dimsX[2] != dimsY[2]) {
            stop("Error: X and Y have different number of latitudes")
        }
        else if (dimsX[1] != dimsY[1]) {
            stop("Error: X and Y have different number of longitudes")
        }
        else {
            nyears <- dimsX[3]
            nlat <- dimsX[2]
            nlon <- dimsX[1]
        }
    }
    if (!is.null(timevals)) {
        if (nyears != timevals) {
            stop("Error: wrong dimension, the time in X and Y does not correspond to timevals")
        }
    }
    out <- list()
    data_name <- c(as.character(substitute(x)), as.character(substitute(y)))
    names(data_name) <- c("Control", "Experiment")
    out$data.name <- data_name
    x_avg <- apply(x, c(1:2), mean)
    y_avg <- apply(y, c(1:2), mean)
    anom_seas <- x_avg - y_avg
    est_L <- estBlockL(x, y)
    L <- est_L$L
    V_prime_exp <- est_L$V_prime_x
    V_prime_ctl <- est_L$V_prime_y
    exp_var <- apply(x, c(1:2), var)
    ctl_var <- apply(y, c(1:2), var)
    sigma_s <- (V_prime_exp * exp_var/nyears + V_prime_ctl * 
        ctl_var/nyears)^(1/2)
    d <- anom_seas/sigma_s
    D <- max(abs(d))
    move_boot <- moveBlocksBoot(x, y, L, nb, verbose)
    d_star <- move_boot$d_star
    D_star <- move_boot$D_star
    for (lon in 1:nlon) {
        for (lat in 1:nlat) {
            I_u <- sum(d_star[lon, lat, 1:nb] >= d[lon, lat], 
                na.rm = T)
            I_l <- sum(d_star[lon, lat, 1:nb] <= d[lon, lat], 
                na.rm = T)
            if (I_u/(nb + 1) <= (as.numeric(siglev)/2)) {
                signif[lon, lat] <- 1
            }
            else if (I_l/(nb + 1) <= (as.numeric(siglev)/2)) {
                signif[lon, lat] <- 1
            }
            else {
                signif[lon, lat] <- 0
            }
        }
    }
    loc_sig_coverage <- sum(signif, na.rm = T)/sum(!is.na(x[, 
        , 1]))
    I_u = sum(D_star[1:nb] >= D)
    I_l = sum(D_star[1:nb] <= D)
    if (I_u/(nb + 1) <= (as.numeric(siglev)/2)) {
        fs = 1
    }
    else if (I_l/(nb + 1) <= (as.numeric(siglev)/2)) {
        fs = 1
    }
    else {
        fs = 0
    }
    method <- paste("Moving blocks bootstrapping for autocorrelated fields")
    out$estimate <- anom_seas
    out$loc.sig <- siglev
    out$field.sig <- alpha
    out$h.value <- signif
    out$h.valueFS <- fs
    out$nr.sigpt <- sum(signif, na.rm = T)
    out$total.test <- sum(!is.na(x[, , 1]))
    out$loc.sig.coverage <- loc_sig_coverage
    out$FS <- if (fs == 1) 
        TRUE
    else FALSE
    out$nboot <- nb
    out$call <- match.call()
    class(out) <- "mbbtest"
    return(out)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
