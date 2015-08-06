print.walkFS <-
function(x, ...) {
	    cat("Call:\n")
	    print(x$call)
	    cat("\nNumber of Significant points per time:",x$nr.sigpt,"out of ", x$total.test," tests.\n")
}
