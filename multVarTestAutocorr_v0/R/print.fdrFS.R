print.fdrFS <-
function(x, ...) {
	    cat("Call:\n")
	    print(x$call)
	    cat("\nNumber of Significant points:",x$nr.sigpt,"out of ", x$total.test," tests.\n")
}
