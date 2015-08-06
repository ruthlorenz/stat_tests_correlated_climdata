summary.fdrFS <-
function(object, ...){
	    cat("\n")
	    msg <- paste("Results for FDR-test",sep='' )
	    print(msg)
	    cat("Field significance level: ",object$field.sig,"\n")
   	    cat("Number of significant points: ", object$nr.sigpt, "\n")
	    cat("Total number of tests: ", object$total.test, "\n")
	    invisible()
}
