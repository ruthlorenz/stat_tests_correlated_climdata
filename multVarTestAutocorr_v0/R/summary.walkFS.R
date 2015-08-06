summary.walkFS <-
function(object, ...){
	    cat("\n")
	    msg <- paste("Results for Walker's test",sep='' )
	    print(msg)
	    cat("Field significance level: ",object$field.sig,"\n")
   	    cat("Number of significant points per time: ", object$nr.sigpt, "\n")
   	    cat("Total number of tests: ", object$total.test, "\n")
	    invisible()
}
