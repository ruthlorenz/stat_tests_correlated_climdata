summary.mbbtest <-
function(object, ...){
   cat("\n")
   msg <- paste("Results for moving blocks bootstrapping ", object$data.name[2], " compared against ", object$data.name[1], sep="")
   print(msg)
   cat("\n")
   cat("Local significance level: ", object$loc.sig, "\n")
   cat("Field significance level: ", object$field.sig, "\n")
   cat("Number of locally significant tests: ",object$nr.sigpt, "\n")
   cat("Total number of tests: ", object$total.test, "\n")
   cat("Local significance coverage: ", object$loc.sig.coverage, "\n")
   cat("Result of field significance test: ", object$FS, "\n")
   invisible()
}
