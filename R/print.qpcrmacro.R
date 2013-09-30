print.qpcrmacro <-
function(x, ...){
  print(x$model, correlation=FALSE)
  if (!is.null(x$df)) if (x$df > 0) cat("\nA degree of freedom of ", round(x$df,2), " is assumed for all comparisons.\n")
  print(x$pvalues)
}
