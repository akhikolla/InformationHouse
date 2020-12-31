print.summary.robmixglm <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  if (!inherits(x, "summary.robmixglm"))
    stop("Use only with 'summary.robmixglm' objects.\n")
  printCoefmat(x$coefficients,digits=digits,na.print="")
  cat("\n")
  print(data.frame(logLik = x$logLik,AIC = x$AIC, BIC = x$BIC,
                   row.names = " "))
  invisible(x)
 }