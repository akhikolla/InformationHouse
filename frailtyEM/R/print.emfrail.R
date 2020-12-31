#' @export
#' @method print emfrail
#' @keywords internal
print.emfrail <- function(x, ...) {

  cat("Call: \n")
  dput(attr(x, "call"))
  cat("\n")

  cat("log-likelihood:", x$loglik[2], "\n")

  with(x$distribution,
       if(dist == "gamma" | (dist == "pvf" & pvfm < 0))
         cat("frailty variance:",
             exp(-x$logtheta),
             "\n"))

  cat("theta:", exp(x$logtheta), "\n")
  cat("\n")
  if(length(x$coef) > 0) {
    coefmat <- list(
      coef = x$coef,
      "exp(coef)" = exp(x$coef),
      "se(coef)" = sqrt(diag(x$var)[seq_along(x$coef)]),
      "adj. se" = sqrt(diag(x$var_adj)[seq_along(x$coef)] ))

    if(all(is.na(coefmat$`adj. se`))) {
      coefmat$`adj. se` <- NULL
      coefmat$z <- coefmat$coef / coefmat$`se(coef)`
    } else {
      coefmat$z <- coefmat$coef / coefmat$`adj. se`
    }

    coefmat$p <-  1 - pchisq(coefmat$z^2, df = 1)

    coefmat <- do.call(cbind, coefmat)

    printCoefmat(coefmat)
  }

  if(!is.null(x$ca_test)) {
    cat("\n")
    cat("Score test for heterogeneity: p-val", format(x$ca_test[3], digits = 3))
  }

  if(!is.null(x$cens_test)) {
    cat("\n")
    cat("Score test for dependent censoring: p-val", format(x$cens_test[2], digits = 3))
  }

  invisible(x)

}




