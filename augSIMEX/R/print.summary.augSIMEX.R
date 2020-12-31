print.summary.augSIMEX <-
function (x, digits = max(3, getOption("digits") - 3), 
          signif.stars = getOption("show.signif.stars"), ...) 
{   
    cat("\nCall:\n",
      paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("\nThe error-prone variable: ")
    cat(x$err.var, sep = ", ")
    cat("\nThe corresponding true variable: ")
    cat(x$err.true, sep = ", ")
    cat("\nNumber of simulations: ", x$B, sep = "")
    cat("\nNumber of iterations in bootstrap: ", x$nBoot, "\n\n", sep = "")
    
    if(length(x$aliased) == 0L) {
      cat("\nNo Coefficients\n")
    } else {
      cat("\nCoefficients:\n")
      coefs <- x$coefficients
      if(!is.null(aliased <- x$aliased) && any(aliased)) {
        cn <- names(aliased)
        coefs <- matrix(NA, length(aliased), 4L,
                        dimnames=list(cn, colnames(coefs)))
        coefs[!aliased, ] <- x$coefficients
      }
      printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
                   na.print="NA", ...)
    }
    cat("\n")
    
    cat("\n(Dispersion parameter for ", x$family$family,
        " family taken to be ", format(x$dispersion), ")\n\n",
        apply(cbind(paste(format(c("Null","Residual"), justify="right"),
                          "deviance:"),
                    format(unlist(x[c("null.deviance","deviance")]),
                           digits = max(5L, digits + 1L)), " on",
                    format(unlist(x[c("df.null","df.residual")])),
                    " degrees of freedom\n"),
              1L, paste, collapse = " "), sep = "")
    if(nzchar(mess <- naprint(x$na.action))) cat("  (", mess, ")\n", sep = "")
    cat("AIC: ", format(x$aic, digits = max(4L, digits + 1L)),"\n\n",sep = "")
    ##
    cat("\n")
    invisible(x)

}
