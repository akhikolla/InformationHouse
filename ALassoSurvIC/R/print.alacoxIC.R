print.alacoxIC <- function(x, ...) {

  xnames <- x$xnames
  n <- x$n
  ind.normalize <- x$normalize.X
  ind.convergence <- x$convergence
  ind.trun <- x$ind.trun
  user.theta <- x$user.theta
  pcon <- ifelse(ind.convergence, "converged", "failed to converge")
  iteration <- x$iteration
  coef <- x$b
  theta <- x$theta
  bic <- x$bic
  ind_zero <- coef == 0
  nb <- length(coef)
  se <- x$se
  z <- coef/se
  pvalue <- pnorm(abs(z), lower.tail = FALSE)
  ci <- cbind(coef - qnorm(0.975) * se, coef + qnorm(0.975) * se)

  digit <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")
  pcoef <- sprintf(digit, coef)
  pz <- sprintf(digit, z)
  ppvalue <- fun_less(sprintf(digit, pvalue))
  pci <- t(apply(ci, 1, function(x, digit) sprintf(digit, x), digit = digit))
  pse <- sprintf(digit, se)
  pse[ind_zero] <- "-"
  pz[ind_zero] <- "-"
  ppvalue[ind_zero] <- "-"
  pci[ind_zero, ] <- "-"

  if(is.null(xnames)) {
    pxnames <- paste("  X", 1:nb, sep = "")
  } else {
    pxnames <- paste("  ", xnames, sep = "")
  }

  results <- cbind(pcoef, pci, pse, pz, ppvalue)
  colnames(results) <- c("coef", "lower.CI", "upper.CI", "se", "z", "p")
  rownames(results) <- pxnames
  presults <- as.data.frame(results)
  presults
  cat("\n=======================================================\n")
  cat("    Variable Selection - adpative lasso estimation\n")
  if(ind.trun) {
    cat("  (Input: interval censored and left truncated data)")
  } else {
    cat("          (Input: interval censored data)")
  }
  cat("\n-------------------------------------------------------\n")
  print(presults)
  cat("-------------------------------------------------------\n")
  cat(paste(" * n = ", n, "\n", sep = ""))
  if(user.theta) {
    cat(paste(" * The user input tuning parameter = ", round(theta, 2), "\n", sep = ""))
  } else {
    cat(paste(" * Tuning parameter selected by BIC = ", round(theta, 2),"\n", sep = ""))
  }
  cat(paste(" * EM algorithm ", pcon, " after ", iteration, " iterations\n", sep = ""))
  if(isTRUE(x$normalize.X)) cat(paste(" * X is normalized", sep = ""))
  cat("\n=======================================================\n")
}
