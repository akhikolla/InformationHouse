print.baseline <- function(x, ...) {

  lower.set <- x$lower.set
  upper.set <- x$upper.set
  lambda <- x$lambda
  clambda <- x$clambda

  digit2 <- paste("%.", max(2, getOption("digits") - 5), "f", sep = "")
  digit4 <- paste("%.", max(3, getOption("digits") - 3), "f", sep = "")

  plower.set <- sprintf(digit2, lower.set)
  pupper.set <- sprintf(digit2, upper.set)
  plambda <- sprintf(digit4, lambda)
  pclambda <- sprintf(digit4, clambda)
  pset <- paste("(", plower.set, ", ", pupper.set, "]", sep = "")

  presults <- data.frame(support = pset, lambda = plambda, cum.lambda = pclambda, stringsAsFactors=FALSE)
  rownames(presults) <- paste(1:nrow(presults), ":", sep = "")

  cat("\n========================================\n")
  cat("  Baseline Hazard Estimation (NPMLE)")
  cat("\n----------------------------------------\n")
  print(presults)
  cat("========================================\n")

}
