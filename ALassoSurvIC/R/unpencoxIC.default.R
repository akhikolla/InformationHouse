unpencoxIC.default <- function(lowerIC, upperIC, X, trunc = NULL, normalize.X = TRUE, cl = NULL, tol = 1e-3, niter = 1e5, string.cen = Inf, string.missing = NA, ...) {

  match.call()

  if(missing(trunc)) {
    trunc <- NULL
    ind.trunc <- FALSE
    smallest.trunc <- 0
  } else {
    ind.trunc <- TRUE
    smallest.trunc <- min(trunc)
  }

  if (!is.null(cl)) {
    if (.Platform$OS.type == "windows") {
      if (!inherits(cl, "cluster"))
        cl <- NULL
    } else {
      if (inherits(cl, "cluster")) {
        if (length(cl) < 2L)
          cl <- NULL
      } else {
        if (cl < 2)
          cl <- NULL
      }
    }
  }

  xnames <- colnames(X)

  arglist <- fun_arglist(lowerIC, upperIC, X, trunc, normalize.X, tol, niter)
  arglist$initial_lambda <- rep(1/nrow(arglist$set), nrow(arglist$set))

  message(" Now: Obtaining the unpenalized nonparametric MLE")
  unpen <- fun_unpenSurvIC(rep(0, ncol(arglist$z)), arglist)
  final.b0 <- unpen$b
  final.lambda <- unpen$lambda

  arglist$initial_lambda <- final.lambda

  message(" Now: calculating the covariance matrix")
  cov <- fun_cov_parallel(b = final.b0, theta = 0, var.h = 5, arglist, cl)

  message(" Done.")

  if (!is.null(cl)) stopCluster(cl)

  if (normalize.X == TRUE) {
    atrue_sd <- (arglist$true_sd)
    final.b <- final.b0/atrue_sd
    final.cov <- cov / (atrue_sd %*%t(atrue_sd))
  } else {
    final.b <- final.b0
    final.cov <- cov
  }

  results <- list()
  results$xnames <- xnames
  results$n <- nrow(X)
  results$b <- final.b
  results$se <- sqrt(diag(final.cov))
  results$cov <- final.cov
  results$lambda <- final.lambda
  results$lambda.set <- arglist$set
  results$convergence <- unpen$convergence
  results$iteration <- unpen$iteration
  results$ind.trunc <- ind.trunc
  results$smallest.trunc <- ifelse(ind.trunc, min(trunc), 0)
  results$normalize.X <- normalize.X

  class(results) <- "unpencoxIC"

  return(results)
}
