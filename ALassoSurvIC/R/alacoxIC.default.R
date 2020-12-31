alacoxIC.default <- function(lowerIC, upperIC, X, trunc, theta, normalize.X = TRUE,  cl = NULL, max.theta = 1e3, tol = 1e-3, niter = 1e5, string.cen = Inf, string.missing = NA, ...) {

  match.call()

  if(missing(trunc)) {
    trunc <- NULL
    ind.trunc <- FALSE
    smallest.trunc <- 0
  } else {
    ind.trunc <- TRUE
    smallest.trunc <- min(trunc)
  }

  if(missing(theta)) {
    theta <- NULL
    user.theta <- FALSE
  } else {
    user.theta <- TRUE
  }

  if(!any(is.numeric(theta), is.null(theta))) stop("The input for 'theta' is not numeric.")

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

  message(" Now: obtaining the unpenalized NPMLE")
  initial <- fun_unpenSurvIC(rep(0, ncol(arglist$z)), arglist)

  tilde_b <- initial$b
  arglist$initial_lambda <- initial$lambda

  if (is.null(theta)) {
    bic_b_cvg <- fun_est_parallel(max.theta, tilde_b, arglist, cl)
  } else {
    message(" Now: estimating beta with the user input theta")
    bic_b_cvg <- fun_penSurvIC(theta = theta, tilde_b, arglist)
    can_b <- bic_b_cvg$b
    log_pen <- log_penlikelihood(can_b, arglist)
    n <- arglist$n
    bic_b_cvg$bic <- -2 * log_pen + log(n) * sum(can_b != 0)
  }

  final.b.BIC <- bic_b_cvg$b
  final.theta <- bic_b_cvg$theta
  final.bic <- bic_b_cvg$bic
  final.lambda <- bic_b_cvg$lambda

  message(" Now: calculating the covariance matrix")
  cov <- fun_cov_parallel(b = final.b.BIC, theta = final.theta, var.h = 5, arglist, cl)

  message(" Done.")

  if (!is.null(cl)) stopCluster(cl)

  if (normalize.X == TRUE) {
    atrue_sd <- (arglist$true_sd)
    final.b <- final.b.BIC/atrue_sd
    final.cov <- cov / (atrue_sd %*%t(atrue_sd))
  } else {
    final.b <- final.b.BIC
    final.cov <- cov
  }

  results <- list()
  results$xnames <- xnames
  results$n <- nrow(X)
  results$b <- final.b
  results$se <- sqrt(diag(final.cov))
  results$cov <- final.cov
  results$theta <- final.theta
  results$user.theta <- user.theta
  results$bic <- final.bic
  results$lambda <- final.lambda
  results$lambda.set <- arglist$set
  results$unpen.b <- tilde_b
  results$convergence <- bic_b_cvg$convergence
  results$iteration <- bic_b_cvg$iteration
  results$ind.trunc <- ind.trunc
  results$smallest.trunc <- smallest.trunc
  results$normalize.X <- normalize.X

  class(results) <- "alacoxIC"

  return(results)

}
