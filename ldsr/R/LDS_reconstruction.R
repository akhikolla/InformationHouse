#' Make initial values for parameters.
#'
#' If init is a vector, make it a list of one element.
#' If init is NULL, randomize it. In this case, the function will randomize the initial value by sampling uniformly within the range for each parameters
#' (A in \[0, 1\], B in \[-1, 1\], C in \[0, 1\] and D in \[-1, 1\]).
#' @param p Dimension of u
#' @param q Dimension of v
#' @param num.restarts Number of randomized initial values
#' @return A list of initial conditions, each element is an object of class `theta`.
#' @examples
#' make_init(5, 5, 1)
#' make_init(5, 5, 2)
#' @export
make_init <- function(p, q, num.restarts) {

  replicate(num.restarts,
            {
              theta <- list(A = matrix(stats::runif(1)),
                            B = matrix(stats::runif(p, -1, 1), nrow = 1),
                            C = matrix(stats::runif(1)),
                            D = matrix(stats::runif(q, -1, 1), nrow = 1),
                            Q = matrix(1),
                            R = matrix(1),
                            mu1 = matrix(0),
                            V1 = matrix(1))
              class(theta) <- 'theta'
              theta
            },
            simplify = FALSE)
}

#' Learn LDS model with multiple initial conditions
#'
#' This is the backend computation for [LDS_reconstruction]. You should not use this directly.
#' @inheritParams LDS_EM
#' @param init A list of initial parameters for the EM search. See [make_init].
#' @param niter Maximum number of iterations, default 1000
#' @param tol Tolerance for likelihood convergence, default 1e-5. Note that the log-likelihood is normalized by dividing by the number of observations.
#' @param return.init Indicate whether the initial condition that results in the highest
#' @return a list as produced by [LDS_EM]. If return.init is true, a vector of initial condition is included in the list as well.
#' @export
LDS_EM_restart <- function(y, u, v, init, niter = 1000, tol = 1e-5, return.init = TRUE) {

  # To prevent global variable error in R CMD check
  theta0 <- NULL
  models <- foreach(theta0 = init, .packages = 'ldsr') %dopar% LDS_EM(y, u, v, theta0, niter, tol)

  # Select the model with highest likelihood
  # Only select models with C > 0 for physical interpretation (if possible)
  liks <- sapply(models, '[[', 'lik')
  allC <- lapply(lapply(models, '[[', 'theta'), '[[', 'C')
  posC <- which(allC > 0) # Positions of positive C
  if (length(posC) > 0) {
    maxInd <- which(liks == max(liks[posC], na.rm = TRUE))
  } else {
    maxInd <- which.max(liks)
  }
  ans <- models[[maxInd]]
  if (return.init) ans$init <- init[[maxInd]]

  ans
}

#' Call a reconstruction method
#'
#' Call a reconstruction method subroutine according to the method required
#' @inheritParams LDS_reconstruction
#' @param y Catchment output, preprocessed from data
#' @return The results as produced by `LDS_EM_restart()` when the default method (EM) is called. Other methods are experimental.
call_method <- function(y, u, v, method, init, num.restarts, return.init,
                        ub, lb, num.islands, pop.per.island,
                        niter, tol) {
  # Reconstruction with the supplied method
  switch(method,
         EM = LDS_EM_restart(y, u, v, init, niter, tol, return.init),
         GA = LDS_GA(y, u, v, ub, lb, num.islands, pop.per.island,
                     niter, parallel = FALSE),
         BFGS = LDS_BFGS(y, u, v, ub, lb, num.restarts, parallel = FALSE),
         BFGS_smooth = {
           fit1 <- LDS_BFGS(y, u, v, ub, lb, num.restarts, parallel = FALSE)
           fit2 <- Kalman_smoother(y, u, v, fit1$theta)
           fit1$fit <- fit2
           fit1$lik <- fit2$lik
           fit1
         })
}

#' Learn LDS model.
#'
#' The initial conditions can either be randomized (specifiled by num.restarts) or provided beforehand.
#' @inheritParams PCR_reconstruction
#' @inheritParams LDS_EM_restart
#' @param u Input matrix for a single-model reconstruction, or a list of input matrices for an ensemble reconstruction.
#' @param v Same as u.
#' @param method By default this is "EM". There are experimental methods but you should not try.
#' @param init A list, each element is a vector of initial values for the parameters. If `NULL`, will be created by `make_init()`. See [make_init] for details.
#' @param num.restarts The number of initial conditions to start the EM search; ignored if `init` is provided.
#' @param return.init If `TRUE`, the list of initial values (`init`) will be returned. This can be useful if you want to reproduce the model from this one set of initial values.
#' @param ub Upper bounds, a vector whose length is the number of parameters
#' @param lb Lower bounds
#' @param num.islands Number of islands (if method is GA; experimental)
#' @param pop.per.island Initial population per island (if method is GA; experimental)
#' @param return.raw If TRUE, state and streamflow estimates without measurement updates will be returned.
#' @return A list of the following elements
#' * rec: reconstruction results, a data.table with the following columns
#'     - year: calculated from Qa and the length of u
#'     - X: the estimated hidden state
#'     - Xl, Xu: lower and upper range for the 95\% confidence interval of X
#'     - Q: the reconstructed streamflow
#'     - Ql, Qu: lower and upper range for the 95\% confidence interval of Q
#' * theta: model parameters
#' * lik: maximum likelihood
#' * init: the initial condition that resulted in the maximum likelihood (if `return.init = TRUE`).
#' @examples
#' # Make a shorter time horizon so that this example runs fast
#' u <- v <- t(NPpc[601:813])
#' # We run this example without parallelism
#' foreach::registerDoSEQ()
#' LDS_reconstruction(NPannual, u, v, start.year = 1800, num.restarts = 1)
#' # Please refer to the vignette for the full run with parallel options. It takes a second or two.
#' @export
LDS_reconstruction <- function(Qa, u, v, start.year, method = 'EM', transform = 'log',
                               init = NULL, num.restarts = 50, return.init = FALSE,
                               ub = NULL, lb = NULL, num.islands = 4, pop.per.island = 250,
                               niter = 1000, tol = 1e-5, return.raw = FALSE) {

  # Preprocessing ------------------------------------------------------------------
  # We don't allow both u and v to be NULL for now
  single <- is.matrix(u) || is.matrix(v)
  if (single) {
    if (is.null(u)) { # v is provided
      u <- matrix(0)
      N <- ncol(v)
    } else if (is.null(v)) { # u is provided
      v <- matrix(0)
      N <- ncol(u)
    } else { # both are provided
      if (ncol(u) != ncol(v)) stop('u and v must have the same number of time steps.')
      N <- ncol(u)
    }
    if (is.null(init)) init <- make_init(nrow(u), nrow(v), num.restarts)
  } else {
    if (is.null(u)) { # v is provided
      u <- replicate(length(v), matrix(0))
      N <- ncol(v[[1]])
    } else if (is.null(v)) { # u is provided
      v <- replicate(length(u), matrix(0))
      N <- ncol(u[[1]])
    } else { # both are provided
      if (!identical(sapply(u, ncol), sapply(v, ncol)))
        stop('u and v must have the same number of time steps.')
      N <- ncol(u[[1]])
    }
    if (is.null(init))
      init <- lapply(seq_along(u), function(i) make_init(nrow(u[[i]]), nrow(v[[i]]), num.restarts))
  }

  Qa <- as.data.table(Qa)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of u is earlier than the last year of the instrumental period.')
  years <- start.year:end.year

  y <- Qa$Qa
  if (transform == 'log') {
    y <- log(y)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(y ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    y <- (y^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, and set transform = "none".')

  if (method != 'EM' && (is.null(ub) || is.null(lb)))
    stop("For GA and BFGS methods, upper and lower bounds of parameters must be provided.")

  # Attach NA and make the y matrix
  mu <- mean(y, na.rm = TRUE)
  y <- t(c(rep(NA, Qa[1, year] - start.year), # Before the instrumental period
           y - mu,     # Instrumental period
           rep(NA, end.year - Qa[.N, year]))) # After the instrumental period

  if (!(method %in% c('EM', 'GA', 'BFGS', 'BFGS_smooth')))
    stop("Method undefined. It has to be either EM, GA, BFGS or BFGS_smooth.")

  # Subroutines ------------------------------------------------------------------

  construct_rec <- function(fit, theta, mu, transform, years) {
    # Convert Y to Q and calculate confidence intervals
    X <- c(fit$X)
    V <- c(fit$V)
    Y <- c(fit$Y) + mu
    C <- c(theta$C)
    # Confidence intervals
    CI.X <- 1.96 * sqrt(V)
    varY <- C * V * C + c(theta$R)
    CI.Y <- 1.96 * sqrt(varY)

    X.out <- data.table(year = years, X = X, Xl = X - CI.X, Xu = X + CI.X)

    Q.out <- switch(transform,
                    log = exp_ci(Y, varY),
                    none = data.table(Q = Y, Ql = Y - CI.Y, Qu = Y + CI.Y),
                    boxcox =
                      if (lambda == 0) exp_ci(y, varY) else {
                        dt <- cbind(Q = Y, Ql = Y - CI.Y, Qu = Y + CI.Y)
                        dt <- as.data.table(inv_boxcox(dt, lambda))
                      })
    cbind(X.out, Q.out)
  }

  format_results <- function(results, u, v) {
    with(results, {
      rec <- construct_rec(fit, theta, mu, transform, years)
      if (ncol(u) > 1) colnames(theta$B) <- rownames(u)
      if (ncol(v) > 1) colnames(theta$D) <- rownames(v)
      ans <- list(rec = rec, theta = theta, lik = lik)
      if (return.raw) {
        raw <- propagate(theta, u, v, y)
        ans$rec2 <- construct_rec(raw, theta, mu, transform, years)
      }

      if (return.init) ans$init <- results$init
      if (method != 'EM') ans$pl <- results$pl
      if (transform == 'boxcox') ans$lambda <- lambda
      ans
    })
  }

  # Reconstruction ----------------------------------------------------------

  if (single) {
    # For single model, we can parallelize at the init level
    results <- call_method(y, u, v, method, init, num.restarts, return.init,
                           ub, lb, num.islands, pop.per.island, niter, tol)
    format_results(results, u, v)
  } else {
    # Otherwise, we can parallelize at the u level
    i <- Q <- X <- NULL
    ensemble <- foreach(i = seq_along(u), .packages = 'ldsr') %dopar% {
      results <- call_method(y, u[[i]], v[[i]], method, init[[i]], num.restarts, return.init,
                             ub, lb, num.islands, pop.per.island, niter, tol)
      format_results(results, u[[i]], v[[i]])
    }
    rec <- rbindlist(lapply(ensemble, function(x) x$rec[, list(year, X, Q)]))
    rec <- rec[, list(X = mean(X), Q = mean(Q)), by = year]
    ans <- list(rec = rec, ensemble = ensemble)

    if (return.raw) {
      rec.raw <- rbindlist(lapply(ensemble, function(x) x$rec2[, list(year, X, Q)]))
      rec.raw <- rec.raw[, list(X = mean(X), Q = mean(Q)), by = year]
      ans$rec.raw <- rec.raw
    }
    ans
  }
}

#' One cross-validation run
#'
#' Make one prediction for one cross-validation run. This is a subroutine that is called by cvLDS, without any checks. You should not need to use this directly.
#' @param z A vector of left-out points, indexed according to the intrumental period
#' @param instPeriod indices of the instrumental period in the whole record
#' @param mu Mean of the observations
#' @inheritParams cvLDS
#' @inheritParams call_method
#' @return A vector of prediction.
#' @export
one_lds_cv <- function(z, instPeriod, mu, y, u, v, method = 'EM', num.restarts = 20,
                       ub = NULL, lb = NULL, num.islands = 4, pop.per.island = 100,
                       niter = 1000, tol = 1e-6, use.raw = FALSE) {

  y[instPeriod][z] <- NA
  init <- make_init(nrow(u), nrow(v), num.restarts)
  result <- call_method(y, u, v, method, init, num.restarts, return.init = FALSE,
                        ub, lb, num.islands, pop.per.island,
                        niter, tol)
  if (use.raw) {
    raw <- propagate(result$theta, u, v, y)
    c(raw$Y[instPeriod]) + mu
  } else {
    c(result$fit$Y[instPeriod]) + mu
  }
}


#' Cross validate LDS model. This is a wrapper for [LDS_reconstruction]
#'
#' @inheritParams LDS_reconstruction
#' @inheritParams cvPCR
#' @param use.raw Whether performance metrics are calculated on the raw time series. Experimental; don't use.
#' @return A list of cross validation results
#' * metrics.dist: distribution of performance metrics across all cross-validation runs; a matrix, one column for each metric, with column names.
#' * metrics: average performance metrics; a named vector.
#' * target: the (transformed) observations used for cross-validation; a data.table with two columns (year, y)
#' * Ycv: the predicted streamflow in each cross validation run; a matrix, one column for each cross-validation run
#' * Z: the cross-validation
#' @examples
#' # Make a shorter time horizon so that this example runs fast
#' u <- v <- t(NPpc[601:813])
#' # We run this example without parallelism.
#' foreach::registerDoSEQ()
#' cvLDS(NPannual, u, v, start.year = 1800, num.restarts = 2,
#'       Z = make_Z(NPannual$Qa, nRuns = 1))
#' # Please refer to the vignette for the full run with parallel options. It takes a minute or two.
#' @export
 cvLDS <- function(Qa, u, v, start.year, method = 'EM', transform = 'log', num.restarts = 50,
                   Z = make_Z(Qa$Qa), metric.space = 'transformed', use.raw = FALSE,
                   ub = NULL, lb = NULL, num.islands = 4, pop.per.island = 100,
                   niter = 1000, tol = 1e-5) {

  # Preprocessing ------------------------------------------------------------------
  single <- is.matrix(u) || is.matrix(v)
  if (single) {
    if (is.null(u)) { # v is provided
      u <- matrix(0)
      N <- ncol(v)
    } else if (is.null(v)) { # u is provided
      v <- matrix(0)
      N <- ncol(u)
    } else { # both are provided
      if (ncol(u) != ncol(v)) stop('u and v must have the same number of time steps.')
      N <- ncol(u)
    }
  } else {
    if (is.null(u)) { # v is provided
      u <- replicate(length(v), matrix(0))
      N <- ncol(v[[1]])
    } else if (is.null(v)) { # u is provided
      v <- replicate(length(u), matrix(0))
      N <- ncol(u[[1]])
    } else { # both are provided
      if (!identical(sapply(u, ncol), sapply(v, ncol))) stop('u and v must have the same number of time steps.')
      N <- ncol(u[[1]])
    }
  }

  Qa <- as.data.table(Qa)
  end.year <- start.year + N - 1
  if (end.year < Qa[.N, year])
    stop('The last year of pc is earlier than the last year of the instrumental period.')

  if (method != 'EM' && (is.null(ub) || is.null(lb)))
    stop("For GA and BFGS methods, upper and lower bounds of parameters must be provided.")

  obs <- Qa$Qa
  if (transform == 'log') {
    obs <- log(obs)
  } else if (transform == 'boxcox') {
    bc <- MASS::boxcox(obs ~ 1, plotit = FALSE)
    lambda <- bc$x[which.max(bc$y)]
    obs <- (obs^lambda - 1) / lambda
  } else if (transform != 'none') stop('Accepted transformations are "log", "boxcox" and "none" only. If you need another transformation, please do so first, and then supplied the transformed variable in Qa, with transform = "none".')

  if (!is.list(Z)) stop("Please provide the cross-validation folds (Z) in a list.")

  years <- start.year:end.year
  instPeriod <- which(years %in% Qa$year)

  # Attach NA and make the y matrix
  mu <- mean(obs, na.rm = TRUE)
  y <- t(c(rep(NA, Qa[1, year] - start.year), # Before the instrumental period
           obs - mu,     # Instrumental period
           rep(NA, end.year - Qa[.N, year]))) # After the instrumental period

  i <- z <- NULL
  Ycv <- if (single) {
    foreach(z = Z, .packages = 'ldsr') %dopar%
      one_lds_cv(z, instPeriod, mu, y, u, v, method, num.restarts,
                 ub, lb, num.islands, pop.per.island, niter, tol, use.raw)
  } else {
    foreach(z = Z, .packages = 'ldsr') %:%
      foreach(i = seq_along(u),
              .combine = cbind, .multicombine = TRUE, .final = rowMeans) %dopar%
        one_lds_cv(z, instPeriod, mu, y, u[[i]], v[[i]], method, num.restarts,
                   ub, lb, num.islands, pop.per.island, niter, tol, use.raw)
  }

  if (metric.space == 'original') {
    if (transform == 'log') {
      Ycv <- lapply(Ycv, exp)
    } else if (transform == 'boxcox') {
      Ycv <- lapply(Ycv, inv_boxcox, lambda = lambda)
    }
    target <- Qa$Qa
  } else {
    target <- obs
  }
  # doing mapply is a lot faster than working on data.table
  metrics.dist <- mapply(calculate_metrics, sim = Ycv, z = Z, MoreArgs = list(obs = target))
  metrics.dist <- data.table(t(metrics.dist))
  metrics <- metrics.dist[, lapply(.SD, mean)]

  names(Ycv) <- seq_along(Ycv)
  setDT(Ycv, check.names = FALSE)
  Ycv[, year := Qa$year]

  list(metrics.dist = metrics.dist,
       metrics = metrics,
       target = data.table(year = Qa$year, y = target),
       Ycv = melt(Ycv, id.vars = 'year', variable.name = 'rep', value.name = 'Y'),
       Z = Z) # Retain Z so that we can plot the CV points when analyzing CV results
}
