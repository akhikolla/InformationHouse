#' Transform the estimates before calculating metrics
#'
#' If you already ran the cross-validation on transformed output and now wanted to calculate performance on the back-transformed one, or vice-versa, you don't have to rerun the whole cross-validation, but just need to transform or back-transform the cross-validation Ycv. This function helps you do that.
#' @param cv Cross-validation output as produced by cvLDS or cvPCR
#' @param transform Either "log", "exp", "boxcox" or "inv_boxcox"
#' @param lambda Lambda value used in Box-Cox or inverse Box-Cox
#' @return A new cv object wit hthe new metrics
#' @examples
#' # Cross-validate with log-transform
#' cv <- cvPCR(NPannual, NPpc, start.year = 1200, transform = 'log', metric.space = 'transformed')
#' # Calculate metrics based on back-transformed values
#' m <- metrics_with_transform(cv, 'exp')
#' @export
metrics_with_transform <- function(cv, transform, lambda = NULL) {
  Ycv <- copy(cv$Ycv)
  Y <- y <- NULL
  Ycv[,
      Y := switch(transform,
                 log = log(Y),
                 exp = exp(Y),
                 boxcox = {if (lambda == 0) log(Y) else (Y^lambda - 1) / lambda},
                 inv_boxcox = {if (lambda == 0) exp(Y) else (Y*lambda + 1)^(1/lambda)})
      ]
  obs <- copy(cv$target)
  obs[,
      y := switch(transform,
                 log = log(y),
                 exp = exp(y),
                 boxcox = {if (lambda == 0) log(y) else (y^lambda - 1) / lambda},
                 inv_boxcox = {if (lambda == 0) exp(Y) else (y*lambda + 1)^(1/lambda)})
      ]

  metrics.dist <- Ycv[, data.table(t(calculate_metrics(Y, obs$y, cv$Z[[rep]]))), by = rep]
  metrics.dist[, rep := NULL]
  metrics <- metrics.dist[, lapply(.SD, mean)]

  list(metrics.dist = metrics.dist,
       metrics = metrics,
       target = obs,
       Ycv = Ycv,
       Z = cv$Z)
}

#' Reconstruction metrics
#'
#' Calculate reconstruction metrics from the instrumental period
#' @param sim A vector of reconstruction output for instrumental period
#' @param obs A vector of all observations
#' @param z A vector of left out indices in cross validation
#' @param norm.fun The function (unquoted name) used to calculate the normalizing constant. Default is `mean()`, but other functions such as `sd()` can also be used. THe function must take a vector as input and return a scalar as output, and must have an argument `na.rm = TRUE`.
#' @return A named vector of performance metrics
#' @examples
#' calculate_metrics(rnorm(100), rnorm(100), z = 1:10)
#' calculate_metrics(rnorm(100), rnorm(100), z = 1:10, norm.fun = sd)
#' @export
calculate_metrics <- function(sim, obs, z, norm.fun = mean) {
    train.obs <- obs[-z]
    obsInd <- which(!is.na(train.obs))
    train.obs <- train.obs[obsInd]
    train.sim <- sim[-z][obsInd]
    val.sim <- sim[z]
    val.obs <- obs[z]

    c(R2    = NSE(train.sim, train.obs), # Use the NSE form of R2
      RE    = RE(val.sim, val.obs, mean(train.obs)),
      CE    = NSE(val.sim, val.obs),
      nRMSE = nRMSE(val.sim, val.obs, norm.fun(obs, na.rm = TRUE)),
      KGE   = KGE(val.sim, val.obs)
    )
}

#' Make cross-validation folds.
#'
#' Make a list of cross-validation folds. Each element of the list is a vector of the cross-validation points for one cross-validation run.
#' @param obs Vector of observations.
#' @param nRuns Number of repetitions.
#' @param frac Fraction of left-out points. For leave-one-out, use `frac = 1`, otherwise use any value less than 1. Default is 0.1 (leave-10%-out).
#' @param contiguous Logical. If `TRUE`, the default, the left-out points are made in contiguous blocks; otherwise, they are scattered randomly.
#' @return A list of cross-validation folds
#' @examples
#' Z <- make_Z(NPannual, nRuns = 30, frac = 0.25, contiguous = TRUE)
#' @export
make_Z <- function(obs, nRuns = 30, frac = 0.1, contiguous = TRUE) {
  obsInd <- which(!is.na(obs))
  if (frac == 1) {
    split(obsInd, obsInd)
  } else {
    n <- length(obsInd)
    k <- floor(n * frac) # leave-k-out
    if (contiguous) {
      maxInd <- n - k # Highest possible position in of obsInd
      if (maxInd < nRuns) { # Not enough samples, reduce k
        maxInd <- nRuns
        k <- n - nRuns
      }
      lapply(sort(sample(1:maxInd, nRuns)), function(x) obsInd[x:(x + k)])
    } else {
      replicate(nRuns, sort(sample(obsInd, k, replace = FALSE)), simplify = FALSE)
    }
  }
}

#' Inverse Box-Cox transform
#'
#' @param x A numeric vector to be transformed
#' @param lambda Lambda parameter
#' @return The inverse Box-Cox
#' @examples
#' inv_boxcox(x = rnorm(10), lambda = 1)
#' inv_boxcox(x = rnorm(10), lambda = 0) # exp(x)
#' @export
inv_boxcox <- function(x, lambda) {
  if (lambda == 0) exp(x) else (x*lambda + 1)^(1/lambda)
}

#' Exponential confidence interval
#'
#' Get the confidence interval of Q = exp(Y) when Y is normal, i.e, Q is log-normal.
#' @param y A vector of model estimates
#' @param sigma2 The variance of y as learned from a model
#' @return A data.table with the updated confidence intervals
exp_ci <- function(y, sigma2) {
  dist <- 1.96 * sqrt((exp(sigma2) - 1) * exp(2 * y + sigma2))
  Q <- exp(y)
  data.table(Q = Q, Ql = Q - dist, Qu = Q + dist)
}
