
# functions -------------------------------------------------------------------
tRoot <- function(delta, test_stat, df, conf_level) {
  #' This function is used to find the root for a t-distribution pivotal quantity 
  #'
  #' This function returns the difference between the lower tail
  #' probability of a non-central t-distribution and a confidence level q
  #' where the t-distribution has df degrees of freedom and 
  #' non-centrality parameter delta.
  #' @param delta Non-centrality parameter
  #' @param test_stat Test statistic at which to evaluate the t-distribution
  #' @param df Degrees of freedom
  #' @param conf_level Confidence level (usually alpha/2 or 1-alpha/2)
  #' @return dif Difference between t-distribution quantile and confidence level
  #' @export

  dif <- pnct(test_stat, df, delta) -  conf_level

  return(dif)
}

getDeltaCI <- function(test_stat, 
                  alpha, 
                  d,
                  n, 
                  interval) {
  #' Confidence intervals for noncentrality parameter of t-distribution
  #'
  #' This function obtains confidence intervals for the non-centrality
  #' parameter of a t-distribution.
  #' @param test_stat Test statistics
  #' @param alpha Significance level
  #' @param d Number of parameters in general linear model
  #' @param n Number of observations in initial study
  #' @param interval Interval within which to search for roots
  #' @return ep Exceedance probability with confidence intervals (vector if cutoff is scalar and matrix otherwise)
  #' @export

  delta_lower <- uniroot(f = tRoot, 
                         interval = interval, 
                         test_stat = test_stat, 
                         df = n - d,
                         conf_level = 1 - alpha / 2)$root

  delta_upper <- uniroot(f = tRoot, 
                         interval = interval, 
                         test_stat = test_stat, 
                         df = n - d,
                         conf_level = alpha / 2)$root

  delta_ci <- c(delta_lower, delta_upper)
  names(delta_ci) <- c("lower", "upper")
  return(delta_ci)
}

getDeltaCI <- Vectorize(getDeltaCI, vectorize.args = "test_stat")


exceedProb <- function(cutoff, 
                  theta_hat, 
                  sd_hat, 
                  alpha, 
                  d,
                  n, 
                  m,
                  interval = c(-100, 100),
                  lower_tail = FALSE) {
  #' Confidence intervals for the exceedance probability
  #'
  #' This function obtains confidence intervals for exceedance probability
  #' @param cutoff Cutoff values (scalar or vector)
  #' @param theta_hat Point estimate for the parameter of interest
  #' @param sd_hat Estimated standard deviation for the parameter of interest (Note: not the standard error)
  #' @param d Number of parameters in the general linear model
  #' @param alpha Significance level
  #' @param n Number of observations in the initial study
  #' @param m Number of observations in the replication study
  #' @param interval Interval within which to search for roots
  #' @param lower_tail If TRUE, reports lower tail probabilities
  #' @return ep Exceedance probability with confidence intervals
  #' @export
  #' @examples
  #' library(exceedProb)
  #' 
  #' # Sample mean -------------------------------------------------------
  #' n <- 100
  #' x <- rnorm(n = n)
  #' 
  #' theta_hat <- mean(x)
  #' sd_hat <- sd(x)
  #' 
  #' cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)
  #' 
  #' exceedProb(cutoff = cutoff, 
  #'            theta_hat = theta_hat, 
  #'            sd_hat = sd_hat, 
  #'            alpha = 0.05, 
  #'            d = 1,
  #'            n = n,
  #'            m = n)
  #' 
  #' # Linear regression -------------------------------------------------
  #' n <- 100
  #' beta <- c(1, 2)
  #' x <-runif(n = n, min = 0, max = 10)
  #' y <- rnorm(n = n, mean = cbind(1, x) %*% beta, sd = 1)
  #' 
  #' j <- 2
  #' fit <- lm(y ~ x)
  #' theta_hat <- coef(fit)[j]
  #' sd_hat <- sqrt(n * vcov(fit)[j, j])
  #' 
  #' cutoff <- seq(from = theta_hat - 0.5, to = theta_hat + 0.5, by = 0.1)
  #' 
  #' exceedProb(cutoff = cutoff, 
  #'            theta_hat = theta_hat, 
  #'            sd_hat = sd_hat, 
  #'            alpha = 0.05, 
  #'            d = length(beta),
  #'            n = n,
  #'            m = n)

  if (length(cutoff) < 1 | !is.numeric(cutoff)) {
    stop("cutoff must be numeric a numeric vector with at least 1 element")
  }

  if (length(theta_hat) != 1 | !is.numeric(theta_hat)) {
    stop("theta_hat must be numeric scalar")
  }

  if (length(sd_hat) != 1 | !is.numeric(sd_hat)) {
    stop("sd_hat must be numeric scalar")
  }

  if (length(d) != 1 | !is.numeric(d)) {
    stop("d must be numeric scalar")
  }

  if (length(n) != 1 | !is.numeric(n)) {
    stop("n must be numeric scalar")
  }

  if (length(m) != 1 | !is.numeric(m)) {
    stop("m must be numeric scalar")
  }

  if (length(interval) != 2 | !is.numeric(interval)) {
    stop("interval must be a numeric vector of length 2")
  }

  if (diff(interval) <= 0) {
    stop("interval must have the form (a, b) for a < b")
  }

  if(!is.logical(lower_tail) | length(lower_tail) != 1) {
    stop("lower_tail must be a logical scalar")
  }

  test_stat <- sqrt(n) * (cutoff - theta_hat) / sd_hat

  delta_ci <- getDeltaCI(test_stat = test_stat,
                         d = d,
                         alpha = alpha, 
                         n = n,
                         interval = interval)

  if (lower_tail) {

    lower <- stats::pnorm(sqrt(m/n) * delta_ci["lower", ], lower.tail = TRUE)
    upper <- stats::pnorm(sqrt(m/n) * delta_ci["upper", ], lower.tail = TRUE)
    point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = TRUE)
  } else {

    lower <- stats::pnorm(sqrt(m/n) * delta_ci["upper", ], lower.tail = FALSE)
    upper <- stats::pnorm(sqrt(m/n) * delta_ci["lower", ], lower.tail = FALSE)
    point <- stats::pnorm(q = sqrt(m) * (cutoff - theta_hat) / sd_hat, lower.tail = FALSE)
  }

  ep <- data.frame(cutoff = cutoff, point = point, lower = lower, upper = upper)
  rownames(ep) <- NULL

  return(ep)
}
