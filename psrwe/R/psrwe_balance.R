#' Distance between two distributions
#'
#' Get balance by different metrics
#'
#' @param cov0 Vector (or matrix for \code{mhb}) of samples from the first
#'     distribution
#' @param cov1 Vector (or matrix for \code{mhb}) of samples from the second
#'     distribution
#' @param metric Metrics of distances with options
#' \describe{ \item{ovl}{Overlapping area}
#'     \item{ksd}{Kullback-Leibler distance} \item{std}{Standardized difference
#'     in mean} \item{abb}{Absolute difference in mean} \item{ley}{Levy
#'     distance} \item{mhb}{Mahalanobis distance}
#' }
#'
#' @return A real value of the distance
#'
#' @examples
#'
#' x <- rnorm(100,  mean = 0, sd = 1)
#' y <- rnorm(1000, mean = 1, sd = 2)
#' get_distance(x, y, "ovl")
#' get_distance(x, y, "abd")
#'
#' @export
#'
get_distance <- function(cov0, cov1,
                         metric = c("ovl", "ksd", "std", "abd",
                                    "ley", "mhb")) {
    metric <- match.arg(metric)
    switch(metric,
           std = {
               s <- sqrt((var(cov1) + var(cov0)) / 2)
               abs(mean(cov1) - mean(cov0)) / s
           },
           abd = abs(mean(cov0) - mean(cov1)),
           ovl = metric_ovl(cov0, cov1),
           ksd = metric_ksd(cov0, cov1),
           ley = metric_ley(cov0, cov1),
           mhb = metric_mhb(cov0, cov1)
           )
}


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
##             PRIVATE FUNCTIONS
## ------------------------------------------------------------------------
## ------------------------------------------------------------------------

## overlapping coefficient
metric_ovl <- function(cov0, cov1) {
  cov <- c(cov0, cov1)
  if (length(unique(cov)) <= 10) {
      all_x <- c(rep(0, length(cov0)),
                 rep(1, length(cov1)))
      pt    <- apply(prop.table(table(cov, all_x), 2),
                     1, min)

    return(sum(pt))
  }

  mn <- min(cov) * 1.25;
  mx <- max(cov) * 1.25;

  f1 <- approxfun(density(cov1, from = mn, to = mx,
                          bw = "nrd"));
  f0 <- approxfun(density(cov0, from = mn, to = mx,
                          bw = "nrd"));

  fn <- function(x)
    pmin(f1(x), f0(x))

  s <- try(integrate(fn, lower = mn, upper = mx,
                     subdivisions = 500)$value)

  ifelse(inherits(s, "try-error"),
         NA,
         s)
}

## K-S distance
metric_ksd <- function(cov0, cov1) {
    cov    <- c(cov0, cov1);
    cdf_1  <- ecdf(cov1);
    cdf_0  <- ecdf(cov0);
    1/max(abs(cdf_1(cov) - cdf_0(cov)))
}

## Levy distance
metric_ley <- function(cov0, cov1) {
    cov   <- c(cov0, cov1);
    cdf_1 <- ecdf(cov1);
    cdf_0 <- ecdf(cov0);
    e     <- max(abs(cdf_1(cov) - cdf_0(cov)))

    if (length(unique(cov)) <= 10)
        return(e)

    x     <- seq(min(cov), max(cov), length.out = 1000);
    check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                 cdf_1(x) <= cdf_0(x + e) + e)

    while (check) {
        e <- e - .01
        check <- all(cdf_0(x - e) - e <= cdf_1(x) &
                     cdf_1(x) <= cdf_0(x + e) + e)
    }

    1/e
}

## mahalanobis balance
##
## covs should be a reduced datset that contains only those covariates
## that will be used for calculating Mahalanobis balance, for example,
## covs = dat[,1:6]
## trt should be the exposure variable,
## for example, trt=dat$X
##
metric_mhb <- function(cov0, cov1) {
  cov01 <- rbind(cov0, cov1)
  sinv  <- solve(cov(cov01))
  x0    <- colMeans(cov0)
  x1    <- colMeans(cov1)

  rst <- sum((t(x1 - x0) %*% sinv) * (x1 - x0))
  1/ rst
}
