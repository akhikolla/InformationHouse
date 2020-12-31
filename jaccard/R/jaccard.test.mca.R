#' Compute p-value using the Measure Concentration Algorithm
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#' @param accuracy an error bound on approximating a multinomial distribution
#' @param error.type an error type on approximating a multinomial distribution ("average", "upper", "lower")
#'
#' @return \code{jaccard.test.mca} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @export jaccard.test.mca
#' @useDynLib jaccard
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.mca(x,y,accuracy = 1e-05)
jaccard.test.mca <- function(x, y, px = NULL, py = NULL, accuracy = 1e-05, error.type = "average", verbose = TRUE) {
  if (length(x) != length(y)) stop("Length mismatch")
  m <- length(x)
  null.p<-FALSE
  x <- as.logical(x)
  y <- as.logical(y)
  if (is.null(px) | is.null(py)) {
    px <- mean(x)
    py <- mean(y)
    null.p <- TRUE
  }

  x <- as.logical(x)
  y <- as.logical(y)
  degenerate <- FALSE

  expectation <- jaccard.ev(x, y, px=px, py=py)
  j.obs <- sum(x & y)/sum(x | y) - expectation

  if(px==1 | py==1 | sum(x) == length(x) | sum(y) == length(y)) {
    warning("One or both input vectors contain only 1's.")
    degenerate <- TRUE
  }
  if(px==0 | py==0 | sum(x) == 0 | sum(y) == 0) {
    warning("One or both input vectors contain only 0's")
    degenerate <- TRUE
  }
  if(degenerate) {
    return(list(statistics = 0, pvalue = 1, expectation = expectation))
  }

  if (!null.p) {
    tan = jaccard_mca_rcpp_known_p(px,py,m,j.obs,accuracy)
    pvalue = tan$pvalue
  } else {
    tan = jaccard_mca_rcpp(px,py,m,j.obs,accuracy)
    pvalue = tan$pvalue
  }
  pvalue <- switch(error.type, lower = pvalue, average = pvalue/tan$accuracy,
                   upper = pvalue + 1 - tan$accuracy)

  return(
    list(
      statistics = j.obs,
      pvalue = pvalue,
      expectation = expectation,
      accuracy = 1 - tan$accuracy,
      error.type = error.type
      )
  )
}
