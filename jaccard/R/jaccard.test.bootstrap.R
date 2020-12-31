#' Compute p-value using the bootstrap procedure
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#' @param fix whether to fix (i.e., not resample) \code{x} and/or \code{y}
#' @param B a total bootstrap iteration
#' @param seed a seed for a random number generator
#'
#' @return \code{jaccard.test.bootstrap} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @importFrom qvalue empPvals
#' @export jaccard.test.bootstrap
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.bootstrap(x,y,B=500)
jaccard.test.bootstrap <- function(x, y, px = NULL, py = NULL, verbose=TRUE, fix="x", B=1000, seed=NULL) {
  if(!is.null(seed)) set.seed(seed)
  # length of fingerprints
  if(length(x) != length(y)) stop("Length mismatch")
  m <- length(x)

  # probabilities of ones
  if(is.null(px) | is.null(py)){
    px <- mean(x)
    py <- mean(y)
  }

  expectation <- jaccard.ev(x, y, px=px, py=py)
  j.obs <- jaccard(x, y, center=TRUE, px=px, py=py)
  degenerate <- FALSE
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

  j.null <- vector("numeric",B)
  if(verbose) message("Bootstrap Procedures : ")
  for(i in 1:B) {
    if(verbose & (i %% floor(B/10)) == 0) {
      message(paste0(i, " "),appendLF=FALSE)
    }

    if(fix == "x") {
      j.null[i] <- jaccard(x, sample(y, replace=TRUE), center=TRUE)
    } else if(fix == "y") {
      j.null[i] <- jaccard(sample(x, replace=TRUE), y, center=TRUE)
    } else {
      j.null[i] <- jaccard(sample(x, replace=TRUE), sample(y, replace=TRUE), center=TRUE)
    }
  }

  pvalue <- getp(abs(j.obs), abs(j.null))

  return(
    list(
      statistics = j.obs,
      statistics.null = j.null,
      pvalue = pvalue,
      expectation = expectation)
  )
}
