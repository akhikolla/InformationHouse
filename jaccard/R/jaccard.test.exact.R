#' Compute p-value using the exact solution
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients.
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#'
#' @return \code{jaccard.test.exact} returns a list consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @importFrom stats rbinom pchisq rnorm runif
#' @export jaccard.test.exact
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test.exact(x,y)
jaccard.test.exact <- function(x, y, px = NULL, py = NULL, verbose = TRUE) {
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
  degenerate<-FALSE
  
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
    
  if(!null.p) {
    pvalue <- jaccard_mca_rcpp_known_p(px,py,m,j.obs,0)$pvalue
  } else {
    pvalue <- jaccard_mca_rcpp(px,py,m,j.obs,0)$pvalue
  }

  return(
    list(
      statistics = j.obs,
      pvalue = pvalue,
      expectation = expectation
    )
  )
}