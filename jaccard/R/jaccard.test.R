#' Test for Jaccard/Tanimoto similarity coefficients
#'
#' Compute statistical significance of Jaccard/Tanimoto similarity coefficients between binary vectors, using four different methods.
#'
#' There exist four methods to compute p-values of Jaccard/Tanimoto similarity coefficients:
#' \code{mca}, \code{bootstrap}, \code{asymptotic}, and \code{exact}. This is simply a wrapper function for
#' corresponding four functions in this package: \link{jaccard.test.mca}, \link{jaccard.test.bootstrap},  \link{jaccard.test.asymptotic}, and \link{jaccard.test.exact}.
#'
#' We recommand using either \code{mca} or \code{bootstrap} methods,
#' since the \code{exact} solution is slow for a moderately large vector and \code{asymptotic} approximation may be inaccurate depending on the input vector size.
#' The {bootstrap} method uses resampling with replacement binary vectors to compute a p-value (see optional arguments).
#' The \code{mca} method uses the measure concentration algorithm that estimates the multinomial distribution with a known error bound (specified by an optional argument \code{accuracy}).
#'
#' @section Optional arguments for \code{method="bootstrap"}:
#' \describe{
#'   \item{fix}{whether to fix (i.e., not resample) \code{x} and/or \code{y}}
#'   \item{B}{a total bootstrap iteration}
#'   \item{seed}{a seed for a random number generator}
#' }
#' @section Optional arguments for \code{method="mca"}:
#' \describe{
#'   \item{accuracy}{an error bound on approximating a multinomial distribution}
#'   \item{error.type}{an error type on approximating a multinomial distribution (\code{"average"}, \code{"upper"}, \code{"lower"})}
#'   \item{seed}{a seed for the random number generator.}
#' }
#'
#' @param x a binary vector (e.g., fingerprint)
#' @param y a binary vector (e.g., fingerprint)
#' @param method a method to compute a p-value (\code{"mca"}, \code{"bootstrap"}, \code{"asymptotic"}, or \code{"exact"})
#' @param px probability of successes in \code{x} (optional)
#' @param py probability of successes in \code{y} (optional)
#' @param verbose whether to print progress messages
#' @param ... optional arguments for specific computational methods
#'
#' @return \code{jaccard.test} returns a list mainly consisting of
#' \item{statistics}{centered Jaccard/Tanimoto similarity coefficient}
#' \item{pvalue}{p-value}
#' \item{expectation}{expectation}
#'
#' @export jaccard.test
#' @seealso \link{jaccard.test.bootstrap} \link{jaccard.test.mca} \link{jaccard.test.exact} \link{jaccard.test.asymptotic}
#'
#' @examples
#' set.seed(1234)
#' x = rbinom(100,1,.5)
#' y = rbinom(100,1,.5)
#' jaccard.test(x,y,method="bootstrap")
#' jaccard.test(x,y,method="mca")
#' jaccard.test(x,y,method="exact")
#' jaccard.test(x,y,method="asymptotic")
jaccard.test <- function(x, y, method="mca", px = NULL, py = NULL, verbose = TRUE, ...) {
  # length of fingerprints
  if(length(x) != length(y)) stop("Length mismatch.")

  if(method=="mca") out <- jaccard.test.mca(x, y, px, py, verbose, ...)
  if(method=="bootstrap") out <- jaccard.test.bootstrap(x, y, px, py, verbose, ...)
  if(method=="exact") out <- jaccard.test.exact(x, y, px, py, verbose)
  if(method=="asymptotic") out <- jaccard.test.asymptotic(x, y, px, py, verbose)

  return(out)
}

#' Pair-wise tests for Jaccard/Tanimoto similarity coefficients
#'
#' Given a data matrix, it computes pair-wise Jaccard/Tanimoto similarity coefficients
#' and p-values among rows (variables). For fine controls, use \code{"jaccard.test"}.
#'
#' @param dat a data matrix
#' @param method a method to compute a p-value (\code{"mca"}, \code{"bootstrap"}, \code{"asymptotic"}, or \code{"exact"})
#' @param compute.qvalue whether to compute q-values
#' @param verbose whether to print progress messages
#' @param ... optional arguments for specific computational methods
#'
#' @return \code{jaccard.test.pairwise} returns a list of matrices
#' \item{statistics}{Jaccard/Tanimoto similarity coefficients}
#' \item{pvalues}{p-values}
#' \item{qvalues}{q-values}
#'
#' @import qvalue
#' @export jaccard.test.pairwise
#' @seealso \link{jaccard.test}
jaccard.test.pairwise <- function(dat, method="mca", verbose = TRUE, compute.qvalue = TRUE, ...) {
  m <- dim(dat)[1]
  n <- dim(dat)[2]

  pvalues <- matrix(NA,nrow=m,ncol=m)
  statistics <- matrix(NA,nrow=m,ncol=m)
  expectation <- matrix(NA,nrow=m,ncol=m)
  for(i in 1:m) {
    if(verbose) print(i)
    if(i < m) {
      for(j in (i+1):m) {
        out <- jaccard.test(dat[i,], dat[j,], method=method, verbose=verbose, ...)
        pvalues[i,j] <- out$pvalue
        statistics[i,j] <- out$statistics
        expectation[i,j] <- out$expectation
      }
    }
  }

  qvalues <- matrix(0,nrow=m,ncol=m)
  if(compute.qvalue) {
    qvalues <- qvalue(pvalues)
    results <- list(statistics=statistics, pvalues=pvalues, expectation=expectation, qvalues=qvalues)
  } else {
    results <- list(statistics=statistics, pvalues=pvalues, expectation=expectation)
  }

  return(results)
}
