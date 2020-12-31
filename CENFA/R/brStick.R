#' Broken-stick method for detection of significant factors
#'
#' This function provides a simple way to determine the number of significant
#' factors in a factor analysis. This is done by comparing the eigenvalues of
#' each factor with those expected from a broken-stick distribution.
#'
#' @param eigs numeric. Vector of eigenvalues
#'
#' @examples
#' mod1 <- enfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' brStick(s.factor(mod1))
#'
#' @return Returns the number of significant factors.
#'
#' @references Jackson, Donald A. "Stopping rules in principal components analysis:
#'   a comparison of heuristical and statistical approaches." Ecology 74.8 (1993):
#'   2204-2214.
#'
#' @export

brStick <- function (eigs) {
  if(max(Im(eigs)) > 1e-5) stop("broken-stick method does not work for complex eigenvalues")
  eigs <- Re(eigs)
  p <- length(eigs)
  a <- NULL
  r <- NULL

  for(j in 1:p){
    a[j] <- 1/p * sum(1/(j:p) )
    r[j] <- eigs[j]/(sum(eigs))
  }
  length(which(r > a))
}

