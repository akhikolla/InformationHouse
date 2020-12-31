##' Extract non-dominated points from a set, or with respect to a reference Pareto front
##' @title Generic non-domination computation
##' @param points matrix (one point per row) from which to extract non-dominated points, or,
##' if a reference \code{ref} is provided, non-dominated points with respect to  \code{ref}
##' @param ref matrix (one point per row) of reference (faster if they are already Pareto optimal)
##' @param return.idx if \code{TRUE}, return indices instead of points
##' @return Non-dominated points from \code{points}, unless a \code{ref} is provided, in which case return points from \code{points} non-dominated by \code{ref}.
##' If \code{return.idx} is \code{TRUE}, only returns indices
##' @export
##' @details Use Kung non-domination sorting
##' @references
##' Kung, H. T., Luccio, F., & Preparata, F. P. (1975). On finding the maxima of a set of vectors. Journal of the ACM (JACM), 22(4), 469-476.
##' @examples
##' \dontrun{
##' d <- 6
##' n <- 1000
##' n2 <- 1000
##'
##' test <- matrix(runif(d * n), n)
##' ref <- matrix(runif(d * n), n)
##' indPF <- nonDom(ref, return.idx = TRUE)
##' all(nonDom(ref) == ref[indPF,])
##'
##' system.time(res <- nonDom(test, ref[indPF,,drop = FALSE], return.idx = TRUE))
##'
##' res2 <- rep(NA, n2)
##' library(emoa)
##' t0 <- Sys.time()
##' for(i in 1:n2){
##'   res2[i] <- !is_dominated(t(rbind(test[i,, drop = FALSE], ref[indPF,])))[1]
##' }
##' print(Sys.time() - t0)
##'
##' all(res == which(res2))
##'
##' all(nonDom(test, ref) == test[res2,])
##'
##' }
nonDom <- function(points, ref = NULL, return.idx = FALSE){
  if(is.null(ref)){
    ordrs <- order(points[,1])
    if(return.idx) return(ordrs[nonDomInd_cpp(points[ordrs, , drop = FALSE])])
    return(points[ordrs[nonDomInd_cpp(points[ordrs, , drop = FALSE])],])
  }

  res <- nonDomSet(points, ref)
  if(return.idx) return(which(res))
  return(points[res, ,drop = FALSE])
}
