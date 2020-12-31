#' Confidence intervals for change-points
#' 
#' Generate bootstrap confidence intervals for change-points.
#' @param object an object of class \code{mosum.cpts}
#' @param parm specification of which parameters are to be given confidence intervals; \code{parm = "cpts"} is supported
#' @param level numeric value in (0, 1), such that the \code{100(1-level)\%} confidence bootstrap intervals are computed
#' @param N_reps number of bootstrap replications
#' @param ... not in use
#' @return S3 object of class \code{cpts.ci}, containing the following fields:
#'    \item{level, N_reps}{input parameters}
#'    \item{CI}{data frame of five columns, 
#'    containing the estimated change-points (column \code{cpts}),
#'    the pointwise confidence intervals 
#'    (columns \code{pw.left} and \code{pw.right})
#'    and the uniform confidence intervals 
#'    (columns \code{unif.left} and \code{unif.right}) for the corresponding change-points}
#' @details See the referenced literature for further details
#' @references A. Meier, C. Kirch and H. Cho (2019)
#' mosum: A Package for Moving Sums in Change-point Analysis. \emph{To appear in the Journal of Statistical Software}.
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 3, 1), sds = rep(1, 3), seed = 1337)$x
#' m <- mosum(x, G = 40)
#' ci <- confint(m, N_reps = 5000)
#' print(ci$CI)
#' @importFrom Rcpp evalCpp
#' @importFrom stats confint
#' @useDynLib mosum, .registration = TRUE
#' @method confint mosum.cpts
#' @export
confint.mosum.cpts <- function(object, parm = "cpts", level=0.05, N_reps=1000, ...) {
  stopifnot(class(object) == 'mosum.cpts')
  stopifnot(parm=='cpts')
  if (object$do.confint) {
    return(object$ci)
  } else {
    G.left <- object$G.left
    G.right <- object$G.right
    cpts <- object$cpts
    q <- length(cpts)
    cpts.info = matrix(NA, ncol=3, nrow=q)
    cpts.info[, 1] <- cpts
    cpts.info[, 2] <- G.left
    cpts.info[, 3] <- G.right
    cpts_bootstrap(cpts.info, N_reps, level, object)
  }
}

#' Confidence intervals for change-points
#' 
#' Generate bootstrap confidence intervals for change-points.
#' @param object an object of class \code{multiscale.cpts}
#' @param parm specification of which parameters are to be given confidence intervals; \code{parm = "cpts"} is supported
#' @param level numeric value in (0, 1), such that the \code{100(1-level)\%} confidence bootstrap intervals are computed
#' @param N_reps number of bootstrap replications
#' @param ... not in use
#' @return S3 object of class \code{cpts.ci}, containing the following fields:
#'    \item{level, N_reps}{input parameters}
#'    \item{CI}{data frame of five columns, 
#'    containing the estimated change-points (column \code{cpts}),
#'    the pointwise confidence intervals 
#'    (columns \code{pw.left} and \code{pw.right})
#'    and the uniform confidence intervals 
#'    (columns \code{unif.left} and \code{unif.right}) for the corresponding change-points}
#' @details See the referenced literature for further details
#' @references A. Meier, C. Kirch and H. Cho (2019)
#' mosum: A Package for Moving Sums in Change-point Analysis. \emph{To appear in the Journal of Statistical Software}.
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 3, 1), sds = rep(1, 3), seed = 1337)$x
#' mlp <-  multiscale.localPrune(x, G = c(8, 15, 30, 70))
#' ci <- confint(mlp, N_reps = 5000)
#' print(ci$CI)
#' @importFrom Rcpp evalCpp
#' @importFrom stats confint
#' @useDynLib mosum, .registration = TRUE
#' @method confint multiscale.cpts
#' @export
confint.multiscale.cpts <- function(object, parm='cpts', level=0.05, N_reps=1000, ...) {
  stopifnot(class(object)=='multiscale.cpts')
  stopifnot(parm=='cpts')
  if (object$do.confint) {
    return(object$ci)
  } else {
    cpts_bootstrap(object$cpts.info, N_reps, level, object)
  }
}

#' Helping/wrapper fuction for C++ calls
#' @importFrom stats quantile
#' @keywords internal
cpts_bootstrap <- function(cpts_info, N_reps, level, mcpts) {
  x <- mcpts$x
  q <- nrow(cpts_info)
  if (q==0) cpts <- pointwise.left <- pointwise.right <- uniform.left <- uniform.right <- integer(0)
  if (q>0) {
    n <- length(x)
    cpts <- cpts_info[,1]
    # Get bootstrap replicates from C++
    tmp <- cpts_bootstrap_help(as.matrix(cpts_info), x, N_reps)
    k_star <- tmp$k_star
    # Pointwise confidence intervals
    C_value_j <- apply(abs(tmp$k_star1), 2, quantile, 1-level/2)
    pointwise.left <- ceiling(pmax(pmax(1, cpts-cpts_info[, 2]+1), cpts - C_value_j))
    pointwise.right <- floor(pmin(pmin(n, cpts+cpts_info[, 3]), cpts + C_value_j))
    # Uniform confidence intervals
    uCI_help <- apply(abs(tmp$k_star2), 1, max)
    C_value <- quantile(uCI_help, 1-level)
    uniform.left <- uniform.right <- rep(NA, q)
    for (j in seq_len(q)) {
      uniform.left[j] <- max(max(1, cpts[j] - cpts_info[j, 2]+1),
                      cpts[j] - C_value * tmp$sigma2_hat[j] / tmp$d_hat[j]^2)
      uniform.right[j] <- min(min(n, cpts[j] + cpts_info[j, 3]),
                      cpts[j] + C_value * tmp$sigma2_hat[j] / tmp$d_hat[j]^2)
    }
  }
  uniform.left <- ceiling(uniform.left)
  uniform.right <- floor(uniform.right)
  CI <- data.frame(cpts=cpts, 
                   pw.left=pointwise.left,
                   pw.right=pointwise.right,
                   unif.left=uniform.left,
                   unif.right=uniform.right)
  return(structure(list(N_reps=N_reps,
                        level=level,
                        CI=CI), class='cpts.ci'))
}
