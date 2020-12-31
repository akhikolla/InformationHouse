################################################################################
#' Mean squared or absolute \eqn{h}-step ahead prediction errors
#'
#' The function \code{MSPE} computes the empirical mean squared prediction
#' errors for a collection of \eqn{h}-step ahead, linear predictors
#' (\eqn{h=1,\ldots,H}) of observations \eqn{X_{t+h}}, where
#' \eqn{m_1 \leq t+h \leq m_2}, for two indices \eqn{m_1} and \eqn{m_2}.
#' The resulting array provides
#' \deqn{\frac{1}{m_{\rm lo} - m_{\rm up} + 1} \sum_{t=m_{\rm lo}}^{m_{\rm up}} R_{(t)}^2,}
#' with \eqn{R_{(t)}} being the prediction errors
#' \deqn{R_t := | X_{t+h} - (X_t, \ldots, X_{t-p+1}) \hat v_{N,T}^{(p,h)}(t) |,}
#' ordered by magnitude; i.e., they are such that \eqn{R_{(t)} \leq R_{(t+1)}}.
#' The lower and upper limits of the indices are
#' \eqn{m_{\rm lo} := m_1-h + \lfloor (m_2-m_1+1) \alpha_1 \rfloor} and
#' \eqn{m_{\rm up} := m_2-h - \lfloor (m_2-m_1+1) \alpha_2 \rfloor}.
#' The function \code{MAPE} computes the empirical mean absolute prediction
#' errors
#' \deqn{\frac{1}{m_{\rm lo} - m_{\rm up} + 1} \sum_{t=m_{\rm lo}}^{m_{\rm up}} R_{(t)},}
#' with \eqn{m_{\rm lo}}, \eqn{m_{\rm up}} and \eqn{R_{(t)}} defined as before.  
#'
#' @name measure-of-accuracy
#' @aliases MSPE
#' @export
#' 
#' @param X the data \eqn{X_1, \ldots, X_T}
#' @param predcoef the prediction coefficients in form of a list of an array
#' 				\code{coef}, and two integer vectors \code{t}	and \code{N}. The two
#' 			  integer vectors provide the information for which indices \eqn{t} and
#' 				segment lengths \eqn{N} the coefficients are to be interpreted;
#'        \code{(m1-H):(m2-1)} has to be a subset of \code{predcoef$t}.
#' 				if not provided the necessary coefficients will be computed using
#' 				\code{\link{predCoef}}.
#' @param m1 first index from the set in which the indices \eqn{t+h} shall lie
#' @param m2 last index from the set in which the indices \eqn{t+h} shall lie
#' @param P maximum order of prediction coefficients to be used;
#'          must not be larger than \code{dim(predcoef$coef)[1]}.
#' @param H maximum lead time to be used;
#'          must not be larger than \code{dim(predcoef$coef)[3]}.
#' @param N vector with the segment sizes to be used, 0 corresponds to
#' 					using 1, ..., t;
#'          has to be a subset of predcoef$N.
#' @param trimLo percentage \eqn{\alpha_1} of lower observations to be trimmed away
#' @param trimUp percentage \eqn{\alpha_2} of upper observations to be trimmed away
#' 
#' @return \code{MSPE} returns an object of type \code{MSPE} that has \code{mspe},
#' 	       an array of size \code{H}\eqn{\times}\code{P}\eqn{\times}\code{length(N)},
#'         as an attribute, as well as the parameters \code{N}, \code{m1},
#'         \code{m2}, \code{P}, and \code{H}.
#'         \code{MAPE} analogously returns an object of type \code{MAPE} that
#'         has \code{mape} and the same parameters as attributes.
#'
#' @example inst/examples/mspe.R
################################################################################
MSPE <- function(X, predcoef, m1 = length(X)/10, m2 = length(X), P = 1, H = 1,
                    N = c(0, seq(P + 1, m1 - H + 1)), trimLo = 0, trimUp = 0) {
  
  
  if (missing(predcoef)) {
    predcoef <- predCoef(X, P, H, (m1-H):(m2-1), N)
  }

  if (!is.list(predcoef)) {
    stop("predcoef needs to be a list, as specified in the documentation.")
  }
  
  if (all(attributes(predcoef)$names != c("coef", "t", "N"))) {
    stop("predcoef needs to be a list, as specified in the documentation.")
  }

  if (length(dim(predcoef$coef)) !=  5) {
    stop("predcoef$coef needs to be an array, as specified in the documentation.")
  }
  
  if (!is.numeric(predcoef$t)) {
    stop("predcoef$t needs to be a vector of integers.")
  }
  
  if (!is.numeric(predcoef$N)) {
    stop("predcoef$N needs to be a vector of integers.")
  }
  
  if (dim(predcoef$coef)[1] !=  dim(predcoef$coef)[2]) {
    stop("dim(predcoef$coef)[1] needs to be = dim(predcoef$coef)[2].")
  }
  
  if (dim(predcoef$coef)[1] < P) {
    stop("dim(predcoef$coef)[1] needs to be >= P.")
  }
  
  if (dim(predcoef$coef)[3] < H) {
    stop("dim(predcoef$coef)[3] needs to be >= H.")
  }
  
  if (dim(predcoef$coef)[4] != length(predcoef$t)) {
    stop("dim(predcoef$coef)[4] needs to be = length(predcoef$t).")
  }
  
  if (dim(predcoef$coef)[5] != length(predcoef$N)) {
    stop("dim(predcoef$coef)[5] needs to be = length(predcoef$N).")
  }
  
  if (!all(predcoef$N %in% N)) {
    stop("predcoef$N needs to be a subset of N.")
  }
  
  if (!all((m1-H):(m2-1) %in% predcoef$t)) {
    stop("(m1-H):(m2-1) has to be a subset of predcoef$t.")
  }

  mspe <- array( 0, dim = c(H, P, length(N)) )                
  T <- length(X) 
  
  
  N.len <- length(predcoef$N)
  N.idx <- which(is.element(predcoef$N, N))
  for (h in 1:H) {
    t.idx <- which(is.element(predcoef$t, (m1-h):(m2-h)))
    coef_a <- predcoef$coef[1:P, 1:P, 1:H, t.idx, N.idx, drop=FALSE]
    mspe[h, , ] <- computeMSPEcpp(X, coef_a, h, (m1-h):(m2-h), 1, trimLo, trimUp ) # 1 is for MSPE
  }
  
  value <- list(mspe=mspe, N=N, m1=m1, m2=m2, P=P, H=H) 
  class(value) <- "MSPE"
  return(value)

}

#' @name measure-of-accuracy
#' @aliases MAPE
MAPE <- function(X, predcoef, m1 = length(X)/10, m2 = length(X), P = 1, H = 1,
    N = c(0, seq(P + 1, m1 - H + 1)), trimLo = 0, trimUp = 0) {
  
  
  if (missing(predcoef)) {
    predcoef <- predCoef(X, P, H, (m1-H):(m2-1), N)
  }
  
  if (!is.list(predcoef)) {
    stop("predcoef needs to be a list, as specified in the documentation.")
  }
  
  if (all(attributes(predcoef)$names != c("coef", "t", "N"))) {
    stop("predcoef needs to be a list, as specified in the documentation.")
  }
  
  if (length(dim(predcoef$coef)) !=  5) {
    stop("predcoef$coef needs to be an array, as specified in the documentation.")
  }
  
  if (!is.numeric(predcoef$t)) {
    stop("predcoef$t needs to be a vector of integers.")
  }
  
  if (!is.numeric(predcoef$N)) {
    stop("predcoef$N needs to be a vector of integers.")
  }
  
  if (dim(predcoef$coef)[1] !=  dim(predcoef$coef)[2]) {
    stop("dim(predcoef$coef)[1] needs to be = dim(predcoef$coef)[2].")
  }
  
  if (dim(predcoef$coef)[1] < P) {
    stop("dim(predcoef$coef)[1] needs to be >= P.")
  }
  
  if (dim(predcoef$coef)[3] < H) {
    stop("dim(predcoef$coef)[3] needs to be >= H.")
  }
  
  if (dim(predcoef$coef)[4] != length(predcoef$t)) {
    stop("dim(predcoef$coef)[4] needs to be = length(predcoef$t).")
  }
  
  if (dim(predcoef$coef)[5] != length(predcoef$N)) {
    stop("dim(predcoef$coef)[5] needs to be = length(predcoef$N).")
  }
  
  if (!all(predcoef$N %in% N)) {
    stop("predcoef$N needs to be a subset of N.")
  }
  
  if (!all((m1-H):(m2-1) %in% predcoef$t)) {
    stop("(m1-H):(m2-1) has to be a subset of predcoef$t.")
  }
  
  mape <- array( 0, dim = c(H, P, length(N)) )                
  T <- length(X) 
  
  
  N.len <- length(predcoef$N)
  N.idx <- which(is.element(predcoef$N, N))
  for (h in 1:H) {
    t.idx <- which(is.element(predcoef$t, (m1-h):(m2-h)))
    coef_a <- predcoef$coef[1:P, 1:P, 1:H, t.idx, N.idx, drop=FALSE]
    mape[h, , ] <- computeMSPEcpp(X, coef_a, h, (m1-h):(m2-h), 2, trimLo, trimUp ) # 2 is for MAPE
  }
  
  value <- list(mape=mape, N=N, m1=m1, m2=m2, P=P, H=H) 
  class(value) <- "MAPE"
  return(value)
  
}


################################################################################
#' Plot a \code{MSPE} or \code{MAPE} object
#'
#' The function \code{plot.MSPE} plots a \code{MSPE} object that is returned by
#' the \code{MSPE} function.
#' The function \code{plot.MAPE} plots a \code{MAPE} object that is returned by
#' the \code{MAPE} function.
#' 
#' @name plot.measure-of-accuracy
#' @aliases plot.MSPE
#' 
#' @export
#' 
#' @importFrom graphics abline lines plot
#' 
#' @param x      The \code{MSPE} or \code{MAPE} object to be plotted.
#' @param vr     parameter to plot a line at level \code{vr}.
#'               Intended to be used to plot the mean squared prediction error
#'               of the trivial, null predictor; optional.
#' @param h      Defines for which \eqn{h}-step predictor the mean squared
#'               prediction errors will be shown; default: 1. 
#' @param N_min  If specified, the mean squared prediction errors with
#'               \eqn{N < N_{\rm min}} will not be shown; integer and optional.
#' @param legend Flag to specify if a legend, indicating which colour of the
#'               lines corresponds to which \eqn{p}, will be shown;
#'               default: \code{TRUE}.
#' @param display.mins Flag to specify if the minima for each \eqn{p}, and the
#'                     minimum accross \eqn{N=0} will be highlighted. 
#' @param add.for.legend add this much extra space for the legend, right of the
#'                       lines.
#' @param ...    Arguments to be passed to the underlying plot method
#' 
#' @seealso \code{\link{MSPE}}, \code{\link{MAPE}}
#' 
#' @return Returns the plot, as specified. 
################################################################################
plot.MSPE <- function(x, vr = NULL, h = 1, N_min = 1, legend = TRUE, display.mins = TRUE, add.for.legend = 0, ...) {
  M <- x$mspe
  N <- x$N
  p <- x$P
  
  if (is.null(vr)) {
    ylim.min <- min(M[h, , N == 0 | N >= N_min])
    ylim.max <- max(M[h, , N == 0 |  N >= N_min])
  } else {
    ylim.min <- min(vr, M[h, , N == 0 | N >= N_min])
    ylim.max <- max(vr, M[h, , N == 0 |  N >= N_min])
  }
    
  
  plot(x = N[N!=0  & N >= N_min], y = M[h, 1, N != 0 & N >= N_min],
      type="l", xlab="N", ylab = "MSPE", main = paste("h =", h),
      xlim=c(N_min, max(N)+add.for.legend),
      ylim=c(ylim.min, ylim.max)) 
  for (j in 2:p) {lines(x=N[N!=0  & N >= N_min], M[h, j, N != 0 & N >= N_min], col=j)}
  for (j in 1:p) {abline(h=M[h, j, N == 0], col = j, lty="dashed")}
  # lines(smooth.spline(x=N[N!=0], y=cppRes[1,2:(L+1)], spar=0.5), col="gray") #
  
  if (!is.null(vr)) {
    abline(h=vr, col="gray")
  }
  
  if (display.mins) {
    ## Find minimums
    idx1_s <- which(M[h, , N == 0] == min(M[h, , N == 0]), arr.ind = TRUE)[1]
    abline(h = M[h, idx1_s, N == 0], col = idx1_s, lty="dashed", lwd = 2)
    
    for (p in 1:x$P) {
      idx1_ls <- which(M[h, , N != 0 & N >= N_min] == min(M[h, , N != 0 & N >= N_min]), arr.ind = TRUE)[1,]
      idx1_ls_p <- which(M[h, p, N != 0 & N >= N_min] == min(M[h, p, N != 0 & N >= N_min]), arr.ind = TRUE)[1]
      abline(v = N[N != 0 & N >= N_min][idx1_ls_p], col = p, lty = "dotted")
    }
  }
  
  if (legend) {
    legend("topright", legend = paste("p =",1:p), lwd = 1, col=1:p) # , bg="white"
  }
  
}

#' @name plot.measure-of-accuracy
#' @aliases plot.MAPE
plot.MAPE <- function(x, vr = NULL, h = 1, N_min = 1, legend = TRUE, display.mins = TRUE, add.for.legend = 0, ...) {
  M <- x$mape
  N <- x$N
  p <- x$P
  
  if (is.null(vr)) {
    ylim.min <- min(M[h, , N == 0 | N >= N_min])
    ylim.max <- max(M[h, , N == 0 |  N >= N_min])
  } else {
    ylim.min <- min(vr, M[h, , N == 0 | N >= N_min])
    ylim.max <- max(vr, M[h, , N == 0 |  N >= N_min])
  }
  
  
  plot(x = N[N!=0  & N >= N_min], y = M[h, 1, N != 0 & N >= N_min],
      type="l", xlab="N", ylab = "MAPE", main = paste("h =", h),
      xlim=c(N_min, max(N)+add.for.legend),
      ylim=c(ylim.min, ylim.max)) 
  for (j in 2:p) {lines(x=N[N!=0  & N >= N_min], M[h, j, N != 0 & N >= N_min], col=j)}
  for (j in 1:p) {abline(h=M[h, j, N == 0], col = j, lty="dashed")}
  # lines(smooth.spline(x=N[N!=0], y=cppRes[1,2:(L+1)], spar=0.5), col="gray") #
  
  if (!is.null(vr)) {
    abline(h=vr, col="gray")
  }
  
  if (display.mins) {
    ## Find minimums
    idx1_s <- which(M[h, , N == 0] == min(M[h, , N == 0]), arr.ind = TRUE)[1]
    abline(h = M[h, idx1_s, N == 0], col = idx1_s, lty="dashed", lwd = 2)
    
    for (p in 1:x$P) {
      idx1_ls <- which(M[h, , N != 0 & N >= N_min] == min(M[h, , N != 0 & N >= N_min]), arr.ind = TRUE)[1,]
      idx1_ls_p <- which(M[h, p, N != 0 & N >= N_min] == min(M[h, p, N != 0 & N >= N_min]), arr.ind = TRUE)[1]
      abline(v = N[N != 0 & N >= N_min][idx1_ls_p], col = p, lty = "dotted")
    }
  }
  
  if (legend) {
    legend("topright", legend = paste("p =",1:p), lwd = 1, col=1:p) # , bg="white"
  }
  
}