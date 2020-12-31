#' MOSUM statistic
#' 
#' Computes the statistical values for the MOSUM test for changes in the mean.
#' @param x input data (\code{numeric} vector or object of class \code{ts})
#' @param G an integer value for the length of the moving sum window; 
#' \code{G} should be less than \code{length(n)/2}.
#' Alternatively a number between \code{0} and \code{0.5} describing the moving sum bandwidth
#' relative to \code{length(x)}. 
#' @param G.right iff \code{!is.na(G.right)}, the asymmetric bandwidth (G,G.right) will be used
#' @param var.est.method how the variance is estimated;
#' possible values are
#' \itemize{
#'    \item{\code{'custom'}}{a vector of \code{length(x)} is to be parsed by the user; use \code{var.custom} in this case to to so}
#'    \item{\code{'mosum'}}{both-sided MOSUM variance estimator}
#'    \item{\code{'mosum.min'}}{minimum of the sample variance estimates from the left and right summation windows}
#'    \item{\code{'mosum.max'}}{maximum of the sample variance estimates from the left and right summation windows}   
#' }
#' @param var.custom a numeric vector (of the same length as \code{x}) containing
#' local estimates of the variance or long run variance; use iff \code{var.est.method=custom}
#' @param boundary.extension a logical value indicating whether the boundary
#' values should be filled-up with CUSUM values
#' @return S3 \code{mosum.stat} object, which contains the following fields:
#'    \item{x}{the numeric input vector provided}
#'    \item{G.left,G.right}{left and right bandwidths}
#'    \item{var.est.method,var.custom,boundary.extension}{input parameters}
#'    \item{stat}{a series of MOSUM statistic values; the first \code{G} and last \code{G.right} values are \code{NA} iff \code{boundary.extension=FALSE}}
#'    \item{rollsums}{a series of MOSUM detector values; equals \code{stat*sqrt(var.estimation)}}
#'    \item{var.estimation}{the local variance estimated according to \code{var.est.method}}
#' @details This class only contains the values for the MOSUM statistic.
#' For statistical evaluation and change-point extraction, use \link[mosum]{mosum}.
#' See also \link[mosum]{multiscale.bottomUp} and \link[mosum]{multiscale.localPrune}.
#' @importFrom Rcpp evalCpp
#' @useDynLib mosum, .registration = TRUE
#' @keywords internal 
mosum.stat <- function(x, G, G.right=NA, var.est.method='mosum', 
                  var.custom=NULL, boundary.extension=TRUE) {
  n <- length(x)
  symmetric <- is.na(G.right)
  
  if(G < 1 & G >= 0.5) stop('Please use relative bandwidth between 0 and 0.5.')
  if(G >= n/2) stop('Please use bandwidth smaller than length(x)/2.')
  if(!symmetric){
    if(G.right < 1 & G.right >= 0.5) stop('Please use relative bandwidth between 0 and 0.5.')
    if(G.right >= n/2) stop('Please use bandwidth smaller than length(x)/2.')
  }
  
  abs.bandwidth <- (G>=1)
  if (!abs.bandwidth) {
    G <- floor(n * G)
    if (!symmetric) G.right <- floor(n * G.right)
  }
  
  # Consistency checks on input
  stopifnot(NCOL(x) == 1) 
  stopifnot(class(x)=='ts' || class(x)=='numeric' || class(x) == 'timeSeries')
  stopifnot(G > 0 && G < n)
  stopifnot(symmetric || !is.na(G.right))
  stopifnot(symmetric || (G.right > 0 && G.right < n))
  
  G.left <- G
  if (symmetric) G.right <- G
  G.min <- min(G.right, G.left)
  G.max <- max(G.right, G.left)
  K <- G.min / G.max
  
  # Calculate value of statistic.
  sums.left <- rolling_sum(x, G.left)
  if (G.left == G.right) {
    sums.right <- sums.left
  } else {
    sums.right <- rolling_sum(x, G.right)
  }
  unscaledStatistic <- c(rep(NA,G.left-1), (G.min/G.right*sums.right[(G.left+1):n] - G.min/G.left*sums.left[1:(n - G.left)]), NA) / (sqrt((K+1)*G.min))
  
  # Calculate variance estimation.
  if (!is.null(var.custom) && var.est.method != 'custom') {
    stop('Please use var.est.method = custom when parsing var.custom.')
  }
  if (var.est.method == 'custom') {
    if (is.null(var.custom)) {
      stop('Expecting var.custom to be not NULL for var.est.method=custom')
    }
    if (length(var.custom) != n) {
      stop('Expecting var.custom to be of length n = length(x)')
    }
    var <- var.custom
  } else if (var.est.method == 'global') {
    # Note: This is Deprecated
    var <- rep((sum(x^2) - (sum(x)^2)/n)/n,n)
  } else { # MOSUM-based variance estimators
    summedSquares.left <- rolling_sum(x^2, G.left) # zoo::rollsum(x^2, k=G.left, fill=NA, align='left')
    squaredSums.left <- sums.left^2
    var.tmp.left <- summedSquares.left[1:(n-G.left+1)] - 1/G.left*(squaredSums.left[1:(n-G.left+1)])
    var.left <- c(rep(NA,G.left-1), var.tmp.left) / G.left
    if (G.left == G.right) {
      summedSquares.right <- summedSquares.left
      squaredSums.right <- squaredSums.left
      var.tmp.right <- var.tmp.left
    } else {
      summedSquares.right <- rolling_sum(x^2, G.right) # zoo::rollsum(x^2, k=G.right, fill=NA, align='left')
      squaredSums.right <- sums.right^2
      var.tmp.right <- summedSquares.right[1:(n-G.right+1)] - 1/G.right*(squaredSums.right[1:(n-G.right+1)])
    }
    var.right <- c(var.tmp.right[2:(n-G.right+1)], rep(NA,G.right)) / G.right
    if (var.est.method == 'mosum') {
      var <- (var.left + var.right) / 2
    } else if (var.est.method == 'mosum.left') {
      # Note: This is Deprecated
      var <- var.left
    } else if (var.est.method == 'mosum.right') {
      # Note: This is Deprecated
      var <- var.right
    } else if (var.est.method == 'mosum.min') {
      var <- pmin(var.left, var.right)
    } else if (var.est.method == 'mosum.max') {
      var <- pmax(var.left, var.right)
    } else {
      stop('unknown variance estimation method')
    }
  }
  var.estimation <- var

  # CUSUM extension to boundary
  if (boundary.extension) {
    if (n > 2*G.left) {
      weights.left <- sqrt( (G.left+G.right) / (1:G.left) / ((G.left+G.right-1):(G.right))) 
      unscaledStatistic[1:G.left] <- cumsum(mean(x[1:(G.left+G.right)])-x[1:G.left]) * 
        weights.left
      var.estimation[1:G.left] <- var.estimation[G.left]
    }
    if (n > 2*G.right) {
      weights.right <- sqrt( (G.left+G.right) / ((G.right - 1):0) / ((G.left + 1):(G.left + G.right))) 
      x.rev <- x[(n - G.left - G.right + 1):n]
      unscaledStatistic[(n-G.right+1):n] <- cumsum(mean(x.rev) - x.rev)[-(1:G.left)] * weights.right
      unscaledStatistic[n] <- 0
      var.estimation[(n-G.right+1):n] <- var.estimation[n-G.right]
    }
  }
  res <- abs(unscaledStatistic) / sqrt(var.estimation)
  structure(
    list(x=x,
         G.left=G.left,
         G.right=G.right,
         var.est.method=var.est.method,
         var.custom=var.custom, 
         boundary.extension=boundary.extension,
         stat=res,
         rollsums=unscaledStatistic,
         var.estimation=var.estimation),
    class='mosum.stat')
}

#' Plotting MOSUM statistics
#'
#' Plotting method for objects of class 'mosum.stat'.
#' @method plot mosum.stat
#' @param m a \code{mosum.stat} object
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}
#' @param critical.value.col a specification for the color of the
#' critical value, see \link[graphics]{par}
#' @param xlab graphical parameter
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' @importFrom stats plot.ts time
#' @keywords internal
plot.mosum.stat <- function(m, alpha=0.05, critical.value.col='blue', 
                       xlab='Time', ...) {
  if (class(m$x)=='ts') {
    x_plot <- as.numeric(time(m$x))
  } else {
    x_plot <- seq_len(length(m$x))
  }
  plot(x=x_plot, y=m$stat, type='l', xlab=xlab, ...)
  G.left <- m$G.left
  G.right <- m$G.right

  abline(h=mosum.criticalValue(n=length(m$x), G.left=G.left,
                               G.right=G.right, alpha=alpha),
         col=critical.value.col, ...)
}
