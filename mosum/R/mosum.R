#' MOSUM procedure for multiple change-point estimation
#' 
#' Computes the MOSUM detector, detects (multiple) change-points and estimates their locations.
#' @param x input data (a \code{numeric} vector or an object of classes \code{ts} and \code{timeSeries})
#' @param G an integer value for the moving sum bandwidth;
#' \code{G} should be less than \code{length(n)/2}.
#' Alternatively, a number between \code{0} and \code{0.5} describing the moving sum bandwidth
#' relative to \code{length(x)} can be given
#' @param G.right if \code{G.right != G}, the asymmetric bandwidth \code{(G, G.right)} will be used;
#' if \code{max(G, G.right)/min(G, G.right) > 4}, a warning message is generated
#' @param var.est.method how the variance is estimated;
#' possible values are
#' \itemize{
#'    \item{\code{"mosum"}}{both-sided MOSUM variance estimator}
#'    \item{\code{"mosum.min"}}{minimum of the sample variance estimates from the left and right summation windows}
#'    \item{\code{"mosum.max"}}{maximum of the sample variance estimates from the left and right summation windows}
#'    \item{\code{"custom"}}{a vector of \code{length(x)} is to be parsed by the user; use \code{var.custom} in this case to do so}
#' }
#' @param var.custom a numeric vector (of the same length as \code{x}) containing
#' local estimates of the variance or long run variance; use iff \code{var.est.method = "custom"}
#' @param boundary.extension a logical value indicating whether the boundary
#' values should be filled-up with CUSUM values  
#' @param threshold string indicating which threshold should be used to determine significance.
#' By default, it is chosen from the asymptotic distribution at the given significance level \code{alpha}.
#' Alternatively it is possible to parse a user-defined numerical value with \code{threshold.custom}
#' @param alpha a numeric value for the significance level with
#' \code{0 <= alpha <= 1}; use iff \code{threshold = "critical.value"}
#' @param threshold.custom a numeric value greater than 0 for the threshold of significance;
#' use iff \code{threshold = "custom"}
#' @param criterion string indicating how to determine whether each point \code{k} at which MOSUM statistic 
#' exceeds the threshold is a change-point; possible values are
#' \itemize{
#'    \item{\code{"eta"}}{there is no larger exceeding in an \code{eta*G} environment of \code{k}}
#'    \item{\code{"epsilon"}}{\code{k} is the maximum of its local exceeding environment, which has at least size \code{epsilon*G}}
#' }
#' @param eta a positive numeric value for the minimal mutual distance of 
#' changes, relative to moving sum bandwidth (iff \code{criterion = "eta"})
#' @param epsilon a numeric value in (0,1] for the minimal size of exceeding
#' environments, relative to moving sum bandwidth (iff \code{criterion = "epsilon"})
#' @param do.confint flag indicating whether to compute the confidence intervals for change-points
#' @param level use iff \code{do.confint = TRUE}; a numeric value (\code{0 <= level <= 1}) with which
#' \code{100(1-level)\%} confidence interval is generated
#' @param N_reps use iff \code{do.confint = TRUE}; number of bootstrap replicates to be generated
#' @return S3 object of class \code{mosum.cpts}, which contains the following fields:
#'    \item{x}{input data}
#'    \item{G.left, G.right}{left and right summation bandwidths}
#'    \item{var.est.method, var.custom,boundary.extension}{input parameters}
#'    \item{stat}{a series of MOSUM statistic values; the first \code{G} and last \code{G.right} values are \code{NA} iff \code{boundary.extension = FALSE}}
#'    \item{rollsums}{a series of MOSUM detector values; equals \code{stat*sqrt(var.estimation)}}
#'    \item{var.estimation}{the local variance estimated according to \code{var.est.method}}
#'    \item{threshold, alpha, threshold.custom}{input parameters}
#'    \item{threshold.value}{threshold value of the corresponding MOSUM test}
#'    \item{criterion, eta, epsilon}{input parameters}
#'    \item{cpts}{a vector containing the estimated change-point locations}
#'    \item{cpts.info}{data frame containing information about change-point estimators including detection bandwidths, asymptotic p-values for the corresponding MOSUM statistics and (scaled) size of jumps}
#'    \item{do.confint}{input parameter}
#'    \item{ci}{S3 object of class \code{cpts.ci} containing confidence intervals for change-points iff \code{do.confint=TRUE}}
#' @references A. Meier, C. Kirch and H. Cho (2019)
#' mosum: A Package for Moving Sums in Change-Point Analysis. \emph{To appear in the Journal of Statistical Software}.
#' @references B. Eichinger and C. Kirch (2018)
#' A MOSUM procedure for the estimation of multiple random change-points.
#' \emph{Bernoulli}, Volume 24, Number 1, pp. 526-564.
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
#' m <- mosum(x, G = 40)
#' plot(m)
#' summary(m)
#' @export
mosum <- function(x, G, G.right=G, 
                      var.est.method=c('mosum', 'mosum.min', 'mosum.max', 'custom')[1], 
                      var.custom=NULL, boundary.extension=TRUE, 
                      threshold=c('critical.value', 'custom')[1], alpha=.1, 
                      threshold.custom=NULL, criterion=c('eta', 'epsilon')[1], 
                      eta=0.4, epsilon=0.2, do.confint=FALSE, 
                      level=0.05, N_reps=1000) {
  
  # Consistency checks on input
  stopifnot(alpha >= 0 && alpha <= 1)
  stopifnot(criterion=='epsilon' || criterion=='eta')
  stopifnot(criterion!='epsilon' || epsilon >= 0)
  stopifnot(criterion!='eta' || eta >= 0)
  stopifnot(!do.confint || N_reps>0)

  n <- length(x)
  m <- mosum.stat(x, G, G.right, var.est.method, var.custom, boundary.extension)
  
  G.left <- m$G.left
  G.right <- m$G.right
  G.min <- min(G.right, G.left)
  G.max <- max(G.right, G.left)
  K <- G.min / G.max
  changePoints <- numeric(0)
  
  if (threshold == 'critical.value' & G.max/G.min > 4) {
    warning('Bandwidths are too unbalanced, \n (G, G.right) satisfying max(G, G.right)/min(G, G.right) <= 4 is recommended')
  }
  
  if (threshold == 'critical.value') {
    threshold_val <- mosum.criticalValue(n, G.left, G.right, alpha)
  } else if (threshold == 'custom') {
    threshold_val <- threshold.custom
  } else {
    stop('threshold must be either \'critical.value\' or \'custom\'')
  }
  
  # get exceeding TRUE/FALSE vector
  exceedings <- (m$stat > threshold_val)
  
  # adjust, in case of no boundary CUSUM extension
  if (!m$boundary.extension) {
    exceedings[n-G.right+1] <- FALSE
  }
  
  if (criterion=='epsilon') {
    # get number of subsequent exceedings
    exceedingsCount <- (exceedings) * unlist(lapply(rle(exceedings)$lengths, seq_len))
    # get exceeding-intervals of fitting length
    minIntervalSize <- max(1, (G.min+G.max) / 2 * epsilon)
    intervalEndPoints <- which(diff(exceedingsCount) <= -minIntervalSize)
    intervalBeginPoints <- intervalEndPoints - exceedingsCount[intervalEndPoints] + 1
    if (!m$boundary.extension) {
      # manually adjust right border
      if (exceedings[n-G.right] && !((n-G.right) %in% intervalEndPoints)) {
        lastBeginPoint <- n - G.right - exceedingsCount[n-G.right] + 1
        stopifnot(exceedings[seq(lastBeginPoint,n-G.right)])
        stopifnot(!(lastBeginPoint %in% intervalBeginPoints))
        highestStatPoint <- which.max(m$stat[seq(lastBeginPoint,n-G.right)]) + lastBeginPoint - 1
        if (highestStatPoint-lastBeginPoint >= minIntervalSize/2) {
          # print(paste0('Found change-point at the right border (G=(', G.left, ',', G.right, ')).'))
          intervalEndPoints <- c(intervalEndPoints, n-G.right)
          intervalBeginPoints <- c(intervalBeginPoints, lastBeginPoint)
        }
      }
      # manually adjust left border
      if (exceedings[G.left] && !(G.left %in% intervalBeginPoints)) {
        firstEndPoint <- which(diff(exceedingsCount) < 0)[1]
        stopifnot(exceedings[seq(G.left,firstEndPoint)])
        stopifnot(!(firstEndPoint %in% intervalEndPoints))
        highestStatPoint <- which.max(m$stat[seq(G.left,firstEndPoint)]) + G.left - 1
        if (firstEndPoint - highestStatPoint >= minIntervalSize/2) {
          # print(paste0('Found change-point at the left border (G=(', G.left, ',', G.right, ')).'))
          intervalEndPoints <- c(firstEndPoint, intervalEndPoints)
          intervalBeginPoints <- c(G.left, intervalBeginPoints)
        }
      }
    }
    numChangePoints <- length(intervalBeginPoints)
    if (numChangePoints > 0) {
      for (i in 1:numChangePoints) {
        changePoint <- intervalBeginPoints[i] + which.max(m$stat[(intervalBeginPoints[i]):(intervalEndPoints[i])]) - 1
        changePoints <- c(changePoints, changePoint)
      }
    }
  } else { # (criterion=='eta')
    localMaxima <- (c((diff.default(m$stat) < 0), NA) & c(NA, diff.default(m$stat) > 0))
    # adjust, in case of no boundary CUSUM extension
    if (!m$boundary.extension) {
      localMaxima[n-G.right] <- TRUE
    }
    p.candidates <- which(exceedings & localMaxima)
    changePoints <- eta_criterion_help(p.candidates, m$stat, eta, G.left, G.right)
  }
  est.cpts.info <- data.frame(cpts = changePoints, 
                              G.left =  rep(G.left, length(changePoints)), 
                              G.right =  rep(G.right, length(changePoints)),
                              p.value = mosum.pValue(m$stat[changePoints], n, G.left, G.right),
                              jump = sqrt((G.left+G.right)/G.left/G.right)*m$stat[changePoints])

  ret <- structure(
    list(x=x,
        G.left=G.left,
        G.right=G.right,
        var.est.method=m$var.est.method,
        var.custom=m$var.custom, 
        boundary.extension=m$boundary.extension,
        stat=m$stat,
        rollsums=m$rollsums,
        var.estimation=m$var.estimation,
        threshold = threshold,
        alpha=alpha,
        threshold.custom = threshold.custom,
        threshold.value=threshold_val,
        criterion=criterion,
        eta=eta,        
        epsilon=epsilon,
        cpts=changePoints,
        cpts.info=est.cpts.info,
        do.confint=FALSE,
        ci=NA), # note
    class='mosum.cpts')
  if (do.confint) {
    ret$ci <- confint.mosum.cpts(ret, level=level, N_reps=N_reps)
    ret$do.confint <- TRUE
  } 
  ret
}

#' Plotting the output from MOSUM procedure
#' 
#' Plotting method for S3 objects of class \code{mosum.cpts}
#' @method plot mosum.cpts
#' @param x a \code{mosum.cpts} object
#' @param display which to be plotted against the change-point estimators; possible values are
#' \itemize{
#'    \item{\code{"data"}}{input time series is plotted along with the estimated piecewise constant signal}
#'    \item{\code{"mosum"}}{scaled MOSUM detector values are plotted}
#' }
#' @param cpts.col a specification for the color of the vertical lines at
#' the change-point estimators, see \link[graphics]{par}
#' @param critical.value.col a specification for the color of the horizontal line
#' indicating the critical value, see \link[graphics]{par}; use iff \code{display = "mosum"}
#' @param xlab graphical parameter
#' @param ... additional graphical arguments, see \link[graphics]{plot}
#' and \link[graphics]{abline}
#' @details 
#' The location of each change-point estimator is plotted as a vertical line
#' against the input time series and the estimated piecewise constant signal (\code{display = "data"})
#' or MOSUM detector values (\code{display = "mosum"}).
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
#' m <- mosum(x, G = 40)
#' par(mfrow = c(2, 1), mar = c(2.5, 2.5, 2.5, .5))
#' plot(m, display = "data")
#' plot(m, display = "mosum")
#' @importFrom graphics abline lines plot
#' @export
plot.mosum.cpts <- function(x, display=c('data', 'mosum')[1], cpts.col='red', critical.value.col='blue', xlab='Time', ...) {
  if (class(x$x)=='ts') {
    x_plot <- as.numeric(time(x$x))
  } else if(class(x$x) == 'timeSeries') {
    x_plot <- time(x$x)
  } else {
    x_plot <- seq_len(length(x$x))
  }
  if (display == 'mosum'){ 
    plot(x=x_plot, y=x$stat, type='l', xlab=xlab, ylab = '', main = 'mosum', ...)
    abline(h=x$threshold.value, col=critical.value.col)
  }
  if (display == 'data'){
    brks <- c(0, x$cpts, length(x$x))
    fhat <- x$x * 0
    for(kk in 1:(length(brks) - 1)){
      int <- (brks[kk] + 1):brks[kk + 1]
      fhat[int] <- mean(x$x[int])
    }
    plot(x=x_plot, y = x$x, type='l', xlab=xlab, ylab = '', main = 'input data', ...)
    lines(x=x_plot, y = fhat, col = 'darkgray', type = 'l', lwd = 2)
  }
  for (p in x$cpts) {
    abline(v=x_plot[p], col='red', ...)
  }
}

#' Summary of change-points estimated by MOSUM procedure
#' 
#' Summary method for objects of class \code{mosum.cpts}
#' @method summary mosum.cpts
#' @param object a \code{mosum.cpts} object
#' @param ... not in use
#' @details Provide information about each estimated change-point, 
#' including the bandwidths used for its estimation, associated p-value and (scaled) jump size;
#' if \code{object$do.confint=TRUE}, end points of the pointwise and uniform confidence intervals
#' are also provided.
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
#' m <- mosum(x, G = 40, do.confint = TRUE)
#' summary(m)
#' @export
summary.mosum.cpts <- function(object, ...) { 
  n <- length(object$x)
  if(length(object$cpts) > 0){
    ans <- object$cpts.info
    ans$p.value <- signif(ans$p.value, 3)
    ans$jump <- round(ans$jump, 3)
  } 
  if(object$do.confint) ans <- cbind(ans, object$ci$CI[, -1, drop=FALSE])

#  cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated at alpha = ', object$alpha, ' according to ', object$criterion, '-criterion', sep=''))
  if(object$criterion=='eta') cat(paste('\n with eta = ', object$eta, sep=''))
  if(object$criterion=='epsilon') cat(paste('\n with epsilon = ', object$epsilon, sep=''))
  cat(paste(' and ', object$var.est.method, ' variance estimate:', sep=''))
  cat('\n')
  cat('\n')
  if(length(object$cpts) > 0) print(ans, print.gap = 3) else cat('no change-point is found') 
  cat('\n')
}

#' Change-points estimated by MOSUM procedure
#' 
#' Print method for objects of class \code{mosum.cpts}
#' @method print mosum.cpts
#' @param x a \code{mosum.cpts} object
#' @param ... not in use
#' @examples 
#' x <- testData(lengths = rep(100, 3), means = c(0, 5, -2), sds = rep(1, 3), seed = 1234)$x
#' m <- mosum(x, G = 40)
#' print(m)
#' @export
print.mosum.cpts <- function(x, ...) {
  #cat(paste('created using mosum version ', utils::packageVersion('mosum'), sep=''))
  cat(paste('change-points estimated with bandwidths (', x$G.left, ', ', x$G.right, ')', 
            ' at alpha = ', x$alpha, sep='')) 
  cat(paste('\n according to ', x$criterion, '-criterion', sep=''))
  if(x$criterion=='eta') cat(paste(' with eta = ', x$eta, sep=''))
  if(x$criterion=='epsilon') cat(paste(' with epsilon = ', x$epsilon, sep=''))
  cat(paste(' and ', x$var.est.method, ' variance estimate:', sep=''))
  cat('\n')
  cat('\n')
  cat('  ')
  if(length(x$cpts)==0) cat('no change-point is found') else cat(x$cpts)
  cat('\n')
}