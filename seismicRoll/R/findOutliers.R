#' @export
#' @title Outlier Detection with a Rolling Hampel Filter
#' @param x an \R numeric vector
#' @param n integer window size
#' @param thresholdMin initial value for outlier detection
#' @param selectivity value between [0-1] used in determining outliers, or \code{NA} if \code{fixedThreshold=TRUE}.
#' @param increment integer shift to use when sliding the window to the next location
#' @param fixedThreshold logical specifying whether outlier detection uses \code{selectivity}  (see below) 
#' @description A wrapper for the roll_hampel() function that counts outliers using either a user specified threshold
#' value or a threshold value based on the statistics of the incoming data.
#' @details   The \code{thresholdMin} level is similar to a sigma value for normally distributed data.
#' Hampel filter values above 6 indicate a data value that is extremely unlikely
#' to be part of a normal distribution  (~ 1/500 million) and therefore very likely to be an outlier. By
#' choosing a relatively large value for \code{thresholdMin} we make it less likely that we
#' will generate false positives. False positives can include high frequency environmental noise.
#' 
#' With the default setting of \code{fixedThreshold=TRUE} any value above the threshold is considered an outlier
#' and the \code{selectivity} is ignored.
#' 
#' The \code{selectivity} is a value between 0 and 1 and is used to generate an appropriate 
#' threshold for outlier detection based on the statistics of the incoming data. A lower value
#' for \code{selectivity} will result in more outliers while a value closer to 1.0 will result in 
#' fewer. If \code{fixedThreshold=TRUE}, \code{selectivity} may have a value of \code{NA}.
#' 
#' When the user specifies \code{fixedThreshold=FALSE}, the \code{thresholdMin} and \code{selectivity}
#' parameters work like squelch and volume on a CB radio: \code{thresholdMin} sets a noise threshold 
#' below which you don't want anything returned while \code{selectivity} adjusts the number of points defined as outliers
#' by setting a new threshold defined by the maximum value of \code{roll_hampel} multiplied by \code{selectivity}.
#' 
#' \code{n}, the windowSize, is a parameter that is passed to \code{roll_hampel()}.
#' 
#' \bold{The default value of \code{increment=1} should not be changed.} Outliers are defined
#' as individual points that stand apart from their neighbors. Applying the Hampel filter to
#' every other point by using \code{increment} > 1 will invariably miss some of the outliers.
#' 
#' @return A vector of indices associated with outliers in the incoming data \code{x}.
#' @seealso \code{\link{roll_hampel}}
#' @examples
#' # Noisy sinusoid with outliers
#' a <- jitter(sin(0.1*seq(1e4)),amount=0.2)
#' indices <- sample(seq(1e4),20)
#' a[indices] <- a[indices]*10
#' 
#' # Outlier detection should identify many of these altered indices
#' sort(indices)
#' findOutliers(a)

findOutliers <- function(x, n=41, thresholdMin=10, selectivity=NA, increment=1,
                         fixedThreshold=TRUE) {
  
  if ( !is.logical(fixedThreshold) ) {
    warning(paste0('findOutliers() argument fixedThreshold=',fixedThreshold,' is not a logical value.  Using default value.'))
    fixedThreshold <- TRUE
  }
  
  h <- roll_hampel(x, n, increment)
  
  #if 50%+ of values in window are the same, then h blows up to Inf. In this case, set to NA.

  h[is.infinite(h)] <- NA  

  if(all(is.na(h))) {
    stop("roll_hampel returns a vector with all NA or NaN (50%+ of values in all windows are identical)")
  }

  maxH <- max(h, na.rm=TRUE)
  
  if ( maxH < thresholdMin ) {
    
    # If no values cross thresholdMin, return an empty vector
    return(which(h > maxH)) # integer vector of length zero
    
  } else {
    
    # Some values cross thresholdMin.
    
    if ( fixedThreshold ) {
      # Pick outliers based on thresholdMin
      outliers <- which(h > thresholdMin)
      return(outliers)
    } else {
      # Pick outliers based on maxH and selectivity
      outliers <- which(h > maxH * selectivity)
      return(outliers)
    }
    
  }
  
}

