#' Auxiliary function for removing trends and mean from observation vector or ensemble matrix.
#'
#' @description Detrend fits a linear function to a time-series of observations or to the time-series of ensemble means of an ensemble matrix. The linear trend is removed, and if option demean is true, the total mean is removed as well.
#' @param x A vector, matrix, or data.frame.
#' @param demean logical; if true, the total mean is removed from x
#' @return The function returns an object of the same dimensions as the argument `x`, but with its linear trend and (possibly) its mean removed.
#' @examples
#'  data(eurotempforecast)
#'  Detrend(ens)
#'  Detrend(obs, demean=FALSE)
#' @export

Detrend <- function(x, demean=TRUE) {
  # detrend using the row means
  if (is(x, 'matrix')) {
    xx <- rowMeans(x, na.rm=TRUE)
  } else {
    xx <- x
  }
  N <- length(xx)
  if (N == 1) {
    # if only one time step is given, detrending amounts to substracting the mean
    trnd <- x
  } else {
    # otherwise estimate a linear function of time and take the fitted values
    # as "the trend" 
    t <- 1:N
    lmod <- lm(xx~t)
    trnd <- drop(cbind(1, t) %*% coef(lmod))
  }
  # if demean is false, add the grand mean back to x minus trend
  m <- ifelse(demean, 0, mean(unlist(x), na.rm=TRUE))
  return(x - trnd + m)
}

