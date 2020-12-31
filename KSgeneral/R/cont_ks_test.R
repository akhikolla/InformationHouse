#################################################
## This function computes the p-value of one sample two-sided Kolmogorov-Smirnov
## test when the distribution under the null hypothesis is continuous

cont_ks_test <- function(x, y, ...)
{
  DNAME <- deparse(substitute(x))
  x <- x[!is.na(x)]
  n <- length(x)
  if (n < 1L){
    stop("not enough 'x' data")
  }

  PVAL <- NULL

  if (is.character(y)){
    y <- get(y, mode = "function", envir = parent.frame())
  }
  if (!is.function(y)){
    stop("'y' must be a function or a string naming a valid function")
  }
  METHOD <- "One-sample Kolmogorov-Smirnov test"
  TIES <- FALSE

  if (length(unique(x)) < n){
    stop("ties should not be present for the continuous Kolmogorov-Smirnov test")
  }

  x <- y(sort(x), ...) - (0 : (n - 1))/n
  STATISTIC <- max(c(x, 1/n - x))

  PVAL <- KSgeneral::cont_ks_c_cdf(STATISTIC, n)

  nm_alternative <- "two-sided"

  names(STATISTIC) <- "D"

  RVAL <- list(statistic = STATISTIC, p.value = PVAL, alternative = nm_alternative, method = METHOD, data.name = DNAME)

  class(RVAL) <- "htest"

  return(RVAL)

}
