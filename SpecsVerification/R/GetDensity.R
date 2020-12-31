#' Calculate density and integrated density function of a dressed ensemble forecast at a matrix of values
#'
#' @param dressed.ens A list returned by the function `DressEnsemble`. See `?DressEnsemble` for details.
#' @param x A matrix with either 1 row or nrow(dressed.ens[["ens"]]) rows and an arbitrary number of columns, holding the arguments at which the forecast distributions are to be evaluated. See Details.
#' @param integrated logical, (default=FALSE): If `integrated` is TRUE, the integrated density (i.e. the value of the cumulative distribution function) is returned, otherwise the value of the density is returned.
#' @return The function returns a matrix, whose rows correspond to the individual ensemble forecasts and whose columns correspond to the values provided by the argument `x`.
#' @details If you want to evaluate each forecast distribution function at the same x-values, a matrix with one row can be provided, e.g. `x = matrix(c(-1, 0, 1), nrow=1)`
#'
#' If the N individual forecast distributions are to be evaluated at different x-values, a matrix with N rows must be provided, where N is the number of time instances. 
#'
#' To calculate the PIT values for the dressed ensemble and observations `obs`, use `GetDensity(dressed.ens, x = matrix(obs, ncol=1), integrated=TRUE)`
#' @examples
#'  data(eurotempforecast)
#'  dressed.ens <- DressEnsemble(ens)
#'  # calculate each density at the same x-values
#'  x1 <- matrix(seq(-3, 3, 0.1), nrow=1)
#'  dens1 <- GetDensity(dressed.ens, x1)
#'  # get the densities that the forecast 
#'  # distributions assign to the observations
#'  x2 <- matrix(obs, ncol=1)
#'  dens2 <- GetDensity(dressed.ens, x2)
#'  # get the integrated densities that the forecast 
#'  # distributions assign to the observations (useful
#'  # for constructing a PIT histogram)
#'  pit <- GetDensity(dressed.ens, x2, integrated=TRUE)
#' @seealso DressEnsemble, DressCrps, DressIgn, PlotDressedEns, FitAkdParameters
#' @export

GetDensity <- function(dressed.ens, x, integrated=FALSE) {

  stopifnot(is.matrix(x))

  if (integrated) {
    kernel.fun <- pnorm
  } else {
    kernel.fun <- dnorm
  }

  N <- nrow(dressed.ens[["ens"]])
  if (nrow(x) == 1) {
    x <- matrix(rep(x, N), nrow=N, byrow=TRUE)
  }
  stopifnot(nrow(x) == N)

  d <-  
  sapply(1:ncol(x), function(j) {
    sapply(1:nrow(x), function(i) {
      with(dressed.ens, 
        mean(kernel.fun(x[i,j], mean=ens[i, ], sd=ker.wd[i, ]), na.rm=TRUE)
      )
    })
  })

  # return
  d
}


