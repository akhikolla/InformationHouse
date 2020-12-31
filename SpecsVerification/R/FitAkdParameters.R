#' Fit the 5 parameters used for affine kernel dressing by minimum CRPS estimation.
#'
#' @param ens a N*R matrix. An archive of R-member ensemble forecasts for N time instances.
#' @param obs a vector of length N. The verifying observations corresponding to the N ensemble forecasts.
#' @return The function returns a list of 5 parameters for affine kernel dressing. 
#'
#' @details
#' Affine Kernel Dressing transforms the discrete K-member forecast ensemble at time instance n, `ens[n, ]`, to a continuous distribution function for the target `y` by the equation:
#' 
#'       p(y|ens) = 1 / K * sum {dnorm(y, z.i, s)} \cr
#' where   s = (4/3/K)^0.4 * (s1 + s2 * a^2 * var(ens)) \cr
#' and   z.i = r1 + r2 * mean(ens) + a * ens
#' 
#' The parameters r1, r2, a, s1, s2 are fitted by minimizing the continuously ranked probability score (CRPS). The optimization is carried out using the R function `optim(...)`.
#' 
#' Since the evaluation of the CRPS is numerically expensive, the optimization can take a long time. Speed can be increased by optimizing the parameters only for a part of the forecast instances.
#'
#' @examples
#' data(eurotempforecast)
#' FitAkdParameters(ens, obs)
#' @seealso DressEnsemble, DressCrps, DressIgn, PlotDressedEns, GetDensity
#' @references
#' Broecker J. and Smith L. (2008). From ensemble forecasts to predictive distribution functions. Tellus (2008), 60A, 663--678. \doi{10.1111/j.1600-0870.2008.00333.x}.
#' @export

FitAkdParameters <- function(ens, obs) {

  #       p(y|x) = 1 / K * sum {dnorm(y, z.i(x), s(x))}
  # where   s(x) = (4/3/K)^0.4 * (s1 + s2 * a^2 * var(x))
  # and   z.i(x) = r1 + r2 * mean(x) + a * x[i]

  # ensemble means
  m.x <- rowMeans(ens, na.rm=TRUE)
  # ensemble variance 
  v.x <-  apply(ens, 1, var, na.rm=TRUE)
  # squared silverman factor 
  K <- max(rowSums(!is.na(ens)))
  sf.2 <- (4 / 3 / K) ^ 0.4

  # initial guesses
  m1 <- lm(obs ~ m.x)
  m2 <- lm(resid(m1) ^ 2 ~ v.x)
  coef1 <- as.vector(coef(m1))
  coef2 <- as.vector(coef(m2))

  r1 <- coef1[1]
  s1 <- coef2[1] / sf.2
  s2 <- 1
  a <- ifelse(test = coef2[2] <= 0,
               yes = 0,
               no = sqrt(coef2[2] / (1 + sf.2)))
  r2 <- coef1[2] - a

  # the objective function: mean crps 
  f <- function(parms) {
    parms <- as.list(parms)
    sigma2 <- with(parms, sf.2 * (s1 + s2*a*a*v.x))
    if (any(sigma2 < 0)) {
      ret <- Inf
    } else {
      d.ens <- DressEnsemble(ens, "akd", parms)
      ret <- mean(DressCrps(d.ens, obs))
    }
    ret
  }
  parms <- c(a=a, r1=r1, r2=r2, s1=s1, s2=s2)

  # if f cannot be evaluated at initial guesses, 
  # use standard silverman as initial guess
  if (!is.finite(f(parms))) { # catches Inf, NA, NaN
    parms <- c(a=1, r1=0, r2=0, s1=0, s2=1)
  }

  # optimize
  opt <- optim(par=parms, fn=f)
  
  # return
  opt[["par"]]
}

