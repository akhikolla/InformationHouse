#' Transform an ensemble forecast to a continuous forecast distribution by kernel dressing.
#'
#' @param ens a N*R matrix representing N time instances of real-valued R-member ensemble forecasts
#' @param dressing.method One of "silverman" (default), "akd", "akd.fit". See Details.
#' @param parameters A list, containing the parameters for the dressing method. See Details.
#' @return The function returns a list with elements `ens` (a N*R matrix, where ens[t,r] is the mean of the r-th kernel at time instance t) and `ker.wd` (a N*R matrix, where ker.wd[t,r] is the standard deviation of the r-th kernel at time t)
#' @details The dressing methods currently implemented and their required parameters are:
#' \describe{
#' \item{"silverman" (default)}{No parameters are given. At time instance `n` each ensemble member is replaced by a Gaussian kernel with mean ens[n, k] and variance (4 / 3 / K)^0.4 * var(ens[n, ]). This method is called "Silverman's rule of thumb" and provides a simple non-parametric method for smoothing a discrete ensemble.}
#' \item{"akd"}{Affine Kernel Dressing. The required parameters are list(r1, r2, a, s1, s2). The `k`-th ensemble member at time instance `n` is dressed with a Gaussian kernel with mean r1 + r2 * mean(ens[n,]) + a * ens[n, k] and variance (4 / 3 / K)^0.4 * (s1 + s2 * a^2 * var(ens[n,])). Negative variances are set to zero. Note that parameters = list(r1=0, r2=0, a=1, s1=0, s2=1) yields the same dressed ensemble as dressing.method="silverman".}
#' \item{"akd.fit"}{Affine Kernel Dressing with fitted parameters. The required parameters is list(obs), where `obs` is a vector of observations which are used to optimize the parameters r1, r2, a, s1, s2 by CRPS minimization. See ?FitAkdParameters for more information.}
#' }
#' 
#' @examples
#' data(eurotempforecast)
#' d.silverman <- DressEnsemble(ens)
#' d.akd <- DressEnsemble(ens, dressing.method="akd", 
#'                        parameters=list(r1=0, r2=0, a=1, 
#'                                        s1=0, s2=0))
#' d.akd.fit <- DressEnsemble(ens, dressing.method="akd.fit", 
#'                            parameters=list(obs=obs))
#' @seealso DressCrps, DressIgn, GetDensity, FitAkdParameters
#' @references 
#' Silverman, B.W. (1998). Density Estimation for Statistics and Data Analysis. London: Chapman & Hall/CRC. ISBN 0-412-24620-1.
#' Broecker J. and Smith L. (2008). From ensemble forecasts to predictive distribution functions. Tellus (2008), 60A, 663--678. \doi{10.1111/j.1600-0870.2008.00333.x}.
#' @export

DressEnsemble <- function(ens, dressing.method="silverman", 
                          parameters=NA) 
{

  # silverman's rule of thumb
  if (dressing.method == "silverman") {
    n.members <- rowSums(!is.na(ens))
    if (any(n.members==1)) {
      warning(c("Some ensembles have only one member. ",
                "The kernel width is set to zero for these."))
    }
    stdevs <- apply(ens, 1, sd, na.rm=TRUE)
    ker.wd <- stdevs * (4. / 3. / n.members) ^ 0.2

    K <- max(n.members)
    ker.wd <- matrix(rep(ker.wd, K), ncol=K)
  }

  #
  # affine kernel dressing
  #
  #       p(y|x) = 1 / K * sum {dnorm(y, z.i(x), s(x))}
  # where   s(x) = (4/3/K)^0.4 * (s1 + s2 * a^2 * var(x))
  # and   z.i(x) = r1 + r2 * mean(x) + a * x[i]
  #
  # dressing parameters a, r1, r2, s1, s2 are provided by the user

  if (dressing.method=="akd") {

    stopifnot(all(names(parameters) %in% c("a", "r1", "r2", "s1", "s2")))

    K <- max(rowSums(!is.na(ens)))
    v.x <- apply(ens, 1, var, na.rm=TRUE)
    m.x <- rowMeans(ens, na.rm=TRUE)
    sf.2 <- (4 / 3 / K) ^ 0.4

    z <- with(parameters, r1 + r2 * m.x  + a * ens)
    ens <- z

    s2 <- with(parameters, sf.2 * (s1 + s2 * a * a * v.x))
    s <- sapply(s2, function(x) sqrt(max(x, 0)))
    ker.wd <- matrix(rep(s, K), ncol=K)

  }


  #
  # affine kernel dressing as above, but parameter are estimated by optimizing
  # crps with respect to obs (which must be provided as parameters)
  #
  if (dressing.method=="akd.fit") {

    stopifnot(any(names(parameters) == "obs"))
    obs <- parameters[["obs"]]

    stopifnot(length(obs) == nrow(ens))

    parms <- as.list(FitAkdParameters(ens, obs))

    K <- max(rowSums(!is.na(ens)))
    v.x <- apply(ens, 1, var, na.rm=TRUE)
    m.x <- rowMeans(ens, na.rm=TRUE)
    sf.2 <- (4 / 3 / K) ^ 0.4

    z <- with(parms, r1 + r2 * m.x  + a * ens)
    ens <- z

    s2 <- with(parms, sf.2 * (s1 + s2 * a * a * v.x))
    s <- sapply(s2, function(x) sqrt(max(x, 0)))
    ker.wd <- matrix(rep(s, K), ncol=K)

  }

  # create object
  dressed.ens <- list(ens=ens, ker.wd=ker.wd)

  # return
  dressed.ens
}

