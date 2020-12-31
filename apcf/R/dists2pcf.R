#' Convert Distances & Ratios to PCF
#'
#' Estimates the adapted Pair Correlation Function (PCF) of a pattern together
#' with a pointwise critical envelope based on distances and ratios calculated
#' by [pat2dists()].
#'
#' Since the pair-correlation function is a density function, we employ the
#' frequently used Epanechnikov kernel (Silverman, 1986; Stoyan and Stoyan,
#' 1994; Nuske et al., 2009). The Epanechnikov kernel is a weight function
#' putting maximal weight to pairs with distance exactly equal to *r* but also
#' incorporating pairs only roughly at distance *r* with reduced weight. This
#' weight falls to zero if the actual distance between the points differs from
#' *r* by at least \eqn{\delta}{delta}, the so-called bandwidth parameter,
#' which determines the degree of smoothness of the function. Penttinen et al.
#' (1992) and Stoyan and Stoyan (1994) suggest to set *c* aka stoyan-parameter
#' of \eqn{c / {\sqrt{\lambda}}}{c/sqrt(lambda)} between 0.1 and 0.2 with
#' \eqn{\lambda}{lambda} being the intensity of the pattern.
#'
#' The edge correction is based on suggestions by Ripley (1981). For each pair
#' of objects \eqn{i} and \eqn{j}, a buffer with buffer distance
#' \eqn{r_{ij}}{r_ij} is constructed around the object \eqn{i}. The object
#' \eqn{j} is then weighted by the inverse of the ratios \eqn{p_{ij}}{p_ij}
#' of the buffer perimeter being within the study area. That way we account for
#' the reduced probability of finding objects close to the edge of the study
#' area.
#'
#' The alpha level of the pointwise critical envelope is
#' \eqn{\alpha = \frac{n\_rank * 2}{n\_sim + 1}}{alpha = (n_rank * 2) / (n_sim + 1)}
#' according to (Besag and Diggle, 1977; Buckland, 1984; Stoyan and Stoyan, 1994).
#'
#'
#' @param dists An object of class [dists]. Usually created by [pat2dists()]
#' @param r A step size or a vector of values for the argument r at which g(r)
#'        should be evaluated.
#' @param r_max maximum value for the argument r.
#' @param kernel String. Choice of smoothing kernel (only the "epanechnikov"
#'        kernel is currently implemented).
#' @param stoyan Bandwidth coefficient (smoothing the Epanechnikov kernel).
#'        Penttinen et al. (1992) and Stoyan and Stoyan (1994) suggest values
#'        between 0.1 and 0.2.
#' @param n_rank Rank of the value amongst the n_sim simulated values
#'        used to construct the envelope. A rank of 1 means that the minimum
#'        and maximum simulated values will be used. Must be >= 1 and < n_sim/2.
#'        Determines together with `n_sim` in [pat2dists()] the alpha level of
#'        the envelope. If `alpha` and `n_sim` are fix, n_rank can be
#'        calculated by `(n_sim+1)*alpha/2` eg. `(199+1) * 0.05/2 = 5`
#'
#' @return An object of class [fv_pcf] containing the function values of the
#'   PCF and the envelope.
#'
#' @references
#' Besag, J. and Diggle, P.J. (1977): Simple Monte Carlo tests for spatial
#' pattern. Appl. Sta. 26, 327–333.
#'
#' Buckland, S.T. (1984). Monte Carlo Confidence Intervals. Biometrics, 40(3),
#' 811-817.
#'
#' Nuske, R.S., Sprauer, S. and Saborowski J. (2009): Adapting the
#' pair-correlation function for analysing the spatial distribution of
#' canopy gaps. Forest Ecology and Management (259): 107–116.
#'
#' Penttinen, A., Stoyan, D., Henttonen, H., 1992. Marked point-processes in
#' forest statistics. For. Sci. 38, 806–824.
#'
#' Ripley, B.D., 1981. Spatial Statistics. Wiley, New York.
#'
#' Silverman, B.W., 1986. Density Estimation for Statistics and Data Analysis.
#' Chapman and Hall, London.
#'
#' Stoyan, D. and Stoyan, H. (1994) Fractals, random shapes and point fields:
#' methods of geometrical statistics. John Wiley and Sons.
#'
#' @seealso [pat2dists()], [plot.fv_pcf()]
#'
#' @examples
#' # it's advised against setting n_sim < 199
#' ds <- pat2dists(area=system.file("shapes/sim_area.shp", package="apcf"),
#'                 pattern=system.file("shapes/sim_pat_reg.shp", package="apcf"),
#'                 max_dist=25, n_sim=3)
#'
#' # derive PCF and envelope
#' pcf <- dists2pcf(ds, r=0.2, r_max=25, stoyan=0.15, n_rank=1)
#'
#' @export
dists2pcf <- function(dists, r, r_max=NULL, kernel="epanechnikov", stoyan,
                      n_rank){
  # check parameter ---------------------------------------------------------
  if(!is.dists(dists))
    stop("dists must be an object of class dists")

  # n_rank
  if(n_rank %% 1 != 0 || n_rank < 1)
    stop("n_rank must be an integer and >= 1")
  # n_rank >= n_sim/2 is checked in C++ function do_env(), first time nsim is known

  # lambda
  if(is.null(area <- attr(dists, 'area', exact=TRUE)) |
     is.null(n_obj <- attr(dists, 'n_obj', exact=TRUE)))
       stop("attributes 'area' and 'n_obj' are missing")
  if(area <= 0)
    stop("Something went wrong: area of study area is <= 0")
  if(n_obj < 2)
    stop("Something went wrong: n_obj < 2")

  # r to vector of steps (rs)
  if(missing(r))
    stop("r must be given")
  if(is.null(max_dist <- attr(dists, 'max_dist', exact=TRUE)) & missing(r_max))
    stop("neither r_max nor the attribute 'max_dist' is available")
  else
    if(r_max > max_dist)
      stop("r_max must be <= max_dist")

  if(length(r) == 1){
    rs <- seq(from=r, to=r_max, by=r)
  } else {
    step <- unique(diff(r))
    if(length(step) != 1)
      stop("r must be equidistant")
    if(max(r) > r_max || min(r) < step){
      warning(paste0("reduced r to the range ", step, "-", r_max))
      rs <- r[r >= step & r <= r_max]
    } else {
      rs <- r
    }
  }

  # kernel
  if(tolower(kernel) != "epanechnikov")
    stop("only the Epanechnikov kernel is currently implemented")

  # create fv_pcf object ----------------------------------------------------
  out <- pcf_envelope(dists$sim, dists$dist, dists$ratio,
                      rs, area, n_obj, stoyan, n_rank)

  return(out)
}
