#' @rdname evaluate
#' @export
evaluate.cmodStd = function(mod, d, e = TRUE, f = TRUE) {
  arg_check_evaluate_cmodStd(mod, d, e, f)
  model = mod$model
  psill = mod$psill
  r = mod$r
  evar = mod$evar
  fvar = mod$fvar
  par3 = mod$par3
  longlat = mod$longlat
  angle = mod$angle
  ratio = mod$ratio
  radians = mod$radians
  # class(d) == "ganisoD" is needed because during optimization,
  # ratio may be set to 1 even when d is class ganisoD
  if (class(d)[1] == "ganisoD") {
    rmajor = r
    rminor = ratio * rmajor
    gangle = angle
    if (!radians) gangle = gangle * pi/180
    r = rmajor * rminor / sqrt(rminor^2 * cos(d$angles - gangle)^2 + rmajor^2 * sin(d$angles - gangle)^2)
    d = d$d
  }

  if (model == "exponential") {
    V = psill*exp(-d/r)

  } else if (model == "gaussian") {
    V = psill*exp(-(d/r)^2)
  } else if (model == "matern") {
    sd = d/r
    V = (d > 0) * psill*(2^(1 - par3)/gamma(par3)*sd^par3*besselK(sd, nu = par3))
    V[is.nan(V)] = psill
  } else if (model == "amatern") {
    sd = 2 * sqrt(par3) * d/r
    V = (d > 0) * psill*(2^(1 - par3)/gamma(par3)*sd^par3*besselK(sd, nu = par3))
    V[is.nan(V)] = psill
  } else if (model == "spherical") {
    sd = d/r
    V = psill*(1 - 1.5*sd + 0.5*(sd)^3)*(d < r)
  } else if (model == "wendland1") {
    sd = d/r
    V = psill*((1 - sd)^4 * (4*sd + 1))*(d < r)
  } else if (model == "wendland2") {
    sd = d/r
    V = psill*((1 - sd)^6 * (35*sd^2 + 18*sd + 3))/3*(d < r)
  } else if (model == "wu1") {
    sd = d/r
    V = psill*((1 - sd)^3 * (1 + 3*sd + sd^2))*(d < r)
  }else if (model == "wu2") {
    sd = d/r
    V = psill*((1 - sd)^4*(4 + 16*sd + 12*sd^2 + 3*sd^3))/4*(d < r)
  } else if (model == "wu3") {
    sd = d/r
    V = psill*((1 - sd)^6 * (1 + 6*sd + 41/3*sd^2 + 12*sd^3 + 5*sd^4 + 5/6*sd^5))*(d < r)
  }
  nugget = 0
  if (f & !e) {
    nugget = fvar
  } else if (e & !f) {
    nugget = evar
  } else if (e & f) {
    nugget = evar + fvar
  }
  return(V + (d == 0)*(nugget))
}

#' Argument check evaluate.cmodStd
#'
#' @param mod A cmodStd
#' @param d distnace matrix or ganisoD object if mod$ganiso != NULL
#' @param e Logical value. Measurement error?
#' @param f Logical value. Finescale error?
#' @noRd
arg_check_evaluate_cmodStd = function(mod, d, e, f) {
  if (mod$ratio < 1 && is.matrix(d)) {
    stop("d must be produced by the ganiso_d function when mod$ganiso is not NULL")
  } else if (is.matrix(d)) {
    arg_check_d(d)
  }
  arg_check_e(e)
  arg_check_f(f)
}

