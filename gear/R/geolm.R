#' Linear model for geostatistical data.
#'
#' \code{geolm} creates a geostatistical linear model object
#' of the appropriate class based on the arguments,
#' especially the \code{cmod} arguments.
#'
#' Note: for the multiresolution Gaussian process model, if
#' \code{cmod$est == "f"} (i.e., if the nugget is finescale
#' instead of measurement error), then the \code{weights}
#' argument is internally set to \code{rep(1, n)}, where
#' \code{n} is the number of observations.
#'
#' \code{formula} should be specified after the form \code{y
#' ~ x1 + x2}, where \code{y} is the response variable and
#' \code{x1} and \code{x2} are the covariates of interest.
#' If \code{mu} is provided, the variables to the right of
#' \code{~} are ignored.
#'
#' @param formula An object of class
#'   \code{\link[stats]{formula}} providing a symbolic
#'   description of the model to be fitted.  See Details of
#'   this function and \code{\link[stats]{lm}}.
#' @param data A data frame containing the response,
#'   covariates, and location coordinates.
#' @param coordnames The columns of \code{data} containing
#'   the spatial coordinates, provided as a formula (e.g.,
#'   \code{~ x + y}), column numbers (e.g., \code{c(1, 2)}),
#'   or column names (e.g., \code{c("x", "y")})
#' @param mod A model object produced by one
#'   of the \code{cmod_*} functions, e.g.,
#'   \code{\link[gear]{cmod_std}}.
#' @param weights An optional vector of weights for the
#'   errors to be used in the fitting process.  A vector
#'   that is proportional to the reciprocal variances of the
#'   errors, i.e., errors are assumed to be uncorrelated
#'   with variances \code{evar/weights}.  Default is
#'   \code{NULL}, meaning that the weights are uniformly 1.
#' @param mu A single numeric value indicating the consant
#'   mean of the spatial process if simple kriging is
#'   desired.  Default is \code{NULL}, meaning that ordinary
#'   or universal kriging should be used.
#' @param cmod Retained for backwards compatibility. A model object produced by one
#'   of the \code{cmod_*} functions, e.g.,
#'   \code{\link[gear]{cmod_std}}.
#' @inheritParams evgram
#' @inheritParams cmod_std
#'
#' @return Returns a \code{geolm_*} object, where \code{*}
#'   depends on \code{mod}.
#'
#' @author Joshua French
#' @export
#' @examples
#' data = data.frame(y = rnorm(10), x1 = runif(10),
#'                  x2 = runif(10))
#' d = geodist(data[,c("x1", "x2")])
#' mod = cmod_man(v = exp(-d), evar = 1)
#' gearmod = geolm(y ~ x1, data = data,
#'                 coordnames = ~ x1 + x2, mod = mod)
geolm = function(formula, data, coordnames,
                 mod,
                 weights = NULL,
                 mu = NULL, longlat = NULL,
                 cmod = NULL) {
  call = match.call()
  if (!is.null(cmod)) {
    mod = cmod
    warning("The cmod argument has been replaced by the mod argument. Please update your code. The cmod argument will be deprecated in the future.")
  }
  if (!is.null(longlat)) {
    mod$longlat = longlat
    warning("The longlat argument shoud now be specified through the mod object. Please update your code. The longlat argument will be deprecated in the future.")
  }
  arg_check_geolm(formula, data, mod, weights, mu)
  if (!is.null(longlat)) {
    warning("The longlat argument should now be specified via the appropriate cmod* or vmod* function. Please update your code. The longlat argument in geolm will be depcrecated in the future.")
    arg_check_longlat(longlat)
    message("Attempting to update mode with longlat specificiation")
    mod$longlat = longlat
  }

  # ensure coordinates are available
  data_colnames = names(data)
  coordnames = arg_check_coordnames(coordnames, data_colnames)
  # determine coordinates for observed data
  coords = as.matrix(data[,coordnames])

  # create model frame
  mf = stats::model.frame(formula, data)
  # extract response
  y = c(stats::model.response(mf))
  n = length(y) # number of observations

  x = NULL # assume simple kriging
  # if not doing simple kriging
  if (is.null(mu)) {
    # create design matrix, remove all attributes
    # to make it a basic matrix
    x = stats::model.matrix(formula, data = mf)
    coeff_names = colnames(x)
    dimnames(x) = NULL
    attr(x, "assign") = NULL
  }

  # create default weights, if NULL
  if (is.null(weights)) weights = rep(1, n)

  geolm_fit(mod = mod, x = x, y = y, coords = coords,
            weights = weights, formula = formula, mu = mu,
            coordnames = coordnames, n = n, call = call,
            coeff_names = coeff_names)
  # if(!is.null(cmod))
  # {
  #   if(class(cmod) == "cmodStd"){
  #     return(geolm.Std(x, y, coords, mu, cmod, weights, formula, coordnames, n))
  #   }else if(class(cmod) == "cmodMan"){
  #     return(geolm.Man(x, y, coords, mu, cmod, weights, formula, coordnames, n))
  #   }else if(class(cmod) == "cmodMpp"){
  #     return(geolm.Mpp(x, y, coords, mu, cmod, weights, formula, coordnames, n))
  #   }else if(class(cmod) == "cmodSre"){
  #     return(geolm.Sre(x, y, coords, mu, cmod, weights, formula, coordnames, n))
  #   }else if(class(cmod) == "cmodMgp"){
  #     return(geolm.Mgp(x, y, coords, mu, cmod, weights, formula, coordnames, n))
  #   }
  # }else{
  #   stop("geolm currently only implemented when cmod is provided")
  # }
}

# check arguments of geolm
arg_check_geolm = function(formula, data, mod, weights, mu) {
  arg_check_geolm_formula(formula)
  if (!is.data.frame(data)) {
    stop("data should be a data frame")
  }
  valid_mod = c("cmodMan", "cmodStd")
  if (!is.element(class(mod), valid_mod)) {
    stop(paste("mod must be of one of the following classes:", paste(valid_mod, collapse = ", ")))
  }
  if (!is.null(weights)) {
    arg_check_weights(weights, nrow(data))
  }
  if (!is.null(mu)) {
    arg_check_mu(mu)
  }
}
