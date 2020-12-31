

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% est_ability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Ability Estimation of an examinee
#' @description
#' \code{est_ability} estimates ability using on various methods such as
#' Owen's Bayesian estimation, Maximum Likelihood estimation,
#' Expected-a-Posteriori.
#'
#' @param ip An \code{\link{Item-class}}, \code{\link{Itempool-class}} or a
#'   \code{\link{Testlet-class}} object.
#' @param resp A vector containing examinee responses. If there are missing
#'   responses, they will not be included in the ability estimation.
#' @param method The method that will be used to estimate the ability.
#'   The default value is \code{"eap"}.
#'
#'   Current methods are:
#'   \describe{
#'     \item{\strong{\code{'sum_score'}}}{Basic sum (raw) score of
#'       responses.}
#'     \item{\strong{\code{'owen'}}}{Owen's Bayesian Ability Estimation.
#'
#'       This estimation method can be used only for dichotomous IRT
#'       models, 'Rasch', '1PL', '2PL', '3PL' and '4PL'.
#'
#'       Formulas were implemented in Owen (1975) and Vale (1977).  Original
#'       formulation does not contain D parameter. If \code{D = 1} original
#'       solution will be obtained. If \code{D = 1.7}  the \code{a} parameter
#'       will be multiplied with this number.
#'
#'       }
#'     \item{\strong{\code{'ml'}}}{Maximum Likelihood Ability Estimation
#'       via Newton-Raphson Algorithm}
#'     \item{\strong{\code{'eap'}}}{Expected-a-Posteriori Ability
#'       Estimation}
#'   }
#' @param ... Additional arguments passed to specific methods
#' @param prior_dist The shape of the prior distribution. Currently following
#'          distributions can be specified:
#'          \describe{
#'            \item{'norm'}{Normal distribution}
#'            \item{'unif'}{Uniform distribution}
#'            \item{'t'}{t distribution}
#'            \item{'cauchy'}{Cauchy distribution}
#'          }
#'          Default value is \code{'norm'}.
#' @param prior_pars Parameters of the prior distribution. Default value is
#'   \code{c(0, 1)} where 0 is the mean and 1 is the standard deviation of the
#'   default prior distribution which is normal distribution. Also, for example,
#'   uniform prior parameter can be set as \code{c(a, b)} where \code{a} is the
#'   minimum value and \code{b} is the maximum value. For \code{t} distribution,
#'   prior parameter can be set as \code{df} to represent the degree of freedom.
#'   For Cauchy distribution, prior parameters can be set as \code{c(location,
#'   scale)}.
#'
#'   If method is \code{"owen"}, provide \code{c(<Prior Mean>, <Prior SD>)}.
#'
#' @param theta_range The limits of the ability estimation scale. The estimation
#'   result will be limited to this interval. The default is \code{c(-5, 5)}.
#' @param number_of_quads Number of quadratures. The default value is 41. As
#'   this number increases, the precision of the estimate will also increase.
#'   The default value is \code{41}.
#' @param tol The precision level of ability estimate. The final ability
#'   estimates will be rounded to remove the precision that is smaller than the
#'   \code{tol} value. The default value is \code{1e-06}.
#'
#' @return \code{est} The ability estimated. If the response vector for a
#'   subject contains all \code{NA}s, then, in order to differentiate all
#'   incorrect and all NA, the \code{est} returned will be NA.
#' @return \code{se} The standard error(s) of the ability estimate(s). For
#'   \code{"sum_score"} method, all of the standard errors will be \code{NA}.
#'
#' @author Emre Gonulates
#' @export
#'
#' @references
#' Owen, R. J. (1975). A Bayesian sequential procedure for quantal response in
#' the context of adaptive mental testing. Journal of the American Statistical
#' Association, 70(350), 351-356.
#'
#' Vale, C. D., & Weiss, D. J. (1977). A Rapid Item-Search Procedure for
#' Bayesian Adaptive Testing. Research Report 77-4. Minneapolis, MN.
#'
#' @examples
#' ip <- generate_ip()
#' resp <- sim_resp(ip, theta = rnorm(1))
#'
#' # EAP estimation
#' est_ability(ip, resp)
#' est_ability(ip, resp, number_of_quads = 81)
#' est_ability(ip, resp, prior_pars = c(0, 3))
#' est_ability(ip, resp, prior_dist = 'unif',  prior_pars = c(-3, 3))
#' est_ability(ip, resp, prior_dist = 't',  prior_pars = 3)
#' est_ability(ip, resp, prior_dist = 'cauchy',  prior_pars = c(0, 1))
#'
#' # Maximum Likelihood estimation
#' est_ability(ip, resp, method = 'ml')
#' est_ability(ip, resp, method = 'ml', tol = 1e-8)
#' est_ability(ip, resp = rep(1, length(ip)), method = 'ml')
#' est_ability(ip, resp = rep(1, length(ip)), method = 'ml',
#'             theta_range = c(-3, 3))
#'
#' # Owen's Bayesian ability estimation
#' est_ability(ip, resp, method = 'owen')
#' est_ability(ip, resp, method = 'owen', prior_pars = c(0, 3))
#'
est_ability <- function(ip, resp, method = "eap", ...,
                        prior_dist = "norm",
                        prior_pars = c(0, 1),
                        theta_range = c(-5, 5),
                        number_of_quads = 41,
                        tol = 0.000001) {
  # args <- list(...)
  if (!is.matrix(resp)) {
    if (is.numeric(resp)) {
      if (length(resp) != length(ip))
        stop(paste0("The length of the response pattern should be equal ",
                    "to the number of items."))
      resp <- matrix(resp, nrow = 1)
    } else if (inherits(resp, "tbl_df")) { # if the response vector is a 'tibble'
      resp <- as.matrix(resp)
    } else if (all(is.na(resp))) { # if the response vector is all NAs
      return(list(est = NA, se = NA))
    } else stop("Invalid response pattern.")
  }
  if (!is(ip, "Itempool") | !is(ip, "Item")) ip <- itempool(ip)
  switch(
    method,
    "sum_score" = {
      se <- est <- rep(NA, nrow(resp))
      est <- rowSums(resp, na.rm = TRUE)
    },
    "owen" = {
      # This method cannot be used for models other than dichotomous IRT models.
      if (!all(ip$model %in% names(Pmodels)[sapply(Pmodels, function(x)
        x$model_family == "UIRT")]))
        stop(paste0("Owen's Bayesian ability estimation method can only ",
                    "be used for dichotomous IRT models: ",
                    paste0(names(Pmodels)[sapply(Pmodels, function(x)
                      x$model_family == "UIRT")], collapse = ", "), "."))
      output <- sapply(apply(resp, 1, function(r) est_ability_owen_cpp(
        ip = ip, resp = r, m0 = prior_pars[1], v0 = prior_pars[2]^2)
        ), function(x) x)
      est <- unlist(output[1, ])
      se <- unlist(output[2, ])
    },
    "ml" = {
      se <- est <- rep(NA, nrow(resp))
      est <- apply(resp, 1, function(r) stats::optimize(f = function(x)
        resp_loglik_itempool_cpp(resp = matrix(r, nrow = 1), theta = x, ip = ip),
        interval = theta_range, maximum = TRUE, tol = tol*.5)$maximum)
      se <- as.vector(1/sqrt(info(ip, theta = est, tif = TRUE, resp = resp)))
    },
    "eap" = {
      output <- est_ability_eap_cpp(
        resp = resp, ip = ip, theta_range = theta_range,
        no_of_quadrature = number_of_quads,
        prior_dist = prior_dist, prior_par = prior_pars)
      est <- output$est
      se <- output$se

    },
    stop("This method has not been implemented yet.")
    )

  # Round numbers up to tolerance
  est <- stats::setNames(round(est, floor(abs(log10(tol)))), rownames(resp))
  se <- stats::setNames(round(se, floor(abs(log10(tol)))), rownames(resp))
  # Convert all NA response strings's est and se to NA
  all_na_rows <- apply(is.na(resp), 1, all)
  est[all_na_rows] <- NA
  se[all_na_rows] <- NA
  return(list(est = est, se = se))
}

