#' Summary for \code{emfrail} objects
#'
#' @importFrom expint expint_En
#' @param object An object of class \code{emfrail}
#' @param lik_ci Logical. Should the confidence intervals for the frailty parameter be calculated based on the likelihood? If not, they are calculated with the delta method.
#' @param print_opts A list with options for printing the summary object. These include \code{coef}, \code{dist}, \code{fit}, \code{frailty}, \code{adj_se}, \code{verbose_frailty}.
#' @param ... Ignored
#'
#' @return An object of class \code{emfrail_summary},
#' with some more human-readable results from an \code{emfrail} object.
#'
#' @details
#' Regardless of
#' the fitted model, the following fields will be present in this object: \code{est_dist} (an object of class \code{emfrail_distribution}) with the estimated
#' distribution, \code{loglik} (a named vector with the log-likelihoods of the no-frailty model, the frailty model,
#' the likelihood ratio test statistic and the p-value of the one-sided likelihood ratio test), \code{theta} (a named vector with
#' the estimated value of the parameter \eqn{\theta}, the standard error, and the limits of a 95% cofindence interval) and \code{frail}, which
#' is a data frame with the following columns: \code{id} (cluster identifier), \code{z} (empirical Bayes frailty estimates), and optional
#' \code{lower_q} and \code{upper_q} as the 2.5% and 97.5% quantiles of the posterior distribution of the frailties (only for gamma distribution).
#'
#' For the the PVF or gamma distributions, the field \code{fr_var} contains a transformation of \code{theta} to correspond to the
#' frailty variance.
#' The fields \code{pvf_pars} and \code{stable_pars} are for quantities that are calculated only when the distribution is PVF or stable.
#' If the model contains covariates, the field \code{coefmat} contains the corresponding estimates. The p-values are based on
#' the adjusted standard errors, if they have been calculated successfully (i.e. if they appear when prining the summary object).
#' Otherwise, they are based on the regular standard errors.
#'
#' @export
#' @method summary emfrail
#'
#' @seealso \code{\link{predict.emfrail}, \link{plot.emfrail}}
#'
#' @examples
#' data("bladder")
#' mod_gamma <- emfrail(Surv(start, stop, status) ~ treatment + cluster(id), bladder1)
#' summary(mod_gamma)
#' summary(mod_gamma, print_opts = list(frailty_verbose = FALSE))
#'
#' # plot the Empirical Bayes estimates of the frailty
#' # easy way:
#' plot(mod_gamma, type = "hist")
#'
#' # a fancy graph:
#' sum_mod <- summary(mod_gamma)
#' library(dplyr)
#' library(ggplot2)
#'
#' # Create a plot just with the points
#' pl1 <- sum_mod$frail %>%
#'   arrange(z) %>%
#'   mutate(x = 1:n()) %>%
#'   ggplot(aes(x = x, y = z)) +
#'   geom_point()
#'
#' # If the quantiles of the posterior distribution are
#' # known, then error bars can be added:
#' if(!is.null(sum_mod$frail$lower_q))
#'   pl1 <- pl1 + geom_errorbar(aes(ymin = lower_q, ymax = upper_q), alpha = 0.5)
#'
#' pl1
#'
#' # The plot can be made interactive!
#' # ggplot2 gives a warning about the "id" aesthetic, just ignore it
#' pl2 <- sum_mod$frail %>%
#'   arrange(z) %>%
#'   mutate(x = 1:n()) %>%
#'   ggplot(aes(x = x, y = z)) +
#'   geom_point(aes(id = id))
#'
#' if(!is.null(sum_mod$z$lower_q))
#'   pl2 <- pl2 + geom_errorbar(aes(ymin = lower_q, ymax = upper_q, id = id), alpha = 0.5)
#'
#' library(plotly)
#' ggplotly(pl2)
#'
#' # Proportional hazards test
#' off_z <- log(sum_mod$frail$z)[match(bladder1$id, sum_mod$frail$id)]
#'
#' zph1 <- cox.zph(coxph(Surv(start, stop, status) ~ treatment + cluster(id), data = bladder1))
#'
#' # no sign of non-proportionality
#' zph2 <- cox.zph(coxph(Surv(start, stop, status) ~ treatment + offset(off_z), data = bladder1))
#'
#' zph2
#' # the p-values are even larger; the frailty "corrects" for proportionality.

summary.emfrail <- function(object,
                            lik_ci = TRUE,
                            print_opts = list(coef = TRUE,
                                              dist = TRUE,
                                              fit = TRUE,
                                              frailty = TRUE,
                                              adj_se = TRUE,
                                              verbose_frailty = TRUE),
                            ...) {

  # Calculate the following: estimated distribution of the frailty at time 0
  # fit <- object
  est_dist <- object$distribution
  est_dist$frailtypar <- exp(object$logtheta)

  loglik <- c(L0 = object$loglik[1],
              L = object$loglik[2],
              LRT = 2 * (object$loglik[2] - object$loglik[1]),
              pval = pchisq(2 * (object$loglik[2]  - object$loglik[1]), df = 1, lower.tail = FALSE)/2)

  # theta
  theta <- exp(object$logtheta)
  if(is.na(object$var_logtheta)) {
    se_theta <- NA
  } else {
    se_theta <- msm::deltamethod(~exp(x1),
                                 mean = object$logtheta,
                                 cov = object$var_logtheta)
  }


  if(!isTRUE(lik_ci) | !isTRUE(object$control$lik_ci)) {
    # CI symmetric on log(theta) scale,
    ci_theta_low <- exp(object$logtheta - 1.96 * sqrt(object$var_logtheta))
    ci_theta_high <-  exp(object$logtheta + 1.96 * sqrt(object$var_logtheta))

    # if theta was at the edge, then CI should show this....
    if(theta > object$control$em_control$upper_tol - 0.1) {
      ci_theta_low <- theta
      ci_theta_high <- Inf
    }
  } else {
    ci_theta_low <- exp(object$ci_logtheta[1])
    ci_theta_high <- exp(object$ci_logtheta[2])
  }

  # for gamma and pvf theta is 1/variance
  # for stable the L.T. is exp(- c^(1 - theta / (theta + 1)))

  gamma_pars <- stable_pars <- pvf_pars <- NULL
  fr_var <- se_fr_var <- ci_frvar_low <- ci_frvar_high <- NULL
  coefmat <- NULL


  # Frailty variance for gamma and pvf

  if(est_dist$dist != "stable") {
    fr_var <- 1/theta
    if(is.na(object$var_logtheta))
      se_fr_var <- NA else
        se_fr_var <- msm::deltamethod(~1/exp(x1),
                                      mean = object$logtheta,
                                      cov = object$var_logtheta)
    ci_frvar_low <- 1/ci_theta_high
    ci_frvar_high <- 1/ci_theta_low
  }


  # gamma-specific
  if (est_dist$dist == "gamma") {

    # Kendall's tau
    tau_gamma <- list(tau = 1 / (1 + 2 * theta),
                      se_tau =  msm::deltamethod(~1 / (1 + 2 * exp(x1)),
                                                      mean = object$logtheta,
                                                      cov = object$var_logtheta),
                      ci_tau_low = 1 / (1 + 2 * ci_theta_high),
                      ci_tau_high = 1 / (1 + 2 * ci_theta_low)
                      )

    # the CI seems to behave quite erratic...
    kappa_gamma <- list(kappa = 4 * (2^(1 + 1/theta) - 1)^(-theta) - 1,
                        se_kappa = msm::deltamethod(~4 * (2^(1 + 1/exp(x1)) - 1)^(-exp(x1)) - 1,
                                                         mean = object$logtheta,
                                                         cov = object$var_logtheta),
                        ci_kappa_low = ifelse(ci_theta_high == Inf,
                                              0,
                                              4 * (2^(1 + 1/ci_theta_high) - 1)^(-ci_theta_high) - 1),
                        ci_kappa_high = 4 * (2^(1 + 1/ci_theta_low) - 1)^(-ci_theta_low) - 1
                        )

    e_log_z_gamma <- list(e_log_z = digamma(theta) - log(theta),
                          se_e_log_z = msm::deltamethod(~digamma(exp(x1)) - log(exp(x1)) ,
                                                             mean = object$logtheta,
                                                             cov = object$var_logtheta),
                          ci_e_log_z_low = digamma(ci_theta_low) - log(ci_theta_low),
                          ci_e_log_z_high = ifelse(ci_theta_high == Inf,
                                                   0,
                                                   digamma(ci_theta_high) - log(ci_theta_high))
                          )

    var_log_z_gamma <- list(var_log_z = trigamma(theta),
                          se_var_log_z = msm::deltamethod(~trigamma(exp(x1)) ,
                                                             mean = object$logtheta,
                                                             cov = object$var_logtheta),
                          ci_var_log_z_low = trigamma(ci_theta_high),
                          ci_var_log_z_high = trigamma(ci_theta_low)
                          )




    shape <- est_dist$frailtypar + object$nevents_id
    rate <- est_dist$frailtypar + object$residuals[[1]]
    var_z <- shape / rate^2
    lower_q <- stats::qgamma(0.025,
                             shape = shape,
                             rate = rate)
    upper_q <- stats::qgamma(0.975,
                             shape = shape,
                             rate = rate)

    gamma_pars <- list(
      tau = tau_gamma,
      kappa = kappa_gamma,
      e_log_z = e_log_z_gamma,
      var_log_z = var_log_z_gamma,
      shape = shape,
      rate = rate,
      var_z = var_z,
      lower_q = lower_q,
      upper_q = upper_q
    )

  }

  if(est_dist$dist == "stable") {

    # Kendall's tau
    tau_stab <- list(
      tau = 1 - theta / (theta + 1),
      se_tau = msm::deltamethod(~ 1 - (exp(x1) / (exp(x1) + 1)),
                                           mean = object$logtheta,
                                           cov = object$var_logtheta),
      ci_tau_low = ifelse(ci_theta_high == Inf,
                                 0,
                                 1 - ci_theta_high / (1 + ci_theta_high)),
      ci_tau_high = 1 - ci_theta_low / (1 + ci_theta_low)
    )

    kappa_stab <- list(
      kappa = 2^(2 - 2^(theta / (theta + 1))) - 1,
      se_kappa = with(object$outer_m,
                    msm::deltamethod(~ 2^(2 - 2^(exp(x1) / (exp(x1) + 1))) - 1,
                                     mean = object$logtheta,
                                     cov = object$var_logtheta)),
      ci_kappa_low = ifelse(ci_theta_high == Inf,
                            0,
                            2^(2 - 2^(ci_theta_high / (ci_theta_high + 1))) - 1),
      ci_kappa_high = 2^(2 - 2^(ci_theta_low / (ci_theta_low + 1))) - 1
    )

    e_log_z_stab <- list(
      e_log_z = (-1) *  ((theta + 1) / theta - 1) * digamma(1),
      se_e_log_z = msm::deltamethod(~ (-1) *  ((exp(x1) + 1) / exp(x1) - 1) * digamma(1),
                                         mean = object$logtheta,
                                         cov = object$var_logtheta),
      ci_e_log_z_low = ifelse(ci_theta_high == Inf,
                              0,
                              (-1) *  ((ci_theta_high + 1) / ci_theta_high - 1) * digamma(1)),
      ci_e_log_z_high = (-1) *  ((ci_theta_low + 1) / ci_theta_low - 1) * digamma(1)
    )

    var_log_z_stab <- list(
      var_log_z = (((1 + theta) / theta)^2 - 1) * trigamma(1),
      se_var_log_z = msm::deltamethod(~ (((1 + exp(x1)) / exp(x1))^2 - 1) * trigamma(1),
                                           mean = object$logtheta,
                                           cov = object$var_logtheta),
      ci_var_log_z_low = ifelse(ci_theta_high == Inf,
                                0,
                                (((1 + ci_theta_high) / ci_theta_high)^2 - 1) * trigamma(1)),
      ci_var_log_z_high = (((1 + ci_theta_low) / ci_theta_low)^2 - 1) * trigamma(1)

    )


    attenuation <- list(
      attenuation = theta / (theta + 1),
      se_attenuation = msm::deltamethod(~ exp(x1) / (exp(x1) + 1),
                                           mean = object$logtheta,
                                           cov = object$var_logtheta),
      ci_attenuation_low = ci_theta_low / (ci_theta_low + 1),
      ci_attenuation_high = ifelse(ci_theta_high == Inf,
                                   1,
                                   ci_theta_high / (ci_theta_high + 1))
    )

    stable_pars <- list(
      tau = tau_stab,
      kappa = kappa_stab,
      e_log_z = e_log_z_stab,
      var_log_z = var_log_z_stab,
      attenuation = attenuation
    )
  }




  # dist_pars <- with(est_dist, dist_to_pars(dist, log(frailtypar), pvfm))
  if(est_dist$dist == "pvf") {

    if(est_dist$pvfm > 0) {
      pvfm <- est_dist$pvfm

      mass_at_0 <-  exp(-(pvfm+1)/pvfm * theta)
      pvf_pars <- list(mass_at_0 = mass_at_0)
    }

    if(est_dist$pvfm < 0) {
      pvfm <- est_dist$pvfm


      tau_pvf <- list(tau = suppressWarnings((1 + pvfm) - 2 * (pvfm + 1) * theta +
        4 * (pvfm + 1)^2 * theta^2 / (- pvfm) * exp(2 * (pvfm + 1) * theta / (- pvfm)) *
        expint::expint_En(2 * (pvfm + 1) * theta / (-pvfm), order = 1 / (- pvfm) - 1)),

         # se_tau = with(object$outer_m,
         #               msm::deltamethod(~ (1 + pvfm) - 2 * (pvfm + 1) * exp(x1) +
         #                                  4 * (pvfm + 1)^2 * exp(x1)^2 / (- pvfm) * exp(2 * (pvfm + 1) * exp(x1) / (- pvfm)) *
         #                                  expint::expint_En(2 * (pvfm + 1) * exp(x1) / (-pvfm), order = 1 / (- pvfm) - 1),
         #                                mean = minimum,
         #                                cov = 1/hess))

        se_tau = NULL,
        ci_tau_low = ifelse(ci_theta_high == Inf,
                            0,
                            (1 + pvfm) - 2 * (pvfm + 1) * ci_theta_high +
          4 * (pvfm + 1)^2 * ci_theta_high^2 / (- pvfm) * exp(2 * (pvfm + 1) * ci_theta_high / (- pvfm)) *
          expint::expint_En(2 * (pvfm + 1) * ci_theta_high / (-pvfm), order = 1 / (- pvfm) - 1)),
        ci_tau_high = (1 + pvfm) - 2 * (pvfm + 1) * ci_theta_low +
          4 * (pvfm + 1)^2 * ci_theta_low^2 / (- pvfm) * exp(2 * (pvfm + 1) * ci_theta_low / (- pvfm)) *
          expint::expint_En(2 * (pvfm + 1) * ci_theta_low / (-pvfm), order = 1 / (- pvfm) - 1)

      )

      # this may happen in the limiting case of theta very large (i.e. no dependence)
      if(is.nan(tau_pvf$tau)) tau_pvf$tau <- 0


      kappa_pvf <- list(kappa = 4 * exp(
        (-1) * (2 * ((pvfm + 1) * theta / (- pvfm) + log(2))^(-1/pvfm) -
                  ((pvfm + 1) * theta / (-pvfm) )^(- 1 / pvfm) )^(-pvfm) +
          (pvfm + 1) * theta / (-pvfm)
      ) - 1,
      # se_kappa = msm::deltamethod(~  4 * exp(
      #                 (-1) * (2 * ((pvfm + 1) * exp(x1) / (- pvfm) + log(2))^(-1/pvfm) -
      #                           ((pvfm + 1) * exp(x1) / (-pvfm) )^(- 1 / pvfm) )^(-pvfm) +
      #                   (pvfm + 1) * exp(x1) / (-pvfm)
      #               ) - 1,
      #                                mean = object$logtheta,
      #                                cov = object$var_logtheta),
      se_kappa = NULL,
      ci_kappa_low = ifelse(ci_theta_high == Inf,
                            0,
                            4 * exp(
                              (-1) * (2 * ((pvfm + 1) * ci_theta_high / (- pvfm) + log(2))^(-1/pvfm) -
                                        ((pvfm + 1) * ci_theta_high / (-pvfm) )^(- 1 / pvfm) )^(-pvfm) +
                                (pvfm + 1) * ci_theta_high / (-pvfm)) - 1),
      ci_kappa_high =    4 * exp(
        (-1) * (2 * ((pvfm + 1) * ci_theta_low / (- pvfm) + log(2))^(-1/pvfm) -
                  ((pvfm + 1) * ci_theta_low / (-pvfm) )^(- 1 / pvfm) )^(-pvfm) +
          (pvfm + 1) * ci_theta_low / (-pvfm)) - 1)

      if(est_dist$pvfm != -1/2)
        pvf_pars <- list(tau = tau_pvf,
                         kappa = kappa_pvf)

      if(est_dist$pvfm == -1/2) {
        # e_log_z_pvf = log(1/2 * theta) + expint::expint_E1(theta) * exp(theta)
        pvf_pars <- list(tau = tau_pvf,
                         kappa = kappa_pvf)
      }
    }


  }


  # The frailty estimates

  z <- data.frame(id = names(object$frail),
                  z = object$frail)

  if(!is.null(gamma_pars)) {
    z$lower_q <- as.numeric(gamma_pars$lower_q)
    z$upper_q <- as.numeric(gamma_pars$upper_q)
  }


  if(length(object$coef) > 0) {
    coefmat <- list(
      coef = object$coef,
      "exp(coef)" = exp(object$coef),
      "se(coef)" = sqrt(diag(object$var)[seq_along(object$coef)]),
      "adj. se" = sqrt(diag(object$var_adj)[seq_along(object$coef)] ))

      if(all(is.na(coefmat$`adj. se`))) {
        coefmat$`adj. se` <- NULL
        coefmat$z <- coefmat$coef / coefmat$`se(coef)`
      } else {
        coefmat$z <- coefmat$coef / coefmat$`adj. se`
      }

    coefmat$p <-  1 - pchisq(coefmat$z^2, df = 1)

    coefmat <- do.call(cbind, coefmat)
  }




  ret <- list(est_dist = est_dist,
       loglik = loglik,
       ca_test = object$ca_test,
       cens_test = object$cens_test,
       lik_ci = lik_ci,
       theta = c(theta = theta,
                 se_theta = se_theta,
                 ci_theta_low = ci_theta_low,
                 ci_theta_high = ci_theta_high),
       fr_var = c(fr_var = fr_var,
                  se_fr_var = se_fr_var,
                  ci_frvar_low = ci_frvar_low,
                  ci_frvar_high = ci_frvar_high),
       gamma_pars = gamma_pars,
       pvf_pars = pvf_pars,
       stable_pars = stable_pars,
       coefmat = coefmat,
       frail = z
       )


  #attr(ret, "class") <- "emfrail_summary"
  attr(ret, "print_opts") <-  print_opts
  attr(ret, "call") <- attr(object, "call")
  class(ret) <- "emfrail_summary"
  ret
}
