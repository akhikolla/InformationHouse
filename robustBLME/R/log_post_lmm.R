##' @title Posterior Distribution for Gaussian LMM with one random effects
##' @description Computes the unormalised log-posterior distribution for the Gaussian LMM with one random effect assuming Gaussian priors for the fiexed effects and halfCuchy priors for the variance components.
##' @usage log_post_lmm(param, y, Xn, ZZt_b_ii, ZZt_eps_ii,
##'              prior_beta_stdev, prior_sig2_scale, n_betas,
##'               m_i, n, n_ind)
##'
##' @param param The parameter vector made by the fixed effects, the log of the variance of random effects and log of error variance.
##' @param y The response variable. Must be a matrix.
##' @param Xn The design matrix for the all statistical units.
##' @param ZZt_b_ii The ii block of the ZZt_b matrix.
##' @param ZZt_eps_ii The ii block of the ZZt_eps matrix.
##' @param prior_beta_stdev Standard deviation of the Gaussian prior for the fixed effects. The prior mean is assumed equal to zero.
##' @param prior_sig2_scale Scale of the halfCauchy prior for the variance components. Here both priors have the same scale.
##' @param n_betas Number of fixed effects.
##' @param m_i Number of replications for each individual
##' @param n Overall number of observations.
##' @param n_ind Similar to m_i.
##' @noRd
log_post_lmm <- function(param, y, Xn, ZZt_b_ii, ZZt_eps_ii,
                         prior_beta_stdev, prior_sig2_scale,
                         n_betas, m_i, n, n_ind){
  ll <- .Call('robustBLME_log_lik_lmm', PACKAGE = 'robustBLME',
              param[1:n_betas], exp(param[n_betas+1]), exp(param[n_betas+2]), y,
              Xn, ZZt_b_ii, ZZt_eps_ii, n_betas,
              m_i, n, n_ind);
  pr <- dPrior_lmm(param[1:n_betas], param[n_betas+1], param[n_betas+2],
                   prior_beta_stdev, prior_sig2_scale, n_betas, TRUE)
  return(ll + pr)
}
