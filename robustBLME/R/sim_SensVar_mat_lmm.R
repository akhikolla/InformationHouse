##' @title Monte Carlo approximation of the sensitivity and sariability matrices
##'
##' @description Computes the sensitivity and the variability matrices for the REML2 estimating function via Monte Carlo. The derivative of the estimating function is done numerically.
##' @usage sim_SensVar_mat_lmm(n_sim, betas, sig2_b, sig2_eps, Xn, Xnt, ZZt_b, ZZt_b_ii,
##'                 ZZt_eps, ZZt_eps_ii, c_hub, c2_hub, K2n, m_i, n, n_ind)
##'
##' @param n_sim Number of Monte Carlo draws from the model.
##' @param betas Vector of fixed effects.
##' @param sig2_b Variance of the random effects.
##' @param sig2_eps Error variance.
##' @param Xn Model design matrix for all units.
##' @param Xnt The transposed version of Xn.
##' @param ZZt_b The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the random effects.
##' @param ZZt_b_ii The block of \code{ZZt_b} corresponding to the i-th individual.
##' @param ZZt_eps The matrix \eqn{$Z \times Z^t$}{Z times Z^t} for the error term.
##' @param ZZt_eps_ii The block of \code{ZZt_eps} corresponding to the i-th individual.
##' @param c_hub Tuning constant of Huber's psi function for location.
##' @param c2_hub Tuning constant of Huber's psi function for scale.
##' @param K2n The K2 matrix for consitency correction of scale parameters.
##' @param m_i The number of fixed effects.
##' @param n Overall sample size.
##' @param n_ind Number of statistical units.
##' @noRd
sim_SensVar_mat_lmm <-
  function(n_sim, betas, sig2_b, sig2_eps, Xn, Xnt, ZZt_b, ZZt_b_ii, ZZt_eps,
           ZZt_eps_ii, c_hub, c2_hub, K2n, m_i, n, n_ind)
  {
    cat('\n#########################################################\n Sensitivity & Variability mat.s with', n_sim, 'Monte Carlo samples\n')
    n_betas <- length(betas)

    blocki <- function(x, y_sim) {
      .Call('robustBLME_Psi_reml2',
            PACKAGE = 'robustBLME',
            x[1:n_betas],
            x[n_betas+1],
            x[n_betas+2],
            y_sim,
            Xn,
            Xnt,
            ZZt_b,
            ZZt_b_ii,
            ZZt_eps,
            ZZt_eps_ii,
            c_hub,
            c2_hub,
            K2n,
            n_betas,
            m_i,
            n,
            n_ind)
    }

    J <- H <- matrix(0, n_betas + 2, n_betas + 2)
    progress_bar <- txtProgressBar(min = 0, max = n_sim, style = 3)

    for (s in 1:n_sim) {
      y_sim <- simData_rlmm(betas = betas,
                            sig2_b = sig2_b,
                            sig2_eps = sig2_eps,
                            Xn = Xn,
                            ZZt_b_ii = ZZt_b_ii,
                            ZZt_eps_ii = ZZt_eps_ii,
                            m_i = m_i,
                            n_ind = n_ind)
      jj <- .Call('robustBLME_Psi_reml2',
                  PACKAGE = 'robustBLME',
                  betas,
                  sig2_b,
                  sig2_eps,
                  y_sim,
                  Xn,
                  Xnt,
                  ZZt_b,
                  ZZt_b_ii,
                  ZZt_eps,
                  ZZt_eps_ii,
                  c_hub,
                  c2_hub,
                  K2n,
                  n_betas,
                  m_i,
                  n,
                  n_ind)
      hh <- jacobian(blocki,
                     x = c(betas, sig2_b, sig2_eps),
                     y_sim = y_sim)

      H = H + hh
      J = J + tcrossprod(jj)

      # if(s%%floor(n_sim/2)==0) setTxtProgressBar(progress_bar, s)
      if(s%%10==0) setTxtProgressBar(progress_bar, s)
    }
    return(list(H = -H/n_sim, J = J/n_sim))
  }
