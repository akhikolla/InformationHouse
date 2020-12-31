##' @title Simulate data from an LMM
##' @description Simulates data from the LMM with one random effects
##' @usage simData_rlmm(betas, sig2_b, sig2_eps, Xn, ZZt_b_ii, ZZt_eps_ii, m_i, n_ind)

##' @param betas Vector of fixed effects.
##' @param sig2_b Variance of the random effects.
##' @param sig2_eps Variance of the error term.
##' @param Xn The design matrix for all statistical units.
##' @param ZZt_b_ii The ii block of the ZZt matrix of the random effects.
##' @param ZZt_eps_ii The ii block of the ZZt matrix of the error term.
##' @param m_i Number of replications.
##' @param n_ind Number of statistical units.
##' @noRd
simData_rlmm <- function(betas, sig2_b, sig2_eps, Xn, ZZt_b_ii, ZZt_eps_ii, m_i, n_ind) {
  .Call('robustBLME_simData_rlmm',
        PACKAGE = 'robustBLME',
        betas,
        sig2_b,
        sig2_eps,
        Xn,
        ZZt_b_ii,
        ZZt_eps_ii,
        m_i,
        n_ind)
}


##' @title Simulate data from a multivariate mixture of Gaussians
##' @description Simulates data from an LME
##' @usage simData_mixture_rlmm(delta, epsilon, betas, sig2_b,
##'               sig2_eps, X, ZZt_b, ZZt_b_ii, ZZt_eps, ZZt_eps_ii,
##'               m_i, n_ind)
##'
##' @param delta The contamination weight.
##' @param epsilon The variance inflation.
##' @param betas Vector of fixed effects.
##' @param sig2_b Variance of the random effects.
##' @param sig2_eps sig2_eps Variance of the error term.
##' @param X The design matrix for all statistical units.
##' @param ZZt_b The ZZt matrix for the random effects.
##' @param ZZt_b_ii The ii block of the ZZt matrix of the random effects.
##' @param ZZt_eps The ZZt matrix for the error term.
##' @param ZZt_eps_ii The ii block of the ZZt matrix of the error term.
##' @param m_i Number of replications.
##' @param n_ind Number of statistical units.
##' @noRd

simData_mixture_rlmm <- function(delta,
                                 epsilon,
                                 betas,
                                 sig2_b,
                                 sig2_eps,
                                 X,
                                 ZZt_b,
                                 ZZt_b_ii,
                                 ZZt_eps,
                                 ZZt_eps_ii,
                                 m_i,
                                 n_ind) {
  Vstuff <- V_list(sig2_b,
                   sig2_eps,
                   ZZt_b,
                   ZZt_b_ii,
                   ZZt_eps,
                   ZZt_eps_ii,
                   m_i)
  S1 <- Vstuff$V_i
  S2 <- S1
  S2 <- epsilon*S1
  # S = list(S1, S2)

  mu <- X%*%betas
  out <- matrix(NA, nrow = m_i, ncol = n_ind)

  for(i in 1:n_ind){
    w <- sample(c(0,1), size = 1, prob = c(1-delta, delta), replace = TRUE)
    out[,i] <- rmvnorm(mu, (1-w)*S1 + w*S2, m_i)
  }
  return(out)
}
