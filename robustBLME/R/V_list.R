##' @title V matrix and related
##' @description Computes the V matrix and related quantities
##' @usage V_list(sig2_b, sig2_eps, ZZt_b, ZZt_b_ii, ZZt_eps, ZZt_eps_ii, m_i)
##'
##' @param sig2_b Variance of the random effects.
##' @param sig2_eps  Variance of the error term.
##' @param ZZt_b The ZZt matrix for the random effects.
##' @param ZZt_b_ii The ii block of ZZt corresponding to the i-th unit.
##' @param ZZt_eps The ZZt matrix for the error term.
##' @param ZZt_eps_ii The ii block of ZZt_eps corresponding to the i-th unit.
##' @param m_i Number of replications.
##' @noRd

V_list <- function(sig2_b,
                   sig2_eps,
                   ZZt_b,
                   ZZt_b_ii,
                   ZZt_eps,
                   ZZt_eps_ii,
                   m_i) {
  .Call(
    'robustBLME_V_list',
    PACKAGE = 'robustBLME',
    sig2_b,
    sig2_eps,
    ZZt_b,
    ZZt_b_ii,
    ZZt_eps,
    ZZt_eps_ii,
    m_i
  )
}
