fun_bic <- function(x, tilde_b, arglist, tol) {

  if(missing(tol)) tol <- arglist$tol

  can_results <- fun_penSurvIC(theta = x, tilde_b, arglist, tol)
  can_b <- can_results$b
  can_convergence <- can_results$convergence
  log_pen <- log_penlikelihood(can_b, arglist)
  n <- arglist$n
  bic_value <- -2 * log_pen + log(n) * sum(can_b != 0)

  return(c(bic_value, can_convergence, can_b))
}
