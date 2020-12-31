log_penlikelihood <- function(b, arglist) {

  n <- arglist$n
  npmle_set <- arglist$set
  u <- npmle_set[,2]
  l <- arglist$l
  r <- arglist$r
  r_cen <- is.infinite(r)
  lr <- cbind(l, r)
  z <- arglist$z
  trunc <- arglist$trunc
  tol <- arglist$tol
  niter <- arglist$niter

  distance <- tol + 1000
  iter <- 1
  old_lambda <- arglist$initial_lambda

  while ((distance > tol) & (iter < niter)) {
    ew <- fun_ew(b, old_lambda, arglist)
    new_lambda <- fun_updatelambda(b, ew, arglist)
    distance <- max(abs(new_lambda - old_lambda))
    old_lambda <- new_lambda

    iter <- iter + 1
  }

  exp_zb <- exp(z %*% b)
  lambda_exp_zb <- exp_zb %x% t(new_lambda)

  if (is.null(trunc)) {
    target_set1 <- fun_subless(u = u, lessthan = l)
  } else {
    target_set1 <- fun_sublr(u = u, l = trunc-1e-10, r = l)
  }

  target_set2 <- fun_sublr(u = u, l = l, r = r)

  value <- sum(-rowSums(target_set1 * lambda_exp_zb)) + sum(log((1 - exp(- rowSums(target_set2 * lambda_exp_zb)))[!r_cen]))

  return(value)

}
