fun_updatelambda <- function(b, ew, arglist) {

  l <- arglist$l
  r <- arglist$r
  u <- arglist$set[,2]
  z <- arglist$z
  trunc <-arglist$trunc

  exp_zb <- exp(z %*% b)

  r_star <- r
  r_star[is.infinite(r)] <- l[is.infinite(r)]

  # target_set <- fun_subless(u = u, lessthan = r_star)

  if (is.null(trunc)) {
    target_set <- fun_subless(u = u, lessthan = r_star)
  } else {
    target_set <- fun_sublr(u = u, l = trunc-1e-10, r = r_star)
  }

  numer <- colSums(target_set * ew)
  denom <- colSums(target_set * as.numeric(exp_zb))

  update_lambda <- numer/denom

  return(update_lambda)
}
