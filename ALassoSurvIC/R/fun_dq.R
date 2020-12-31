fun_dq <- function(b, ew, arglist) {

  l <- arglist$l
  r <- arglist$r
  lr  <- cbind(l, r)
  u <- arglist$set[,2]
  n <- arglist$n
  z <- arglist$z
  trunc <-arglist$trunc

  zb <- z %*% b
  exp_zb <- exp(z %*% b)

  r_star <- r
  r_star[is.infinite(r)] <- l[is.infinite(r)]

  # target_set <- fun_subless(u = u, lessthan = r_star)

  if (is.null(trunc)) {
    target_set <- fun_subless(u = u, lessthan = r_star)
  } else {
    target_set <- fun_sublr(u = u, l = trunc-1e-10, r = r_star)
  }

  zero_part <- target_set * ew
  first_part <- t(target_set * as.numeric(exp_zb)) %*% z
  second_part <- colSums(target_set * as.numeric(exp_zb))

  value <- colSums(zero_part %*% (- first_part/second_part) + rowSums(zero_part)*z)

  return(value)

}
