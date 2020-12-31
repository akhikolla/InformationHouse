fun_ew <- function(b, lambda, arglist) {

  l <- arglist$l
  r <- arglist$r
  lr <- cbind(l, r)
  u <- arglist$set[,2]
  z <- arglist$z

  exp_zb <- exp(z %*% b)
  lambda_exp_zb <- exp_zb %*% t(lambda)
  target_set <- fun_sublr(u = u, l = l, r = r)

  denom <- 1 - exp(-rowSums(lambda_exp_zb * target_set))
  ew <- target_set * lambda_exp_zb/denom
  ew[is.infinite(r), ] <- 0

  return(ew)

}
