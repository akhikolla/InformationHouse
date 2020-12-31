fun_ddq <- function(b, ew, arglist) {

  l <- arglist$l
  r <- arglist$r
  lr  <- cbind(l, r)
  u <- arglist$set[,2]
  n <- arglist$n
  z <- arglist$z
  len_z <- ncol(z)

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
  first_part1 <- t(fun_hcross(z) %*% (target_set * as.numeric(exp_zb)))
  first_part2 <- fun_hcross(t(target_set * as.numeric(exp_zb)) %*% z)
  second_part <- colSums(target_set * as.numeric(exp_zb))

  cvalue_temp1 <- - first_part1/second_part + t(first_part2)/(second_part^2)
  cvalue <- zero_part %*% cvalue_temp1
  value <- matrix(0, len_z, len_z)
  value[upper.tri(value, diag = TRUE)] <- colSums(cvalue)
  value[lower.tri(value)] <- t(value)[lower.tri(value)]

  return(value)

}
