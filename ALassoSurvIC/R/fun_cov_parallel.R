fun_cov_parallel <- function(b, theta, var.h, arglist, cl) {

  n <- arglist$n

  length_b <- length(b)
  profile_hess <- matrix(NA, length_b, length_b)

  bone <- matrix(NA, nrow = length_b*(length_b+1)/2, ncol = 2)
  com_ele <- 1

  for ( ej in 1:length_b) {
    for (ek in 1:ej) {
      bone[com_ele,] <- c(ej, ek)
      com_ele <- com_ele + 1
    }
  }

  h <- var.h * n^(-1/2)

  if (is.null(cl)) {
    vbone <- t(apply(bone, 1, fun_covij, b = b, length_b = length_b, h = h, arglist = arglist))
  } else {
    lbone <- split(bone, row(bone))
    if (inherits(cl, "cluster")) {
      parallel_fun <- if (isTRUE(getOption("pboptions")$use_lb)) parLapplyLB else parLapply
      vbone0 <- parallel_fun(cl, lbone, fun_covij, b = b, length_b = length_b, h = h, arglist = arglist)
    } else {
      vbone0 <- mclapply(lbone, fun_covij, b = b, length_b = length_b, h = h, arglist = arglist)
    }
    vbone <- t(sapply(vbone0, function(x) x))
  }

  for (i in 1:nrow(bone)) {
    profile_hess[vbone[i,1], vbone[i,2]] <- vbone[i,3]
  }

  profile_hess[upper.tri(profile_hess)] <- t(profile_hess)[upper.tri(profile_hess)]

  #variance computation
  vbeta <- abeta <- numeric(length_b)

  abeta[abs(b) > 0] <- 1/abs((b[abs(b) > 0])^2)
  abeta[b == 0] <- 10e10

  part1 <- profile_hess + diag(n*theta*abeta)
  inv_part1 <- solve(part1)

  vbeta[abs(b) > 0] <- 1/abs((b[abs(b) > 0])^2)
  vbeta[b == 0] = 0.0
  part2 <- profile_hess + diag(n*theta*vbeta)

  inv_hess <- solve(profile_hess)

  cov <- inv_part1 %*% part2 %*% inv_hess %*% part2 %*% inv_part1

  return(cov)

}
