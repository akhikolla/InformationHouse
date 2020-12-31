loglik.M23 <-
function (pam, paf, x, y, fhatm) { # x males y females
  
  alles <- c(x,y)
  paaf <- paf^2
  pabf <- 2 * paf * (1 - paf)
  pbbf <- (1 - paf)^2
  paam <- pam^2 + pam * (1 - pam) * fhatm
  pabm <- 2 * pam * (1 - pam) * (1 - fhatm)
  pbbm <- (1 - pam)^2 + pam * (1 - pam) * fhatm
  
  pvec <- c(paam, pabm, pbbm, paaf, pabf, pbbf)
  ind <- !(alles == 0)
  logvec <- log(pvec[ind])
  loglik0 <- sum(alles[ind] * logvec)
  nparam0 <- 3
  res <- c(loglik0, nparam0)
  return(res)
}
