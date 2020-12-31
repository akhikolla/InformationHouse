loglik.M22 <-
function (pam, paf, x, y, fhatf) { # x males y females
  
  alles <- c(x,y)
  paaf <- paf^2 + paf * (1 - paf) * fhatf
  pabf <- 2 * paf * (1 - paf) * (1 - fhatf)
  pbbf <- (1 - paf)^2 + paf * (1 - paf) * fhatf
  paam <- pam^2
  pabm <- 2 * pam * (1 - pam)
  pbbm <- (1 - pam)^2
  
  pvec <- c(paam, pabm, pbbm, paaf, pabf, pbbf)
  ind <- !(alles == 0)
  logvec <- log(pvec[ind])
  loglik0 <- sum(alles[ind] * logvec)
  nparam0 <- 3
  res <- c(loglik0, nparam0)
  return(res)
}
