loglik.M12 <-
function (x, y, tol=0.000001) { # x males y females
  
  results <- fitmodel12(x,y,tol=tol)
  
  pa <- results$p
  f <- results$fhat
  
  alles <- c(x,y)
  paaf <- pa^2 + pa * (1 - pa) * f
  pabf <- 2 * pa * (1 - pa) * (1 - f)
  pbbf <- (1 - pa)^2 + pa * (1 - pa) * f
  paam <- pa^2
  pabm <- 2 * pa * (1 - pa)
  pbbm <- (1 - pa)^2
  
  pvec <- c(paam, pabm, pbbm, paaf, pabf, pbbf)
  ind <- !(alles == 0)
  logvec <- log(pvec[ind])
  loglik0 <- sum(alles[ind] * logvec)
  nparam0 <- 2
  res <- c(loglik0, nparam0)
  
  return(res)
}
