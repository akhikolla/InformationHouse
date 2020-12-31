fun_shooting_algorithm <- function(x, y, theta, tilde_b, arglist) {

  xx <- t(x)%*%x
  xy <- t(x)%*%y

  z <- arglist$z
  num_p <- ncol(z)
  tol <- arglist$tol
  niter <- arglist$niter

  old_b <- tilde_b
  temp_b <- old_b

  weighted_theta <- theta/abs(tilde_b)

  distance <- tol + 1000
  iter <- 1

  while ((distance > tol) & (iter < niter)) {

    for (p in 1:num_p) {
      ss2  <- 2*(sum(temp_b*xx[,p]) - temp_b[p]*xx[p,p] - xy[p])

      if (ss2 > weighted_theta[p]) temp_b[p] <- (weighted_theta[p] - ss2)/(2*xx[p,p])
      if (ss2 < -weighted_theta[p]) temp_b[p] <- (-weighted_theta[p] - ss2)/(2*xx[p,p])
      if (abs(ss2) <= weighted_theta[p]) temp_b[p] <- 0
    }

    distance <- max(abs(temp_b - old_b))

    old_b <- temp_b
    iter <- iter + 1

  }

  return(temp_b)



}
