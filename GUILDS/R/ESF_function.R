logLikelihood.ESF <- function(theta, m, abund) {
  if (theta < 1 ||
     m > (1 - .Machine$double.eps) ||
     m <= 0) {
     return(-Inf)
  }

  J <- sum(abund)
  S <- length(abund)
  I <- m * (J-1) / (1 - m)

  kda <- calcKDA(abund)  #confirmed in PARI

  sumkda <- calc_sum_kda(S, J, I, theta, kda)

  x <- c(table(abund))
  freq_x <- c()
  for (i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor1 <- -1 * ( sum(log(abund)) + sum(lgamma(1 + freq_x)))

  #J!/[prod(n1)prod(Sx!)]  #confirmed in PARI
  factor1 <- lgamma(J+1) + prefactor1

  factor2 <- S*log(theta)   -  (lgamma(I + J) - lgamma(I))

  ll <- factor1 + factor2 + sumkda
  return(ll)
}

esf_local <- function(v, abund, prefactor, kda) {
  theta <- v[1]
  m <- v[2]

  if (theta < 1 ||
     m <= 0 ||
     m > (1 - .Machine$double.eps)) {
    return(-Inf)
  }

  J <- sum(abund)
  S <- length(abund)
  I <- m * (J - 1) / (1 - m)

  sumkda <- calc_sum_kda(S, J, I, theta, kda)

  factor2 <- S * log(theta)   -  (lgamma(I + J) - lgamma(I))

  ll <- prefactor + factor2 + sumkda
  return(ll)
}

maxLikelihood.ESF <- function(init_vals, abund, verbose = TRUE) {
  if (init_vals[1] < 1) {
     stop("maxLikelihood.ESF: ",
          "initial theta can not be below one")
   }
  if (init_vals[2] < 0) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be below zero")
  }
  if (init_vals[2] > 1) {
    stop("maxLikelihood.ESF: ",
         "initial m can not be above 1 (did you mean to enter I?)")
  }
  if (length(abund) < 2) {
    stop("maxLikelihood.ESF: ",
         "Need more than 1 species in the dataset")
  }

  kda <- calcKDA(abund)  #confirmed in PARI

  J <- sum(abund)
  x <- c( table(abund))
  freq_x <- c()
  for (i in seq_along(x)) freq_x[i] <- x[[i]]
  prefactor <- lgamma(J+1) - ( sum(log(abund)) + sum(lgamma(1 + freq_x)) )

  g <- function(x) {
	  out <- -1 * esf_local(x, abund, prefactor, kda)
		return(out)
  }

  optimum <- simplex(init_vals, g, verbose)
  return(optimum)
}