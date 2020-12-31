calc_sum_kda <- function(S, N, I, theta, kda) {
  numbers <- S:N
  s <- length(I)
  results <- matrix(0, nrow = s, ncol = 1)

  for (i in seq_len(s)) {
    outcomes <- kda[numbers + 1] + (log(I[i]) * numbers)  -
      (lgamma(theta + numbers) - lgamma(theta))

    log_max <- max(outcomes)
    sumlist <- sum(exp(outcomes - log_max))
    results[i] <- log_max + log(sumlist)
  }

  return(results)
}