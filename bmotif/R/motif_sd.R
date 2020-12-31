motif_sd <- function(W, mc = NA, meanw = NA) {
  # the input is the weighted adjacency matrix
  # output is a vector with standard deviations of the motif weight mean
  # v[i] gives sd for motif i
  # we do the main computation in C++, so use package Rcpp

  # library(Rcpp)

  W <- as.matrix(W)
  M <- W
  M[M > 0] <- 1

  if (any(is.na(mc))) {
    # user did not give any input for this, so we compute it again
    mc <- mcount (M, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  }

  if (any(is.na(meanw))) {
    # user did not give any input for this, so we compute it again
    meanw <- mean_weight(W, mc = mc, six_node = FALSE)
  }

  NZ <- nrow(W)
  NP <- ncol(W)

  N <- 1 - M

  DZ <- matrix(rep(1, NZ * NZ), nrow = NZ)
  for (i in 1:NZ) {
    DZ[i,i] <- 0
  }
  DP <- matrix(rep(1, NP * NP), nrow = NP)
  for (i in 1:NP) {
    DP[i,i] <- 0
  }

  # sourceCpp('weights/cppcode/sd.cpp')

  msd <- rep(NA, 17)

  # for motif 1, need standard deviation over all non-zero entries in W
  v <- as.vector(W)
  v <- v[v > 0]
  msd[1] <- pop_sd(v)

  # for motifs 2-17, calculation is done in C++
  for (i in 2:17) {
    num <- get(paste('sd_m', i, sep = ""))(NZ, NP, W, meanw[i])
    msd[i] <- sqrt(num / mc[i])
  }

  # we now could have NaN values in the msd vector, in case some motifs don't occur
  msd <- replace(msd, which(is.nan(msd)),NA)
  msd
}
