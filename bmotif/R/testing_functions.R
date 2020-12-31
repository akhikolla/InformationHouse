block_matrix <- function(mlist) {
  # input: list of matrices
  # output: one big block matrix with those matrices on the diagonal

  N <- sum(sapply(mlist, nrow)) # total number of rows
  M <- sum(sapply(mlist, ncol)) # total number of columns

  len <- length(mlist)
  bm <- matrix(rep(0,N*M), nrow = N)

  r <- 1
  c <- 1

  for (i in 1:len) {
    A <- mlist[[i]]
    n <- nrow(A)
    m <- ncol(A)
    bm[r: (r+n-1), c: (c+m-1)] <- A
    r <- r+n
    c <- c+m
  }
  bm
}

rbm <- function (m, n, one_prob = 0.5) {
  matrix(sample(0:1, m * n, replace = TRUE, prob = c(1- one_prob,one_prob)), m, n)
}

rwm <- function(n, m, edge_prob = 0.5, min = 0, max = 1) {
  # first do a random binary matrix
  bin <- sample(0:1, m * n, replace = TRUE, prob = c(1- edge_prob,edge_prob))
  ind <- which(bin == 1)
  # now assign random weights to the edges present in the binary matrix
  weights <- stats::runif(length(ind), min, max)
  bin[ind] <- weights
  M <- matrix(bin, n, m)
  M
}

motif_sd_minors <- function(W) {
  # meanw is a vector of mean weights
  # mc is a vector of motif counts
  M <- W
  M[M > 0] <- 1
  mcdf <- mcount(W, six_node = FALSE, mean_weight = TRUE, standard_dev = FALSE, normalisation = FALSE)
  meanw <- mcdf$mean_weight
  mc <- mcdf$frequency

  # compute all_mot_perm
  all_mot_perm <- gen_all_mot_perm()
  m <- nrow(W)
  n <- ncol(W)

  sd_sum <- rep(0, 17)

  for (i in 1:m) {
    # i is now nrow of the minor
    rp <- gtools::combinations(m, i, 1:m, repeats = FALSE)
    for (j in 1:n) {
      # j is now ncol of the minor
      if (i + j > 5) {
        next
      }
      cp <- gtools::combinations(n, j, 1:n, repeats = FALSE)
      nrp <- nrow(rp)
      ncp <- nrow(cp)
      for (k in 1:nrp) {
        for (l in 1:ncp) {
          minor <- W[rp[k,], cp[l,], drop = FALSE]
          bin_minor <- minor
          bin_minor[bin_minor > 0] <- 1

          mt <- check_motif(bin_minor, all_mot_perm)

          if (mt == -1) {
            next
          }

          minor_nozero <- as.vector(minor[minor != 0])
          # compute numerator for the sd
          d <- mean(minor_nozero)
          sd_sum[mt] <- sd_sum[mt] + (d - meanw[mt])^2
        }
      }
    }
  }

  msd <- sqrt(sd_sum / mc)
  msd <- replace(msd, which(is.nan(msd)), NA)
  msd
}

mean_weight_minors <- function(W, six_node = FALSE) {
  # mc is a vector of motif counts
  M <- W
  M[M > 0] <- 1
  mc <- mcount(M, six_node = six_node, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency

  ###### compute all_mot_perm
  all_mot_perm <- gen_all_mot_perm()
  m <- nrow(W)
  n <- ncol(W)

  if (six_node) {
    sum_weights <- rep(0, 44)
  }
  if (!six_node) {
    sum_weights <- rep(0, 17)
  }


  for (i in 1:m) {
    # i is now nrow of the minor
    rp <- gtools::combinations(m, i, 1:m, repeats = FALSE)
    for (j in 1:n) {
      # j is now ncol of the minor
      if (i + j > (5 + six_node)) {
        next
      }
      cp <- gtools::combinations(n, j, 1:n, repeats = FALSE)
      nrp <- nrow(rp)
      ncp <- nrow(cp)
      for (k in 1:nrp) {
        for (l in 1:ncp) {
          minor <- W[rp[k,], cp[l,], drop = FALSE]
          bin_minor <- minor
          bin_minor[bin_minor > 0] <- 1

          mt <- check_motif(bin_minor, all_mot_perm)
          if (mt == -1) {
            next
          }
          minor_nozero <- as.vector(minor[minor != 0])
          # compute numerator for the sd
          d <- mean(minor_nozero)
          sum_weights[mt] <- sum_weights[mt] + d

        }
      }
    }
  }

  mw <- sum_weights / mc
  mw <- replace(mw, which(is.nan(mw)), NA)
  mw
}

mcount_minors <- function(W, six_node = FALSE) {
  # mc is a vector of motif counts
  M <- W
  M[M > 0] <- 1

  ###### compute all_mot_perm
  all_mot_perm <- gen_all_mot_perm()
  m <- nrow(W)
  n <- ncol(W)

  if (six_node) {
    mc <- rep(0, 44)
  }
  if (!six_node) {
    mc <- rep(0, 17)
  }


  for (i in 1:m) {
    # i is now nrow of the minor
    rp <- gtools::combinations(m, i, 1:m, repeats = FALSE)
    for (j in 1:n) {
      # j is now ncol of the minor
      if (i + j > (5 + six_node)) {
        next
      }
      cp <- gtools::combinations(n, j, 1:n, repeats = FALSE)
      nrp <- nrow(rp)
      ncp <- nrow(cp)
      for (k in 1:nrp) {
        for (l in 1:ncp) {
          minor <- W[rp[k,], cp[l,], drop = FALSE]
          bin_minor <- minor
          bin_minor[bin_minor > 0] <- 1

          mt <- check_motif(bin_minor, all_mot_perm)
          if (mt == -1) {
            next
          }
          mc[mt] <- mc[mt] + 1
        }
      }
    }
  }

  mc
}

gen_all_mot_perm <- function(all_mot = -1) {
  # all_mot should be a list of all motifs that we want to permute
  # the function generates all permutations of the motifs
  # returns a nested list
  # all_mot_perm[[i]] is a list of all permutations of the element i in all_mot
  # if called without any argument, will load in data

  if (all_mot == -1) {
    # load('all_motifs.RData')
    all_mot <- motifs
  }

  # library(gtools)
  all_mot_perm <- list()
  v <- list()
  v[[1]] <- NA
  v[[2]] <- 1
  v[[3]] <- 2:3
  v[[4]] <- 4:7
  v[[5]] <- 8:17
  v[[6]] <- 18:44

  # loop through possible motif sizes
  for (n in 2:6) {
    perm <- list()

    # loop through poss number of vertices in first class
    for (i in 1:(n-1)) {
      # now have i vertices in first class, n-i in second
      pr <- gtools::permutations(i, i, 1:i, repeats = FALSE)
      pc <- gtools::permutations(n-i, n-i, 1:(n-i), repeats = FALSE)
      perm[[i]] <- list(pr, pc)
    }

    # now permute all the motifs
    for (j in v[[n]]) {
      mot <- all_mot[[j]]
      k <- nrow(mot)
      pr <- perm[[k]][[1]]
      pc <- perm[[k]][[2]]
      all_mot_perm[[j]] <- bipgraph_isom_class(mot, pr, pc)
    }
  }
  all_mot_perm
}

bipgraph_isom_class <- function(M, poss_perm_rows, poss_perm_cols ) {
  # M is adjacency matrix of a bipartite graph (i.e. not nec. square matrix)
  # now we want to find all bipartite graphs obtained by just permututing rows and columns
  # poss_perm_rows is a matrix with possible permutations of the row indices
  # poss_perm_cols ... same for column indices
  # we return a list with all isomorphic bipartite graphs

  # poss_perm_rows <- permutations(N, N, 1:N, repeats = FALSE)
  # poss_perm_cols <- permutations(M, M, 1:M, repeats = FALSE)

  # these are given as arguments to this function so we don't have to compute them again, every time the function is called

  I <- nrow(poss_perm_rows)
  J <- nrow(poss_perm_cols)

  res <- list()

  for (i in 1:I) {
    M_pr <- M[poss_perm_rows[i,], , drop = FALSE]
    for (j in 1:J) {
      M_pc <- M_pr[,poss_perm_cols[j,], drop = FALSE]
      len <- length(res)
      res[[len + 1]] <- M_pc
    }
  }

  res <- unique(res)
  return(res)

}

check_motif <- function(minor, all_mot_perm) {
  # all_mot_perm[[i]] is a list of all permutations of motif i
  # minor is a subgraph
  # function returns the motif type of the given minor

  # delete row and column names
  dimnames(minor) <- NULL

  # don't have to try all motifs, only the ones of same dimension
  poss_mot <- motif_dimensions(nrow(minor), ncol(minor))

  for (i in poss_mot) {
    iso_class <- all_mot_perm[[i]]
    for (j in 1:length(iso_class)) {
      candidate <- iso_class[[j]]
      dimnames(candidate) <- NULL
      if (identical(minor, candidate)) {
        return(i)
      }
    }
  }
  return(-1)
}

motif_dimensions <- function(m, n) {
  # This function takes as arguments the dimensions of a subgraph
  # It returns a vector with all possible motif types of the same dimension
  if (m == 1 & n == 1) {
    return(1)
  }
  if (m == 1 & n == 2) {
    return(2)
  }
  if (m == 2 & n == 1) {
    return(3)
  }
  if (m == 3 & n == 1) {
    return(4)
  }
  if (m == 2 & n == 2) {
    return(5:6)
  }
  if (m == 1 & n == 3) {
    return(7)
  }
  if (m == 4 & n == 1) {
    return(8)
  }
  if (m == 3 & n == 2) {
    return(9:12)
  }
  if (m == 2 & n == 3) {
    return(13:16)
  }
  if (m == 1 & n == 4) {
    return(17)
  }
  if (m == 5 & n == 1) {
    return(18)
  }
  if (m == 4 & n == 2) {
    return(19:24)
  }
  if (m == 3 & n == 3) {
    return(25:37)
  }
  if (m == 2 & n == 4) {
    return(38:43)
  }
  if (m == 1 & n == 5) {
    return(44)
  }

}
