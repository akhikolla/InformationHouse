mean_weight <- function(W, mc = NA, six_node = FALSE) {
  # W is the weighted adjacency matrix
  # mc a vector of motif counts
  # returns a vector
  # position i is the mean weight for motif i

  # assume W is the weighted adj matrix
  # M is the binary adj
  M <- W
  M[M > 0] <- 1

  if (any(is.na(mc))) {
    # user did not give any input for this, so we compute it again
    mc <- mcount (M, six_node = six_node, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  }
  # make sure mc has the right length
  if (!six_node) {
    mc <- mc[1:17]
  }

  # give M default rownames, otherwise the reordering later does not work
  rownames(M) <- paste0("r", 1:nrow(M))
  colnames(M) <- paste0("c", 1:ncol(M))

  rownames(W) <- paste0("r", 1:nrow(M))
  colnames(W) <- paste0("c", 1:ncol(M))

  lr <- link_positions(M, normalisation = "none", weights = FALSE, six_node = six_node)
  colnames(lr) <- 1:ncol(lr)

  # now create a list of edges with weights
  # have to transpose so the list is ordered by rows, not by columns
  mw <- reshape2::melt(t(W))
  # note: now Var2 gives me the row name
  # Var1 gives me the column name
  mw[, 'link'] <- paste(mw$Var2, mw$Var1, sep = " -- ")
  # mw[, 'link'] <- paste("(", mw$Var2, ", ", mw$Var1, ")", sep = "")
  mw2 <- mw [,c('link', 'value')]

  # now clean the links that are not present in the network
  # wl (weighted links) is now the list we want
  wl <- subset(mw2, mw2$value != 0)

  # potentially reorder
  if (!identical(rownames(lr),wl$link)) {
    mtch <- match(rownames(lr), wl$link)
    wl <- wl[mtch, ]
  }

  # want: matrix which sums over the motifs, i.e.
  # each row: one link
  # each column: one motif
  # entry gives the count for that link in that motif (in any position)
  nl <- nrow(lr) # number of links

  if (six_node) {
    nc <- 44
  } else {
    nc <- 17
  }

  lr_sum <- matrix(rep(NA, nc*nl), nrow = nl)

  for (i in 1:nc) {
    sel <- lr[, motif_info(i, node = FALSE), drop = FALSE]
    lr_sum[,i] <- rowSums(sel)
  }

  rownames(lr_sum) <- rownames(lr)
  # multiply each frequency with the weight of the link
  lr_w <- apply(lr_sum, 2, function(x) {x * wl$value})
  if(inherits(lr_w,"numeric") == TRUE) {
    lr_w <- t(as.matrix(lr_w))
  }

  # now we are interested in the column sums
  # sum over column i gives the total weight of all submotifs of motif-type i

  tot_w <- colSums(lr_w)
  mean_w <- tot_w / mc
  # If the motif does not occur, we want the mean weight to be NA
  mean_w <- replace(mean_w, which(is.nan(mean_w)), NA)

  # now divide by the number of links in each network
  num_link <- c(1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 5, 6, 4, 4, 5, 6, 4)
  if (six_node) {
    num_link <- c(num_link, 5, 5, 5, 6, 6, 7, 8, 5, 5, 5, 5, 6, 6, 7, 6, 7, 6, 7, 8, 9, 5, 5, 6, 6, 7, 8, 5)
  }
  return(mean_w / num_link)
}
