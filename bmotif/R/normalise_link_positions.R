normalise_link_positions <- function(lp, type, six_node) {

  n <- 29 + six_node *(106 - 29)
  if(inherits(lp,"matrix") != TRUE) {
    lp <- as.matrix(lp, ncol = n)
  }

  if (type == "sum") {
    # divide by total number of times that each edge occurs in any position
    # the required number is exactly the rowsum

    # if lp only has one row, R does weird things later, so we save the row names
    if (nrow (lp) == 1) {
      rn <- rownames(lp)
    }

    s <- rowSums(lp)
    lp <- as.matrix(apply(lp, 2, function(x) {x / s}))

    # if lp had one row previously, R will now have transposed the result and killed the row names
    if (ncol (lp) == 1) {
      lp <- t(lp)
      rownames(lp) <- rn
    }

    lp <- replace(lp, which(is.nan(lp)), NA)
    return(lp)
  } else if (type == "sizeclass") {

    # note: matrix / vector divides each column of the matrix by the vector elementwise
    lp[,2:3] <- lp[,2:3, drop = FALSE] / rowSums(lp[,2:3, drop = FALSE])
    lp[, 4:9] <- lp[,4:9, drop = FALSE] / rowSums(lp[,4:9, drop = FALSE])
    lp[, 10:29] <- lp[,10:29, drop = FALSE] / rowSums(lp[,10:29, drop = FALSE])

    if (six_node) {
      lp[, 30:106] <- lp[,30:106, drop = FALSE] / rowSums(lp[,30:106, drop = FALSE])
    }

    lp <- replace(lp, which(is.nan(lp)), NA)
    return(lp)

  } else if (type == "sizeclass_plus1") {

    lp <- lp + 1
    # note: matrix / vector divides each column of the matrix by the vector elementwise
    lp[,2:3] <- lp[,2:3, drop = FALSE] / rowSums(lp[,2:3, drop = FALSE])
    lp[, 4:9] <- lp[,4:9, drop = FALSE] / rowSums(lp[,4:9, drop = FALSE])
    lp[, 10:29] <- lp[,10:29, drop = FALSE] / rowSums(lp[,10:29, drop = FALSE])

    if (six_node) {
      lp[, 30:106] <- lp[,30:106, drop = FALSE] / rowSums(lp[,30:106, drop = FALSE])
    }

    return(lp)

  } else if (type == "sizeclass_NAzero") {

    # note: matrix / vector divides each column of the matrix by the vector elementwise
    lp[,2:3] <- lp[,2:3, drop = FALSE] / rowSums(lp[,2:3, drop = FALSE])
    lp[, 4:9] <- lp[,4:9, drop = FALSE] / rowSums(lp[,4:9, drop = FALSE])
    lp[, 10:29] <- lp[,10:29, drop = FALSE] / rowSums(lp[,10:29, drop = FALSE])

    if (six_node) {
      lp[, 30:106] <- lp[,30:106, drop = FALSE] / rowSums(lp[,30:106, drop = FALSE])
    }

    lp <- replace(lp, which(is.na(lp)), 0)
    return(lp)

  } else if (type == "levelsize") {

    lp[,2:4] <- ifelse(lp[,2:4] == 0, yes = NA, no = 1)
    lp[,5:8] <- lp[,5:8, drop = FALSE] / rowSums(lp[,5:8, drop = FALSE])
    lp[,9:10] <- ifelse(lp[,9:10] == 0, yes = NA, no = 1)
    lp[,11:19] <- lp[,11:19, drop = FALSE] / rowSums(lp[,11:19, drop = FALSE])
    lp[,20:28] <- lp[,20:28, drop = FALSE] / rowSums(lp[,20:28, drop = FALSE])
    lp[,29] <- ifelse(lp[,29] == 0, yes = NA, no = 1)

    if (six_node) {
      lp[,30] <- ifelse(lp[,30] == 0, yes = NA, no = 1)
      lp[,31:46] <- lp[,31:46, drop = FALSE] / rowSums(lp[,31:46, drop = FALSE])
      lp[,47:89] <- lp[,47:89, drop = FALSE] / rowSums(lp[,47:89, drop = FALSE])
      lp[,90:105] <- lp[,90:105, drop = FALSE] / rowSums(lp[,90:105, drop = FALSE])
      lp[,106] <- ifelse(lp[,106] == 0, yes = NA, no = 1)
    }

    lp <- replace(lp, which(is.nan(lp)), NA)
    return(lp)

  } else if (type == "levelsize_plus1") {

    lp <- lp + 1
    lp[,2:4] <- ifelse(lp[,2:4] == 0, yes = NA, no = 1)
    lp[,5:8] <- lp[,5:8, drop = FALSE] / rowSums(lp[,5:8, drop = FALSE])
    lp[,9:10] <- ifelse(lp[,9:10] == 0, yes = NA, no = 1)
    lp[,11:19] <- lp[,11:19, drop = FALSE] / rowSums(lp[,11:19, drop = FALSE])
    lp[,20:28] <- lp[,20:28, drop = FALSE] / rowSums(lp[,20:28, drop = FALSE])
    lp[,29] <- ifelse(lp[,29] == 0, yes = NA, no = 1)

    if (six_node) {
      lp[,30] <- ifelse(lp[,30] == 0, yes = NA, no = 1)
      lp[,31:46] <- lp[,31:46, drop = FALSE] / rowSums(lp[,31:46, drop = FALSE])
      lp[,47:89] <- lp[,47:89, drop = FALSE] / rowSums(lp[,47:89, drop = FALSE])
      lp[,90:105] <- lp[,90:105, drop = FALSE] / rowSums(lp[,90:105, drop = FALSE])
      lp[,106] <- ifelse(lp[,106] == 0, yes = NA, no = 1)
    }

    return(lp)

  } else if (type == "levelsize_NAzero") {

    lp[,2:4] <- ifelse(lp[,2:4] == 0, yes = NA, no = 1)
    lp[,5:8] <- lp[,5:8, drop = FALSE] / rowSums(lp[,5:8, drop = FALSE])
    lp[,9:10] <- ifelse(lp[,9:10] == 0, yes = NA, no = 1)
    lp[,11:19] <- lp[,11:19, drop = FALSE] / rowSums(lp[,11:19, drop = FALSE])
    lp[,20:28] <- lp[,20:28, drop = FALSE] / rowSums(lp[,20:28, drop = FALSE])
    lp[,29] <- ifelse(lp[,29] == 0, yes = NA, no = 1)

    if (six_node) {
      lp[,30] <- ifelse(lp[,30] == 0, yes = NA, no = 1)
      lp[,31:46] <- lp[,31:46, drop = FALSE] / rowSums(lp[,31:46, drop = FALSE])
      lp[,47:89] <- lp[,47:89, drop = FALSE] / rowSums(lp[,47:89, drop = FALSE])
      lp[,90:105] <- lp[,90:105, drop = FALSE] / rowSums(lp[,90:105, drop = FALSE])
      lp[,106] <- ifelse(lp[,106] == 0, yes = NA, no = 1)
    }

    lp <- replace(lp, which(is.na(lp)), 0)
    return(lp)

  } else if (type == "motif") {
    nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
    mp <- sapply(1:nm, function(x) motif_info(x, node = FALSE)) # motif positions
    for (i in mp) {
      lp[,i] <- lp[,i, drop = FALSE] / rowSums(lp[,i, drop = FALSE])
    }
    lp <- replace(lp, which(is.nan(lp)), NA)
    return(lp)

  } else if (type == "motif_plus1") {

    lp <- lp + 1
    nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
    mp <- sapply(1:nm, function(x) motif_info(x, node = FALSE)) # motif positions
    for (i in mp) {
      lp[,i] <- lp[,i, drop = FALSE] / rowSums(lp[,i, drop = FALSE])
    }

    return(lp)

  } else if (type == "motif_NAzero") {

    nm <- ifelse(six_node, yes = 44, no = 17) # number of motifs
    mp <- sapply(1:nm, function(x) motif_info(x, node = FALSE)) # motif positions
    for (i in mp) {
      lp[,i] <- lp[,i, drop = FALSE] / rowSums(lp[,i, drop = FALSE])
    }
    lp <- replace(lp, which(is.nan(lp)), 0)
    return(lp)
  } else if(type == "position"){
    # divide by total number of times this edge position occurs
    # i.e. columnSums
    s <- colSums(lp)
    lp <- t(apply(lp, 1, function(x) {x / s}))
    lp <- replace(lp, which(is.nan(lp)), NA)
    return(lp)
  } else {
    stop("No valid normalisation method specified")
  }

}
