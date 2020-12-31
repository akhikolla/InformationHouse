Lmoments <- function(data, rmax = 4, na.rm = FALSE, returnobject = FALSE, 
                     trim = c(0, 0)) {
  
  if (!identical(trim, c(0, 0)) & !identical(trim, c(1, 1))) {
    stop("The current version of Lmoments supports only ordinary L-moments (trim=c(0,0)) and T1L-moments (trim=c(1,1))")
  }
  if (identical(trim, c(1, 1)) & rmax > 4) {
    warning("The current version of t1lmoments uses rmax=4.")
    rmax <- 4
  }
  data <- as.matrix(data)
  rmax <- min(rmax, nrow(data))
  p <- dim(data)[2]
  if (!na.rm) {
    if (identical(trim, c(0, 0))) 
      L <- Lmoments_calc(data, rmax)
    if (identical(trim, c(1, 1))) 
      L <- t1lmoments(data, rmax)
  }
  else {
    L <- array(0.0, c(p, rmax))
    for (i in 1:p) {
      xi <- data[, i]
      xi <- xi[!is.na(xi)]
      if (identical(trim, c(0, 0))) 
        L[i, ] <- Lmoments_calc(xi, rmax)
      if (identical(trim, c(1, 1))) 
        L[i, ] <- t1lmoments(xi, rmax)
    }
  }
  colnames(L) <- paste("L", 1:rmax, sep = "")
  if (p > 1) 
    rownames(L) <- names(data)
  if (returnobject) {
    Lobject <- list()
    Lobject$lambdas <- L
    Lobject$trim <- trim
    if (identical(trim, c(0, 0))) 
      Lobject$source <- "Lmoments"
    if (identical(trim, c(1, 1))) 
      Lobject$source <- "t1lmoments"
    if (rmax > 2) {
      ratios <- L[, 3:rmax]/(cbind(L[, 2]) %*% rep(1, (rmax - 2)))
      colnames(ratios) <- paste("tau", 3:rmax, sep = "")
      if (p == 1) 
        ratios <- cbind(rbind(L[, 1:2]), cbind(ratios))
      if (p > 1) 
        ratios <- cbind(L[, 1:2], ratios)
      if (p > 1) 
        rownames(ratios) <- names(data)
      Lobject$ratios <- ratios
    }
    return(Lobject)
  }
  else {
    return(L)
  }
}

Lmomcov <- function (data, rmax = 4, na.rm = FALSE) {
  if (rmax > 5) warning(sprintf("The numerical accuracy of the results is likely less than %s with rmax > 5.",
                                formatC(sqrt(.Machine$double.eps), format = "g", digits = 3)))
  if (rmax <= 1) return(NA)
  data <- as.matrix(data)
  rmax <- min(rmax, nrow(data))
  if (!na.rm) {
    res <- Lmomcov_calc(data, rmax)
    if (length(res) == 1) {
      rownames(res[[1]]) <- paste("L", 1:rmax, sep = "")
      colnames(res[[1]]) <- paste("L", 1:rmax, sep = "")
      return(res[[1]])
    }
    return(lapply(res, function(covmatrix) {
      rownames(covmatrix) <- paste("L", 1:rmax, sep = "")
      colnames(covmatrix) <- paste("L", 1:rmax, sep = "")
      covmatrix
    }))
  }
  else {
    p <- dim(data)[2]
    if (p == 1) {
      covmatrix <- Lmomcov_calc(data[!is.na(data)], rmax)[[1]]
      rownames(covmatrix) <- paste("L", 1:rmax, sep = "")
      colnames(covmatrix) <- paste("L", 1:rmax, sep = "")
      return(covmatrix)
    }
    else {
      covmatrixlist <- list()
      for (i in 1:p) {
        xi <- data[, i]
        xi <- xi[!is.na(xi)]
        covmatrixlist[[i]] <- Lmomcov_calc(xi, rmax)
        rownames(covmatrix[[i]]) <- paste("L", 1:rmax, sep = "")
        colnames(covmatrix[[i]]) <- paste("L", 1:rmax, sep = "")
      }
      names(covmatrixlist) <- names(data)
      return(covmatrixlist)
    }
  }
}

Lcoefs<-function(data, rmax=4, na.rm=FALSE, trim=c(0,0))
{
  Lobject<-Lmoments(data,rmax=rmax,na.rm=na.rm,returnobject=TRUE,trim=trim);
  return(Lobject$ratios)
}

t1lmoments <- function(data, rmax = 4) {
  rmax_out <- min(rmax, 4)
  if (rmax != 4) {
    if (rmax > 4) 
      warning("The current version of t1lmoments uses rmax = 4.")
    rmax <- 4
  }
  data <- as.matrix(data)
  if (nrow(data) < rmax) stop("Insufficient sample size to estimate trimmed L-moments.")
  t1lmoments_calc(data, rmax)[, seq_len(rmax_out)]
}
