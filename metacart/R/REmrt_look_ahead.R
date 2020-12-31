#' A function to fit the tree with look-ahead option
#'
#' @param mf the data.frame to grow the tree
#' @param maxL the maximum number of splits
#' @param minbucket the minimum number of the studies in a terminal node
#' @param minsplit the minimal number of studies in a parent node to be split
#' @param cp the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @param lookahead an argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @return a list including a tree, the split points, the data, and the nodes after each split
#' @keywords internal
#' @importFrom stats terms model.response
REmrt_GS_ <- function(mf, maxL, minbucket, minsplit, cp, lookahead){
  #===================  Error message  ======================#
  if (minbucket >= minsplit) {
    stop("minbucket should be smaller than minsplit")
  }
  if (!is.logical(lookahead)) {
    stop("lookahead should be TRUE or FALSE")
  }
  if (lookahead & maxL < 2) {
    stop("the maximum number of splits should be larger than two when using lookahead")
  }
  
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(rep(1, nrow(mf)))
  res.Qb = 0
  res.tau2 = (sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)   - length(y)+1)/(sum(1/vi)-sum(1/vi^2)/sum(1/vi)) #VERIFIED
  res.split = NA
  res.mod = NA
  res.pleaf = NA
  delta.Q <- Inf
  if (lookahead) {
    #-------           Make the first two splits          -------#
    first.splits <- lapply(mods, function(x) make_first_split(x, y))
    second.splits <- lapply(1:nmod, function(x) find_second_split(mods[ ,x], first.splits, y, vi, minbucket, minsplit))
    tempQs <- sapply(1:length(second.splits), function(x) second.splits[[x]]$Q)
    x2.inx <- which.max(tempQs)  # choose the combination with largest Q
    x1.inx <- second.splits[[x2.inx]]$split1[1]
    res.mod <- c(NA, mods.names[c(x1.inx, x2.inx)])
    res.pleaf <- c(NA, 1, second.splits[[x2.inx]]$split2$pleaf)
    if (first.splits[[x1.inx]]$is.num) {  # Numeric variable
      cstar1 <- second.splits[[x2.inx]]$split1[2]
      msplit1 <- paste(mods.names[x1.inx], "<", signif(cstar1,2), collapse = " ")
    } else {
      cstar1 <- names(first.splits[[x1.inx]]$rank[first.splits[[x1.inx]]$rank < second.splits[[x2.inx]]$split1[2]])
      msplit1 <- paste(mods.names[x1.inx], "=", paste(cstar1, collapse = "/"), collapse = " ")
    }
    node1 <- first.splits[[x1.inx]]$childNodes[ ,second.splits[[x2.inx]]$split1[3]]
    names(node1) <- NULL
    pleaf.inx <- node1 == second.splits[[x2.inx]]$split2$pleaf
    node2 <- node1
    if (is.null(second.splits[[x2.inx]]$split2$rank)) {
      cstar2 <- second.splits[[x2.inx]]$split2$cstar[1]
      node2[pleaf.inx] <- ifelse(mods[pleaf.inx, x2.inx] < cstar2, 4, 5)
      msplit2 <- paste(mods.names[x2.inx], "<", signif(cstar2,2) , collapse = " ")
    } else {
      cstar2 <- names(second.splits[[x2.inx]]$split2$rank[second.splits[[x2.inx]]$split2$rank < second.splits[[x2.inx]]$split2$cstar[1]])
      node2[pleaf.inx] <- ifelse(mods[pleaf.inx, x2.inx] %in% cstar2, 4, 5)
      msplit2 <- paste(mods.names[x2.inx], "=", paste(cstar2, collapse = "/"), collapse = " ")
    }
    nodemark <- cbind(nodemark, node1, node2)
    cpt[[1]] <- cstar1
    cpt[[2]] <- cstar2
    Q.split1 <- compute_rebetQ(y, vi, node1)
    res.Qb <- c(0, Q.split1[1], second.splits[[x2.inx]]$Q)
    res.tau2 <- c(res.tau2, Q.split1[2], second.splits[[x2.inx]]$split2$cstar[3])
    res.split <- c(NA, msplit1, msplit2)
    Dev <- res.Qb[3]
    delta.Q <- res.Qb[3] - res.Qb[2]
    i = 2
  } else{
    for (i in 1) {
      Dev<- -Inf
      TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
      cnode <- nodemark[ ,i]
      len.node <- tapply(vi, cnode, length)
      nodes <- names(len.node) [len.node >= minsplit]
      for (pl in nodes) {
        pleaf.inx <- cnode == pl
        for (k in 1:nmod) {
          xk <- mods[pleaf.inx, k]
          c.splits <- unique(xk)
          if (length(c.splits) < 2) next
          if (is.numeric(xk)) {
            # NUMERIC VARIABLE
            temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- temp[1]
              msplit <- paste(mods.names[k], "<", signif(c.star,2), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk < c.star, 2*i, 2*i+1)
            }
          } else {
            xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
            xk.ordinal <- xk.rank[as.character(xk)]
            temp <- re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- names(xk.rank[xk.rank < temp[1]])
              msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
            }
          }
        }
      }
      if (is.null(TQb)) { 
        delta.Q <- -Inf
      } else {
        delta.Q <- abs(TQb - res.Qb[i])
      }
    }
    nodemark <- cbind(nodemark, new.node)
    res.Qb <- c(res.Qb, TQb)
    res.tau2 <- c(res.tau2, Ttau2)
    res.split <- c(res.split, Tsplit)
    res.mod <- c(res.mod,Tmod)
    res.pleaf <- c(res.pleaf, Tpleaf)
    cpt[[i]] <- c.star
  }
  #-------          Continue with greedy search          -------#
  while(delta.Q >= cp & i < maxL) {
      i <- i+1
      TQb <- Ttau2 <- Tsplit <- Tmod <- Tpleaf <- NULL
      cnode <- nodemark[ ,i]
      len.node <- tapply(vi, cnode, length)
      nodes <- names(len.node) [len.node >= minsplit]
      for (pl in nodes) {
        pleaf.inx <- cnode == pl
        for (k in 1:nmod) {
          xk <- mods[pleaf.inx, k]
          c.splits <- unique(xk)
          if (length(c.splits) < 2) next
          if (is.numeric(xk)) {
            # NUMERIC VARIABLE
            temp <- re.cutoff_cpp(y, vi, xk, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- temp[1]
              msplit <- paste(mods.names[k], "<", signif(c.star,2), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk < c.star, 2*i, 2*i+1)
            }
          } else {
            xk.rank <- rank(tapply(y[pleaf.inx], xk, mean))
            xk.ordinal <- xk.rank[as.character(xk)]
            temp <-  re.cutoff_cpp(y, vi, xk.ordinal, pleaf.inx, cnode, minbucket)
            if (is.null(temp)) {
              Dev.new <- -Inf
            } else {
              Dev.new <- temp[2]
            }
            if (Dev.new > Dev) {
              Dev <- temp[2]
              c.star <- names(xk.rank[xk.rank < temp[1]])
              msplit <- paste(mods.names[k], "=", paste(c.star, collapse = "/"), collapse = " ")
              TQb = temp[2]
              Ttau2 = temp[3]
              Tsplit = msplit
              Tmod = mods.names[k]
              Tpleaf = as.numeric(pl)
              new.node <- cnode
              new.node[pleaf.inx] <- ifelse( xk.ordinal < temp[1], 2*i, 2*i+1)
            }
          }
        }
      }
      if (is.null(TQb)) { 
        delta.Q <- -Inf
      } else {
        nodemark <- cbind(nodemark, new.node)
        res.Qb <- c(res.Qb, TQb)
        res.tau2 <- c(res.tau2, Ttau2)
        res.split <- c(res.split, Tsplit)
        res.mod <- c(res.mod,Tmod)
        res.pleaf <- c(res.pleaf, Tpleaf)
        cpt[[i]] <- c.star
        delta.Q <- abs(TQb - res.Qb[i])
      }
      
    }
  list(tree = data.frame(Qb = res.Qb, tau2 = res.tau2, split = res.split,
                         mod = res.mod, pleaf = res.pleaf, stringsAsFactors = FALSE),
       node.split = nodemark, cpt = cpt, data = mf)
  
}

