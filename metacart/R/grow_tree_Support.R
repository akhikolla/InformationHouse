#' A function to list all possible split points for the first split
#' 
#' @param xk moderator
#' @param y effect size
#' @return childNodes: child nodes membership
#' @return rank: the rank used to order categorical moderator 
#' @return is.num: if the moderator is numeric
#' @keywords internal
make_first_split <- function(xk, y) {
  if (is.numeric(xk)) {
    split.points <- sort(unique(xk))
    res <- list(childNodes = sapply(split.points[-1], function(x) ifelse(xk < x, 2, 3)),
                rank1 = split.points, is.num = TRUE)
  } else {
    split.points <- sort(rank(tapply(y, xk, mean)))  # rank the categories by group mean
    res <- list(childNodes = sapply(split.points[-1], function(x) ifelse(split.points[as.character(xk)] < x, 2, 3)),
                rank1 = split.points, is.num = FALSE)
  }
}

#'  A function to find the best triplets of parent, moderator, and split point.
#'  
#' @param xk moderator vector
#' @param nodeMbrship node membership vector
#' @param y effect size
#' @param vi sampling variance
#' @param minbucket the minimal number of studies in child nodes
#' @param minsplit the minimal number of studies in parent node
#' @return pleaf: the parent node
#' @return cstar: the split point
#' @return rank: the rank used to order categorical moderator 
#' @keywords internal
find_triplet <- function(xk, nodeMbrship, y, vi, minbucket, minsplit) {
  leaves <- as.numeric(names(table(nodeMbrship)[table(nodeMbrship) >= minsplit]))
  if (is.numeric(xk)) {
    xk.rank <- NULL
    tempQ <-  sapply(leaves, function(x)
      re.cutoff_cpp(y, vi, xk[nodeMbrship == x], nodeMbrship == x, nodeMbrship, minbucket)
    )
  } else {
    xk.rank <- lapply(leaves,function(x) rank(tapply(y[nodeMbrship == x], xk[nodeMbrship == x], mean)))
    names(xk.rank) <- leaves
    tempQ <- sapply(leaves, function(x) 
      re.cutoff_cpp(y, vi, xk.rank[[as.character(x)]][as.character(xk[nodeMbrship == (x)])], 
                    nodeMbrship == x, nodeMbrship, minbucket) )
  }
  if(class(tempQ) == "list") {  # At least one split is not eligible
    if (all(sapply(tempQ, is.null))) {  # no eligible splits
      list(pleaf = NA, cstar = c(NA, -Inf, NA), rank = NULL)
    } else {  # if one split is not eligible (NULL), make its Q-between = -Inf
      pleaf <- which.max(sapply(tempQ, function(x) if (is.null(x)) -Inf else x[2]))
      list(pleaf = leaves[pleaf], cstar = tempQ[[pleaf]], rank = xk.rank[[pleaf]])
    }
    
  } else {  # all splits are eligible 
    pleaf <- which.max(tempQ[2,])
    list(pleaf = leaves[pleaf], 
         cstar = tempQ[ ,pleaf],
         rank = xk.rank[[pleaf]])
  }
}



#' A function to find the optimal combination of first two splits, 
#' and the corresponding Q-between given the first split
#' @param xk moderator vector
#' @param first.splits possible first splits
#' @param y effect size
#' @param vi sampling variance
#' @param minbucket the minimal number of studies in child nodes
#' @param minsplit the minimal number of studies in parent node
#' @return a list including all possible combinations of the triplet
#' @keywords internal
find_second_split <- function(xk, first.splits, y, vi, minbucket, minsplit) {
  # Given a moderator and all possible first splits
  # Return the optimal combination of first two splits, and the corresponding Q-between
  # split1: the first element is the index of the first moderator
  #         the second element is the split point of the first moderator
  # split2: the first element is the selected parent leaf for the 2nd split
  #         the second element contains the split point of the 2nd split, Q, the tau-square
  Q.temp <- -Inf
  Csplit.temp <- list(Q = Q.temp)
  for (i in 1:length(first.splits)) {
    temp.first.split <- first.splits[[i]]$childNodes
    if (length(temp.first.split) == 0) next
    res <- lapply(1:ncol(temp.first.split), function(x) find_triplet(xk, temp.first.split[ ,x], y, vi, minbucket, minsplit))
    tempQs <- sapply(1:ncol(temp.first.split), function(x) {if (is.null(res[[x]])) {-Inf} else {res[[x]]$cstar[2]}})
    inx.max <- which.max(tempQs)
    Q.max <- tempQs[inx.max]
    Csplit1.max <- c(i, (first.splits[[i]]$rank1[inx.max] + first.splits[[i]]$rank1[inx.max+1])/2, inx.max)
    Csplit2.max <- res[[inx.max]]
    if (Q.max > Q.temp) {
      Q.temp <- Q.max
      Csplit.temp <- list(Q = Q.temp,
                          split1 = Csplit1.max, 
                          split2 = Csplit2.max)
    }
  }
  Csplit.temp
}

#' A function to compute RE Q-between
#' @param y effect size
#' @param vi sampling variance
#' @param xk moderator 
#' @return Q-between and tau2
#' @keywords internal
compute_rebetQ <- function(y, vi, xk){
  wts = 1/vi
  wy = wts*y
  wy2 = wts * y^2
  Q <- tapply(wy2, xk, sum) - tapply(wy, xk, function(x) (sum(x))^2)/tapply(wts, xk, sum)
  df <- tapply(wy, xk, length)-1
  C <- tapply(wts, xk, sum) - tapply(wts, xk, function(x) sum(x^2))/ tapply(wts, xk, sum)
  tau2 <- (sum(Q) - sum(df))/sum(C)
  tau2 <- max(0, tau2)
  wstar = 1/(vi+tau2)
  wystar = wstar*y
  wy2star = wstar*y^2
  Qstar <- tapply(wy2star, xk, sum) - tapply(wystar, xk, function(x) (sum(x))^2)/tapply(wstar, xk, sum)
  Qstar.total <- sum(wy2star) - (sum(wystar))^2/sum(wstar)
  Qbet <- Qstar.total - sum(Qstar)
  if (is.na(Qbet)) {
    Qbet <- Inf
  }
  return(c(Qbet, tau2))
  
}
