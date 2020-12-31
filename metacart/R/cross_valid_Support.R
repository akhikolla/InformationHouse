#' A function to parition newdata based on the fitted model
#' @param x a fitted metaCART model
#' @param newdata new data for prediction
#' @keywords internal
#' @importFrom stats delete.response model.frame
prednode_cpp <- function(x, newdata) {
  tt <- terms(x$data)
  new.ms <- model.frame(delete.response(tt), newdata) # the moderators in the new data
  old.ms <- model.frame(delete.response(tt), x$data) # the moderators in the original data set
  tree <- x[["tree"]]
  nameM <- 1:ncol(new.ms)
  names(nameM) <- colnames(new.ms)
  inxM <- nameM[as.character(tree$mod)]
  boolName <- sapply(old.ms, is.numeric)  #tell if a moderator is numeric
  names(boolName) <- colnames(old.ms)
  boolNumeric <- boolName[as.character(tree$mod)]
  .partition(tree, new.ms, boolNumeric, inxM, x$cpt, old.ms)
  # Note: return NA values if the new data contain categories that does not exist in the original data
}

#' A function to compute cross-validation errors
#' @param func the function used to fit the model
#' @param mf the data set with formula specified
#' @param maxL the maximum number of splits
#' @param n.fold the number of folds
#' @param minbucket the minimum number of the studies in a terminal node
#' @param minsplit the minimal number of studies in a parent node to be split
#' @param cp the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @param lookahead an argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @keywords internal
#' @importFrom stats terms model.response
Xvalid_all <- function(func, mf, maxL, n.fold, minbucket, minsplit, cp, lookahead){
  N <- nrow(mf)
  if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
  if (maxL > N) {
    warning("The maximum number of split is too large")
    maxL <- N-1
  }
  
  pred <- matrix(NA, nrow = N, ncol = maxL+1)
  inx <- sample(1:N)
  inx.xvalid <- as.numeric(cut(1:N, n.fold))
  for (i in 1:n.fold){
    inx.test <- inx[inx.xvalid == i]
    test <- mf[inx.test, ]  # test set
    train <- mf[-inx.test, ]  # training set
    fit.train <- do.call(func, list(mf = train, maxL, minbucket, minsplit, cp, lookahead))
    yi.train <- model.response(fit.train$data)
    vi.train <- c(t(fit.train$data["(vi)"]))
    tau2 <- fit.train$tree$tau2
    nsplt <- nrow(fit.train$tree)
    # compute the effect size for all internal and terminal nodes
    train.y <- .ComputeY(fit.train$node.split, yi.train, vi.train, tau2) 
    # predict the node membership for the test set
    test.nodes <- prednode_cpp(fit.train, test)
    # obtain the predicted effect size
    test.y <- .PredY(train.y, test.nodes)
    # missing values are replaced by the overall weighted mean
    if (any(is.na(test.y))) {
      inx.NA <- which(is.na(test.y), arr.ind=TRUE)
      test.y <- .ReplaceNA(inx.NA, test.y, yi.train, vi.train, tau2)
    }
    pred[inx.test, 1:nsplt] <- test.y
  }
  y <- model.response(mf)
  if (!is.null(dim(pred))) { 
    x.error <- apply(pred,2, function(x) sum((y-x)^2)/sum((y-mean(y))^2))
    sdx.error <- apply(pred, 2, function(x)  sqrt(sum(((y-x)^2-mean((y-x)^2))^2))/sum((y-mean(y))^2))
  } else {  # root node with no splits
    x.error <- sum((y-pred)^2)/sum((y-mean(y))^2)
    sdx.error <- sqrt(sum(((y-pred)^2-mean((y-pred)^2))^2))/sum((y-mean(y))^2)
  }
  cbind(x.error, sdx.error)
}
