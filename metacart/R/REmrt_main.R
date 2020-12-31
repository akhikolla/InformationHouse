#' Random effects meta-tree
#'
#'A function to fit a random effects meta-tree
#' @name REmrt
#' @aliases REmrt
#' @param formula A formula, with a response variable (usually the effect size) and the potential moderator variables but no interaction terms.
#' @param vi sampling variance of the effect size.
#' @param data A data frame of a meta-analytic data set, including the study effect sizes, sampling variance, and the potential moderators.
#' @param c A non-negative scalar.The pruning parameter to prune the initial tree by the "c*standard-error" rule.
#' @param maxL the maximum number of splits
#' @param minsplit the minimum number of studies in a parent node before splitting
#' @param cp the stopping rule for the decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @param xval the number of folds to perform the cross-validation
#' @param minbucket the minimum number of the studies in a terminal node
#' @param lookahead an argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @param ... Additional arguments to be passed.
#' @return \strong{If (a) moderator effect(s) is(are) detected, the function will return a list including the following objects:}
#' @return tree: A data frame that represents the tree, with the Q-between and the residual heterogeneity (tau^2) after each split.
#' @return n: The number of the studies in each subgroup
#' @return moderators: the names of identified moderators
#' @return Qb: The between-subgroups Q-statistic
#' @return tau2: The estimate of the residual heterogeneity
#' @return df: The degrees of freedom of the between-subgroups Q test
#' @return pval.Qb: The p-value of the between-subgroups Q test
#' @return g: The subgroup summary effect size, based on Hedges'g
#' @return se: The standard error of subgroup summary effect size
#' @return zval: The test statistic of the subgroup summary effect size
#' @return pval: The p-value of the test statistic of the subgroup summary effect size
#' @return ci.lb: The lower bound of the confidence interval
#' @return ci.ub: The upper bound of the confidence interval
#' @return call: The matched call
#' @return cv.res: The cross-validation table
#' @return data: the data set subgrouped by the fitted tree
#' @return \strong{If no moderator effect is detected, the function will return a list including the following objects:}
#' @return n: The total number of the studies
#' @return Q: The Q-statistics for the heterogeneity test
#' @return df: The degree of freedoms of the heterogeneity test
#' @return pval.Q: The p-value for the heterogeneity test
#' @return g: The summary effect size for all studies (i.e., the overall effect size)
#' @return se: The standard error of the summary effect size
#' @return zval: The test statistic of the summary effect size
#' @return pval: The p-value for the test statistic of the summary effect size
#' @return ci.lb: The lower bound of the confidence interval for the summary effect size
#' @return ci.ub: The upper bound of the confidence interval for the summary effect size
#' @return call: The matched call
#' @examples data(dat.BCT2009)
#' library(Rcpp)
#' REtree <- REmrt(g ~ T1 + T2+ T4 +T25, vi = vi, data = dat.BCT2009, c = 0)
#' summary(REtree)
#' plot(REtree)
#' @importFrom stats terms model.response pchisq pnorm qnorm
#' @import Rcpp
#' @seealso \code{\link{summary.REmrt}}, \code{\link{plot.REmrt}}
#' @export

REmrt <- function(formula, data, vi, c = 1, maxL = 5, minsplit = 6, cp = 1e-5, minbucket = 3, xval = 10, lookahead = FALSE, ...){
  Call <- match.call()
  indx <- match(c("formula", "data", "vi"),
                names(Call), nomatch = 0L)
#------  Check if all arguments are correctly specified  -------#  
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  if (indx[3] == 0L)
    stop("The sampling variances need to be specified")
  if (!is.logical(lookahead))
    stop("The 'lookahead' argument needs to be a logical value")
  if (maxL < 2 & (lookahead == TRUE) )
    stop("The 'maxL' should be at least 2 when applying look-ahead strategy")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(temp)
  
#---    use x-validation to decide the size of the tree    ---#
  cv.res <- Xvalid_all(REmrt_GS_, mf, maxL = maxL, n.fold = xval, minbucket = minbucket, minsplit = minsplit, cp = cp, lookahead = lookahead)
  mindex <- which.min(cv.res[, 1])
  cp.minse <- cv.res[mindex,1] + c*cv.res[mindex,2]
  cp.row <- min(which(cv.res[,1]<= cp.minse))
  
  if (cp.row == 1) {  # if no tree was detected
    warning("no moderator effect was detected")
    y <- model.response(mf)
    vi <- c(t(mf["(vi)"]))
    wts <- 1/vi
    wy <- wts*y
    wy2 <- wts * y^2
    n <- length(y)
    Q <- sum(wy2) - (sum(wy))^2/sum(wts)
    df <- nrow(mf) - 1
    C <- sum(wts) - sum(wts^2)/sum(wts)
    tau2 <- max(0, (Q-df)/C)
    vi.star <- vi + tau2
    g <- sum(y/vi.star)/sum(1/vi.star)
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/vi.star))
    zval <- g/se
    pval <- pnorm(abs(zval), lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    res <- REmrt_GS_(mf, maxL = maxL, minsplit = minsplit, cp = cp, minbucket = minbucket, lookahead = lookahead)
    res.f <- list(n = n ,  Q = Q,
                df = df, pval.Q = pval.Q, tau2 = tau2, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, call = Call, data = mf, cv.res = cv.res, initial.tree = res$tree)
  } else{  
    y <- model.response(mf)
    vi <- c(t(mf["(vi)"]))
    res <- REmrt_GS_(mf, maxL = maxL, minsplit = minsplit, cp = cp, minbucket = minbucket, lookahead = lookahead)
    prunedTree <- res$tree[1:cp.row, ]
    tau2 <- res$tree$tau2[cp.row]
    vi.star <- vi + tau2
    subnodes <- res$node.split[[cp.row]]
    wy.star <- y/vi.star
    n <- tapply(y, subnodes, length)
    g <- tapply(wy.star, subnodes, sum)/tapply(1/vi.star, subnodes, sum)
    df <- cp.row - 1
    Qb <- res$tree$Qb[cp.row]
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(vi.star, subnodes, function(x) sqrt(1/sum(1/x)))
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    mod.names <- unique(prunedTree$mod[!is.na(prunedTree$mod)])
    mf$term.node <- subnodes
    res.f<- list(tree =  prunedTree, n = n, moderators =  mod.names, Qb = Qb, tau2 = tau2, df = df, pval.Qb = pval.Qb,
                 g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
                 ci.ub = ci.ub, call = Call, cv.res = cv.res, cpt = res$cpt,
                 data = mf, initial.tree = res$tree)
  }
  class(res.f) <- "REmrt"
  res.f
  
}
