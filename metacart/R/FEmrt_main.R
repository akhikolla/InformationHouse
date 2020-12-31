#' Fixed effect meta-tree
#'
#' A function to fit fixed effect meta-trees to meta-analytic data.
#' The model is assuming a fixed effect within subgroups and between subgroups.
#' The tree growing process is equivalent to the approach described in Li et al. (2017) using fixed effect weights in the function \pkg{rpart()} developed by Therneau, Atkinson & Ripley (2014).
#' @name FEmrt
#' @aliases FEmrt
#' @inheritParams rpart
#' @param formula A formula, with an outcome variable (usually the effect size) and the potential moderator variables but no interaction terms.
#' @param vi sampling variance of the effect size.
#' @param data A data frame of a meta-analytic data set, including the study effect sizes, sampling variance, and the potential moderators.
#' @param subset optional expression that selects only a subset of the rows of the data.
#' @param c A non-negative scalar.The pruning parameter to prune the initial tree by the "c*standard-error" rule.
#' @param control the control object (similar to rpart.control from the rpart package) that is used in the tree algorithm
#' @param ... Additional arguments passed to the tree growing algorithm based on \pkg{rpart}.
#' @return \strong{If (a) moderator effect(s) is(are) detected, the function will return a \code{FEmrt} object including the following components:}
#' @return tree: The pruned tree that represents the moderator effect(s) and interaction effect(s) between them.
#' @return n: The number of the studies in each subgroup
#' @return Qb: The between-subgroups Q-statistic
#' @return df: The degree of freedoms of the between-subgroups Q test
#' @return pval.Qb: The p-value of the between-subgroups Q test
#' @return Qw: The within-subgroup Q-statistic in each subgroup
#' @return g: The subgroup summary effect size, based on Hedges'g
#' @return se: The standard error of the subgroup summary effect size
#' @return zval: The test statistic of the subgroup summary effect size
#' @return pval: The p-value for the test statistics of the subgroup summary effect size
#' @return ci.lb: The lower bound of the confidence interval
#' @return ci.ub: The upper bound of the confidence interval
#' @return call: The matched call
#' @return \strong{If no moderator effect is detected, the function will return a \code{FEmrt} object including the following components:}
#' @return n: The total number of the studies
#' @return Q: The Q-statistic of the heterogeneity test
#' @return df: The degrees of freedom of the heterogeneity test
#' @return pval.Q: The p-value of the heterogeneity test
#' @return g: The summary effect size for all studies
#' @return se: The standard error of the summary effect size
#' @return zval: The test statistic of the summary effect size
#' @return pval: The p-value of the test statistic of the summary effect size
#' @return ci.lb: The lower bound of the confidence interval for the summary effect size
#' @return ci.ub: The upper bound of the confidence interval for the summary effect size
#' @return call: The matched call
#' @examples data(dat.BCT2009)
#' library(Rcpp)
#' FEtree <- FEmrt(g ~ T1 + T2+ T4 + T25, vi = vi, data = dat.BCT2009, c = 0)
#' print(FEtree)
#' summary(FEtree)
#' plot(FEtree)
#' @references Dusseldorp, E., van Genugten, L., van Buuren, S., Verheijden, M. W., & van Empelen, P. (2014). Combinations of techniques that effectively change health behavior: Evidence from meta-cart analysis.  \emph{Health Psychology, 33(12)}, 1530-1540. doi:
#'      10.1037/hea0000018.
#' @references Li, X., Dusseldorp, E., & Meulman, J. J. (2017). Meta-CART: A tool to identify interactions between moderators in meta-analysis. \emph{ British Journal of Mathematical and Statistical Psychology, 70(1)}, 118-136. doi: 10.1111/bmsp.12088.
#'
#' @references Therneau, T., Atkinson, B., & Ripley, B. (2014) rpart: Recursive partitioning and regression trees. R package version, 4-1.
#' @seealso  \code{\link{summary.FEmrt}}, \code{\link{plot.FEmrt}}, \code{\link{rpart}},\code{\link{rpart.control}}
#' @importFrom rpart rpart.control
#' @importFrom rpart rpart
#' @importFrom rpart prune.rpart
#' @importFrom stats model.response
#' @importFrom stats pchisq
#' @importFrom stats pnorm
#' @importFrom stats qnorm
#' @import Rcpp
#' @export
FEmrt <- function(formula, data, vi, subset, c = 1,
                  control = rpart.control(xval = 10, minbucket = 3, minsplit = 6, cp = 0.0001),
                  ...) {
  Call <- match.call()
  wts.metacart <- NULL  # the weights to be used
  indx <- match(c("formula", "data", "vi", "subset"),
                names(Call), nomatch = 0L)
  if (indx[1] == 0L)
    stop("a 'formula' argument is required")
  if (indx[3] == 0L)
    stop("The sampling variances need to be specified")
  temp <- Call[c(1L, indx)]
  temp[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(temp)
  m$wts.metacart <- c(t(1/m["(vi)"]))
  tree <- rpart(formula, weights = wts.metacart, data = m, control = control, ...)
  prunedtree <- treepruner(tree, c*sqrt(mean(m$wts.metacart)))
  tree$cptable[ ,5] <- tree$cptable[ ,5]*sqrt(mean(m$wts.metacart))  # adjust the scale change due to the weights
  if (length(unique(prunedtree$where)) < 2) {
    warning("no moderator effect was detected")
    y <- model.response(m)
    v <- c(t(m["(vi)"]))
    n <- length(y)
    g <- sum(y/v)/sum(1/v)
    Q <- sum((y-g)^2/v)
    df <- n - 1
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/v))
    zval <- g/se
    pval <- pnorm(abs(zval), lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    res <- list(n = n ,  Q = Q,
                df = df, pval.Q = pval.Q, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, call = Call, cv.res = tree$cptable, data = m)
  } else {
    y <- model.response(m)
    v <- c(t(m["(vi)"]))
    treeframe <- prunedtree$frame
    n <- treeframe[treeframe$var == "<leaf>", 2]
    Qw <- treeframe[treeframe$var == "<leaf>", 4]
    g <- treeframe[treeframe$var == "<leaf>", 5]
    Qb <- treeframe[1,4] - sum(Qw)
    df <- length(unique(prunedtree$where))-1
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- tapply(v, prunedtree$where, function(x) sqrt(1/sum(1/x)))
    names(se) <- rownames(treeframe)[as.numeric(names(se))]
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    mod.names <- unique(prunedtree$frame$var[prunedtree$frame$var != "<leaf>"])
    res <- list(tree =  prunedtree, n = n, moderators =  mod.names, Qb = Qb, df = df, pval.Qb = pval.Qb,
                Qw = Qw, g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
                ci.ub = ci.ub, call = Call, cv.res = tree$cptable, data = m)

  }
  class(res) <- "FEmrt"
  res
}
