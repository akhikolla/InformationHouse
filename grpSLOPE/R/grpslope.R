###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2016 Alexej Gossmann
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#' @importFrom stats lm pchisq qchisq qnorm sd uniroot coef predict
#' @importFrom utils head
#' @importFrom Rcpp sourceCpp
#' @useDynLib grpSLOPE, .registration = TRUE
NULL
#> NULL


#' Regularizing sequence for Group SLOPE
#'
#' Generate the regularizing sequence \code{lambda} for the Group SLOPE
#' problem according to one of multiple methods (see Details).
#'
#' Multiple methods are available to generate the regularizing sequence \code{lambda}:
#' \itemize{
#'   \item "max" -- lambdas as in Theorem 2.5 in Brzyski et. al. (2016).
#'     Provalby controls gFDR in orthogonal designs.
#'   \item "mean" -- lambdas of equation (2.16) in Brzyski et. al. (2016).
#'     Applicable for gFDR control in orthogonal designs. Less conservative than "max".
#'   \item "corrected" -- lambdas of Procedure 1 in Brzyski et. al. (2016);
#'     in the special case that all group sizes are equal and \code{wt} is a constant vector, 
#'     Procedure 6 of Brzyski et. al. (2016) is applied.
#'     Applicable for gFDR control when predictors from different groups are stochastically independent.
#' }
#'
#' @param method Possible values are "max", "mean",
#'    and "corrected". See under Details.
#' @param fdr Target group false discovery rate (gFDR)
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#' @param wt A named vector of weights, one weight per group of predictors
#'    (named according to names as in vector \code{group})
#' @param n.obs Number of observations (i.e., number of rows in \code{A});
#'    required only if method is "corrected"
#'
#' @examples
#' # specify 6 groups of sizes 2, 3, and 4
#' group <- c(1, 1, 2, 2, 2, 3, 3, 3, 3,
#'            4, 4, 5, 5, 5, 6, 6, 6, 6)
#' # set the weight for each group to the square root of the group's size
#' wt <- rep(c(sqrt(2), sqrt(3), sqrt(4)), 2)
#' names(wt) <- 1:6
#' # compute different lambda sequences
#' lambda.max <- lambdaGroupSLOPE(method="max", fdr=0.1, group=group, wt=wt) 
#' lambda.mean <- lambdaGroupSLOPE(method="mean", fdr=0.1, group=group, wt=wt) 
#' lambda.corrected <- lambdaGroupSLOPE(method="corrected", fdr=0.1,
#'                                      group=group, wt=wt, n.obs=1000)
#' rbind(lambda.max, lambda.mean, lambda.corrected)
#' #                      [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#' # lambda.max       2.023449 1.844234 1.730818 1.645615 1.576359 1.517427
#' # lambda.mean      1.880540 1.723559 1.626517 1.554561 1.496603 1.447609
#' # lambda.corrected 1.880540 1.729811 1.637290 1.568971 1.514028 1.467551
#'
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2016) \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1610.04960}
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2019) \emph{Group SLOPE -- adaptive selection of groups of predictors}. Journal of the American Statistical Association 114 (525): 419–33.
#' @export
lambdaGroupSLOPE <- function(method, fdr, group, wt, n.obs=NULL)
{
  # Prepare grouping information
  group.id    <- getGroupID(group)
  n.group     <- length(group.id)
  group.sizes <- sapply(group.id, FUN=length)

  # make sure weights are in the same order as group.id
  wt <- wt[names(group.id)]

  # compute the lambda sequence according to 'method'
  if (method=="max" | method=="mean") {
    lambda <- lambdaChiOrtho(fdr=fdr, n.group=n.group, wt=wt,
                             group.sizes=group.sizes, method=method)

  } else if (method=="corrected") {

    if (is.null(n.obs)) {
      stop("'n.obs' needs to be passed as an argument when method is 'corrected'")
    }

    # Check for equal group sizes and equal weights
    if ( (length(unique(group.sizes))==1) & (length(unique(wt))==1) ) {
      # lambdas of Procedure 6 in Brzyski et. al. (2016)
      m <- unique(group.sizes)
      w <- unique(wt)
      lambda <- lambdaChiEqual(fdr=fdr, n.obs=n.obs, 
                               n.group=n.group, m=m, w=w)
    } else {
      # lambdas of Procedure 1 in Brzyski et. al. (2016)
      lambda <- lambdaChiMean(fdr=fdr, n.obs=n.obs, n.group=n.group,
                              group.sizes=group.sizes, wt=wt)
    }

  } else {
    stop(paste(method, "is not a valid method."))
  }

  return(lambda)
}

#' Group SLOPE (Group Sorted L-One Penalized Estimation)
#' 
#' Performs selection of significant groups of predictors and estimation of the
#' corresponding coefficients using the Group SLOPE method (see Brzyski et. al., 2016).
#'
#' Multiple methods are available to generate the regularizing sequence \code{lambda},
#' see \code{\link{lambdaGroupSLOPE}} for detail.
#' The model matrix is transformed by orthogonalization within each group (see Section 2.1
#' in Brzyski et. al., 2016), and penalization is imposed on \eqn{\| X_{I_i} \beta_{I_i} \|}.
#' When \code{orthogonalize = TRUE}, due to within group orthogonalization,
#' the solution vector \code{beta} cannot be computed, if a group submatrix does not have full
#' column rank (e.g., if there are more predictors in a selected group than there are observations).
#' In that case only the solution vector \code{c} of the transformed (orthogonalized) model is returned.
#' Additionally, in any case the vector \code{group.norms} is returned with its \eqn{i}th entry
#' being \eqn{\| X_{I_i} \beta_{I_i} \|}, i.e., the overall effect of each group.
#' Note that all of these results are returned on the scale of the normalized versions of \code{X} and \code{y}.
#' However, \code{original.scale} contains the regression coefficients transformed to correspond to 
#' the original (unaltered) \code{X} and \code{y}. In that case, an estimate for the intercept term is also
#' returned with the other coefficients in \code{original.scale} (while on the normalized scale the estimate
#' of the intercept is always equal to zero, and is not explicitly provided in the \code{grpSLOPE} output).
#'
#' @param X The model matrix
#' @param y The response variable
#' @param group A vector describing the grouping structure. It should 
#'    contain a group id for each predictor variable.
#' @param fdr Target group false discovery rate (gFDR)
#' @param lambda Method used to obtain the regularizing sequence lambda. Possible
#'    values are "max", "mean", and "corrected" (default).
#'    See \code{\link{lambdaGroupSLOPE}} for detail. Alternatively, any
#'    non-increasing sequence of the correct length can be passed.
#' @param sigma Noise level. If ommited, estimated from the data, using Procedure 2 in Brzyski et. al. (2016).
#' @param verbose A \code{logical} specifying whether to print output or not
#' @param orthogonalize Whether to orthogonalize the model matrix within each group.
#'    Do not set manually unless you are certain that your data is appropriately pre-processed.
#' @param normalize Whether to center the input data and re-scale the columns
#'    of the design matrix to have unit norms. Do not disable this unless you
#'    are certain that your data are appropriately pre-processed.
#' @param max.iter See \code{\link{proximalGradientSolverGroupSLOPE}}.
#' @param dual.gap.tol See \code{\link{proximalGradientSolverGroupSLOPE}}.
#' @param infeas.tol See \code{\link{proximalGradientSolverGroupSLOPE}}.
#' @param x.init See \code{\link{proximalGradientSolverGroupSLOPE}}.
#' @param ... Options passed to \code{\link{prox_sorted_L1}}
#'
#' @return A list with members:
#'   \describe{
#'     \item{beta}{Solution vector. See Details.}
#'     \item{c}{Solution vector of the transformed model. See Details.}
#'     \item{group.norms}{Overall effect of each group. See Details.}
#'     \item{selected}{Names of selected groups (i.e., groups of predictors with at least one non-zero coefficient estimate)}
#'     \item{optimal}{Convergence status}
#'     \item{iter}{Iterations of the proximal gradient method}
#'     \item{lambda}{Regularizing sequence}
#'     \item{lambda.method}{Method used to construct the regularizing sequence}
#'     \item{sigma}{(Estimated) noise level}
#'     \item{group}{The provided grouping structure (corresponding to \code{beta})}
#'     \item{group.c}{Grouping structure of the transformed model (corresponding to \code{c})}
#'     \item{original.scale}{A list containing the estimated intercept and regression coefficients on the original scale. See Details.}
#'   }
#'
#' @examples
#' # generate some data
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(1:20), each=5)
#' b   <- c(runif(20), rep(0, 80))
#' # (i.e., groups 1, 2, 3, 4, are truly significant)
#' y   <- A %*% b + rnorm(10) 
#' fdr <- 0.1 # target false discovery rate
#' # fit a Group SLOPE model
#' result <- grpSLOPE(X=A, y=y, group=grp, fdr=fdr)
#' result$selected
#' # [1] "1"  "2"  "3"  "4"  "14"
#' result$sigma
#' # [1] 0.7968632
#' head(result$group.norms)
#' #         1         2         3         4         5         6 
#' #  2.905449  5.516103  8.964201 10.253792  0.000000  0.000000 
#'
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2016) \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1610.04960}
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2019) \emph{Group SLOPE -- adaptive selection of groups of predictors}. Journal of the American Statistical Association 114 (525): 419–33.
#'
#' @export
grpSLOPE <- function(X, y, group, fdr, lambda = "corrected", sigma = NULL,
                     verbose = FALSE, orthogonalize = NULL, normalize = TRUE,
                     max.iter=1e4, dual.gap.tol=1e-6, infeas.tol=1e-6,
                     x.init=NULL, ...) {

  # check inputs for NA or NaN
  if (anyNA(X)) { stop("Some entries of X are NA or NaN!") }
  if (anyNA(y)) { stop("Some entries of y are NA or NaN!") }
  if (anyNA(group)) { stop("Some entries of group are NA or NaN!") }
  
  # extract some additional info about the inputs
  group.id <- getGroupID(group)
  n.group  <- length(group.id)
  n <- nrow(X)

  # normalize X and y -------------------------------------------------
  # (save scaling values in order to obtain parameter estimates on the original scale later on)
  if (normalize) {
    # (1) center X, so that columns have means equal to zero
    X.mean <- apply(X, 2, mean)
    X <- X - matrix(rep(X.mean, each = n), nrow = n)
    # (2) scale X, so that columns have norms equal to one
    X.scaling <- apply(X, 2, function(col.of.X) { 1 / sqrt(sum(col.of.X^2)) })
    # check if division by 0 has occured
    if (!all(is.finite(X.scaling))) { 
      stop("X cannot be normalized (probably some column has sample variance equal to 0).")
    }
    X <- X %*% diag(X.scaling)
    # (3) center y
    y.mean <- mean(y)
    y <- y - y.mean
  }

  # within group orthogonalization ------------------------------------
  if (is.null(orthogonalize)) {
    if (is.numeric(lambda)) {
      stop("If lambda is numeric, the argument orthogonalize must be set manually.")
    }
    orthogonalize <- TRUE
  }

  if (orthogonalize) {
    ortho <- orthogonalizeGroups(X, group.id)
    # determine sizes of orthogonalized groups:
    ortho.group.length <- rep(NA, n.group)
    names(ortho.group.length) <- names(ortho)
    for (i in 1:n.group) {
      ortho.group.length[i] <- ncol(ortho[[i]]$Q)
    }
    # overwrite X with a matrix that contains orthogonalizations of the original blocks of X;
    # and set up grouping info for the orthogonalized version of X to be used in the optimization:
    X <- matrix(nrow = n, ncol = sum(ortho.group.length))
    grp <- rep(NA, sum(ortho.group.length))
    block.end <- cumsum(ortho.group.length)
    block.start <- head(c(1, block.end + 1), n.group)
    for (i in 1:n.group) {
      ind <- block.start[i]:block.end[i]
      grp[ind] <- names(ortho)[i]
      X[ , ind] <- ortho[[i]]$Q
    }
    ortho.group.id <- getGroupID(grp)
    # set prior weights per group:
    wt <- sqrt(ortho.group.length)
    wt.per.coef <- rep(NA, ncol(X))
    for (i in 1:n.group) { wt.per.coef[ ortho.group.id[[i]] ] <- wt[i] }
  } else {
    # grouping structure to be used in the optimization:
    grp <- group
    # set prior weights per group:
    wt <- sapply(group.id, length)
    wt <- sqrt(wt)
    wt.per.coef <- rep(NA, ncol(X))
    for (i in 1:n.group) { wt.per.coef[ group.id[[i]] ] <- wt[i] }
  }

  # regularizing sequence ---------------------------------------------
  if (is.character(lambda)) { 
    lambda.seq <- lambdaGroupSLOPE(method=lambda, fdr=fdr, group=grp, 
                                   wt=wt, n.obs=n)
  } else if (is.numeric(lambda)) {
    lambda.seq <- lambda
  } else {
    stop("Invalid input for lambda.")
  }


  # optimization ------------------------------------------------------
  if (!is.null(sigma)) {
    sigma.lambda <- sigma * lambda.seq
    optim.result <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=grp,
                                                     wt=wt.per.coef, 
                                                     lambda=sigma.lambda,
                                                     verbose=verbose,
                                                     max.iter=max.iter,
                                                     dual.gap.tol=dual.gap.tol, 
                                                     infeas.tol=infeas.tol,
                                                     x.init=x.init, ...)
  } else {
    # sigma needs to be estimated
    sigma <- sd(y)
    sigma.lambda <- sigma * lambda.seq
    optim.result <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=grp,
                                                     wt=wt.per.coef, 
                                                     lambda=sigma.lambda,
                                                     verbose=verbose,
                                                     max.iter=max.iter,
                                                     dual.gap.tol=dual.gap.tol, 
                                                     infeas.tol=infeas.tol,
                                                     x.init=x.init, ...)
    S.new <- which(optim.result$x != 0)
    S <- c()
    while( !isTRUE(all.equal(S, S.new)) && (length(S.new) > 0) ) {
      S <- S.new
      if (length(S) > n) {
        stop("Sigma estimation fails because more predictors got selected than there are observations.")
      }
      OLS <- lm(y ~ 0 + X[ , S])
      if (normalize) {
        sigma <- sqrt( sum(OLS$res^2) / (n - length(S) - 1) )
      } else {
        sigma <- sqrt( sum(OLS$res^2) / (n - length(S)) )
      }

      sigma.lambda <- sigma * lambda.seq
      optim.result <- proximalGradientSolverGroupSLOPE(y=y, A=X, group=grp,
                                                       wt=wt.per.coef, 
                                                       lambda=sigma.lambda,
                                                       verbose=verbose,
                                                       max.iter=max.iter,
                                                       dual.gap.tol=dual.gap.tol, 
                                                       infeas.tol=infeas.tol,
                                                       x.init=x.init, ...)
      S.new <- which(optim.result$x != 0)
    }
  }

  # create Group SLOPE solution object ---------------------------------
  sol <- list()
  sol$lambda <- lambda.seq
  sol$lambda.method <- lambda
  sol$sigma <- sigma
  sol$iter <- optim.result$iter

  # compute beta
  sol$c <- as.vector(optim.result$x)
  if (!orthogonalize) {
    sol$beta <- sol$c
  } else {
    # compute beta only if all groups have have the same number of columns
    # after orthogonalization as they did before
    group.length <- sapply(group.id, length)
    if (all(group.length == ortho.group.length)) {
      sol$beta <- rep(NA, length(sol$c))
      for (i in 1:n.group) {
        # c corresponds to the (reordered) group structure in orthogonalized version of X
        ci <- sol$c[ortho.group.id[[i]]]
        li <- ortho.group.length[i]

        bi <- tryCatch({ 
          backsolve(ortho[[i]]$R, ci) 
        }, error = function(err) {
          warning(paste("grpSLOPE caught an error:", err))
          return(rep(NA, li))
        })

        or <- rep(NA, li)
        for (j in 1:li) { or[j] <- which(ortho[[i]]$P == j) }
        betai <- bi[or]

        # beta corresponds to the group structure in the original matrix
        sol$beta[group.id[[i]]] <- betai
      }
    } else {
      sol$beta <- NULL
    }
  }

  # compute group norms ||X_I beta_I||
  sol$group.norms <- rep(NA, n.group)
  for (i in 1:n.group) {
    if (orthogonalize) {
      sol$group.norms[i] <- norm(as.matrix(sol$c[ortho.group.id[[i]]]), "f") 
    } else {
      Xbetai <- X[ , group.id[[i]]] %*% as.matrix(sol$beta[group.id[[i]]])
      sol$group.norms[i] <- norm(as.matrix(Xbetai), "f")
    }
  }
  group.names <- names(group.id)
  names(sol$group.norms) <- group.names 

  # selected groups
  sol$selected <- group.names[which(sol$group.norms != 0)]

  if (optim.result$status == 1) {
    sol$optimal <- TRUE
  } else {
    sol$optimal <- FALSE
  }

  # estimates on original scale
  if (is.null(sol$beta)) {
    sol$original.scale <- list("intercept" = NULL, "beta" = NULL)
  } else {
    if (normalize) {
      beta.original.scale <- diag(X.scaling) %*% sol$beta
      intercept <- y.mean - crossprod(X.mean, beta.original.scale)
    } else {
      intercept <- 0
      beta.original.scale <- sol$beta
    }
    sol$original.scale <- list("intercept" = as.double(intercept), 
                               "beta" = as.vector(beta.original.scale))
  }

  sol$group <- group
  sol$group.c <- grp

  class(sol) <- "grpSLOPE"
  return(sol)
}
