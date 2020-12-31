###############################################################################
#
#    grpSLOPE: Group SLOPE (Group Sorted L1 Penalized Estimation)
#    Copyright (C) 2018 Alexej Gossmann
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

#' @importFrom stats rnorm
NULL
#> NULL

#' Prox for group SLOPE
#'
#' Evaluate the proximal mapping for the group SLOPE problem.
#'
#' \code{proxGroupSortedL1} evaluates the proximal mapping of the group SLOPE
#' problem by reducing it to the prox for the (regular) SLOPE and then applying
#' the fast prox algorithm for the Sorted L1 norm.
#'
#' @param y The response vector
#' @param group Either a vector or an object of class \code{groupID} (e.g.,
#'   as produced by \code{\link{getGroupID}}), which is describing the
#'   grouping structure. If it is a vector, then it should contain
#'   a group id for each predictor variable.
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}
#' @param ... Options passed to \code{\link{prox_sorted_L1}}
#'
#' @examples
#' grp <- c(0,0,0,1,1,0,2,1,0,2)
#' proxGroupSortedL1(y = 1:10, group = grp, lambda = c(10, 9, 8))
#' #  [1] 0.2032270 0.4064540 0.6096810 0.8771198 1.0963997 1.2193620 1.3338960
#' #  [8] 1.7542395 1.8290430 1.9055657
#'
#' @references M. Bogdan, E. van den Berg, C. Sabatti, W. Su, E. Candes (2015), \emph{SLOPE -- Adaptive variable selection via convex optimization}, \url{http://arxiv.org/abs/1407.3824}
#'
#' @export
proxGroupSortedL1 <- function(y, group, lambda, ...) {
  if (inherits(group, "groupID")) {
    n.group <- length(group)
    group.id <- group
  } else {
    n.group <- length(unique(group))
    group.id <- getGroupID(group)
  }

  if (length(lambda) != n.group) {
    stop("Length of lambda should be equal to the number of groups.")
  }

  # compute Euclidean norms for groups in y
  group.norm <- rep(NA, n.group)
  for (i in 1:n.group){
    selected <- group.id[[i]]
    group.norm[i] <- norm(as.matrix(y[selected]), "f")
  }

  # get Euclidean norms of the solution vector
  prox.norm <- prox_sorted_L1(group.norm, lambda, ...)

  # compute the solution
  prox.solution <- rep(NA, length(y))
  for (i in 1:n.group){
    selected <- group.id[[i]]
    prox.solution[selected] <- prox.norm[i] / group.norm[i] * y[selected]
  }

  return(prox.solution)
}

# Compute duality gap for Group SLOPE (eq. (B.20) in
# Appendix B in Brzyski et al. (2016))
#
getDualityGapGroupSLOPE <- function(group.id, n.group, b, lambda, g) {
  b.norms <- rep(NA, n.group)
  for (i in 1:n.group) {
    selected <- group.id[[i]]
    b.norms[i] <- norm(as.matrix(b[selected]), "f")
  }
  b.norms.sorted <- sort(b.norms, decreasing = TRUE)
  duality.gap <- crossprod(b, g) + crossprod(lambda, b.norms.sorted)
  return(duality.gap)
}

# Compute the infeasibility for Group SLOPE.
# A comment on the derivation of this infeasibility criterion can be found at:
# https://github.com/agisga/grpSLOPE/blob/gh-pages/_posts/outdated_2015-12-4-Infeasibility.md
#
getInfeasibilityGroupSLOPE <- function(group.id, n.group, lambda, g) {
  g.norms <- rep(NA, n.group)
  for (i in 1:n.group) {
    selected <- group.id[[i]]
    g.norms[i] <- norm(as.matrix(g[selected]), "f")
  }
  g.norms.sorted  <- sort(g.norms, decreasing = TRUE)
  infeasibility <- max(max(cumsum(g.norms.sorted - lambda)), 0)
  return(infeasibility)
}

#' Proximal gradient method for Group SLOPE
#'
#' Compute the coefficient estimates for the Group SLOPE problem.
#'
#' \code{proximalGradientSolverGroupSLOPE} computes the coefficient estimates
#' for the Group SLOPE model. The employed optimization algorithm is FISTA with
#' backtracking Lipschitz search.
#'
#' @param y the response vector
#' @param A the model matrix
#' @param group A vector describing the grouping structure. It should
#'   contain a group id for each predictor variable.
#' @param wt A vector of weights (per coefficient)
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}
#' @param max.iter Maximal number of iterations to carry out
#' @param verbose A \code{logical} specifying whether to print output or not
#' @param dual.gap.tol The tolerance used in the stopping criteria for the duality gap
#' @param infeas.tol The tolerance used in the stopping criteria for the infeasibility
#' @param x.init An optional initial value for the iterative algorithm
#' @param ... Options passed to \code{\link{prox_sorted_L1}}
#'
#' @return A list with the entries:
#'   \describe{
#'     \item{x}{Solution (n-by-1 matrix)}
#'     \item{status}{Convergence status: 1 if optimal, 2 if iteration limit reached}
#'     \item{L}{Approximation of the Lipschitz constant (step size)}
#'     \item{iter}{Iterations of the proximal gradient method}
#'     \item{L.iter}{Total number of iterations spent in Lipschitz search}
#'   }
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(runif(100, 0, 1), 10, 10)
#' grp <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
#' wt  <- c(2, 2, 2, 2, 5, 5, 5, 5, 5, 1)
#' x   <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
#' y   <- A %*% x
#' lam <- 0.1 * (10:7)
#' result <- proximalGradientSolverGroupSLOPE(y=y, A=A, group=grp, wt=wt, lambda=lam, verbose=FALSE)
#' result$x
#' #           [,1]
#' #  [1,] 0.000000
#' #  [2,] 0.000000
#' #  [3,] 3.856005
#' #  [4,] 2.080736
#' #  [5,] 0.000000
#' #  [6,] 0.000000
#' #  [7,] 0.000000
#' #  [8,] 0.000000
#' #  [9,] 0.000000
#' # [10,] 3.512833
#'
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2016) \emph{Group SLOPE -- adaptive selection of groups of predictors}, \url{https://arxiv.org/abs/1610.04960}
#' @references D. Brzyski, A. Gossmann, W. Su, and M. Bogdan (2019) \emph{Group SLOPE -- adaptive selection of groups of predictors}. Journal of the American Statistical Association 114 (525): 419–33. \url{http://dx.doi.org/10.1080/01621459.2017.1411269}
#' @references A. Gossmann, S. Cao, Y.-P. Wang (2015) \emph{Identification of Significant Genetic Variants via SLOPE, and Its Extension to Group SLOPE}. In Proceedings of ACM BCB 2015. \url{http://dx.doi.org/10.1145/2808719.2808743}
#'
#' @export
proximalGradientSolverGroupSLOPE <- function(y, A, group, wt, lambda,
                                             max.iter = 1e4, verbose = FALSE,
                                             dual.gap.tol = 1e-6,
                                             infeas.tol = 1e-6,
                                             x.init = NULL, ...) {
  # This function is based on the source code available from
  # http://statweb.stanford.edu/~candes/SortedL1/software.html
  # under the GNU GPL-3 licence.
  #
  # Original implementation: Copyright 2013, M. Bogdan, E. van den Berg, W. Su, and E.J. Candes
  # Modifications: Copyright 2015, Alexej Gossmann

  group.id <- getGroupID(group)
  n.group  <- length(group.id)
  p <- ncol(A)
  if (is.null(x.init)) x.init <- matrix(0, p, 1)

  # Ensure that lambda is non-increasing and of right length
  n.lambda <- length(lambda)
  if ( (n.lambda > 1) & any(lambda[2:n.lambda] > lambda[1:(n.lambda - 1)]) ) {
    stop("Lambda must be non-increasing.")
  }
  if (n.lambda != n.group) {
    stop("Lambda must have exactly as many entries as there are groups.")
  }

  # Adjust matrix for prior weights
  Dinv <- diag(1 / wt)
  A    <- A %*% Dinv

  # Set constants
  STATUS_RUNNING    <- 0
  STATUS_OPTIMAL    <- 1
  STATUS_ITERATIONS <- 2
  STATUS_MSG        <- c('Optimal','Iteration limit reached')

  # Auxilliary function to record output information
  recordResult <- function(b, status, L, iter, L.iter) {
    result           <- list()
    result$x         <- Dinv %*% b
    result$status    <- status
    result$L         <- L
    result$iter      <- iter
    result$L.iter    <- L.iter
    return(result)
  }

  # Run regular SLOPE if all groups are singletons
  if (length(group.id) == p) {
    sol <- SLOPE_solver(A = A, b = y, lambda = lambda,
                        initial = x.init, max_iter = max.iter,
                        tol_infeas = infeas.tol,
                        tol_rel_gap = dual.gap.tol)
    status <- STATUS_ITERATIONS
    if (sol$optimal) status <- STATUS_OPTIMAL
    result <- recordResult(b = matrix(sol$x, c(p, 1)), status = status,
                           L = sol$lipschitz, iter = sol$iter, L.iter = NULL)
    return(result)
  }

  # Get initial lower bound on the Lipschitz constant
  rand.mat <- matrix(rnorm(p), c(p, 1))
  rand.mat <- rand.mat / norm(rand.mat, "f")
  rand.mat <- t(A) %*% (A %*% rand.mat)
  L <- norm(rand.mat, "f")

  # Initialize parameters and iterates
  tt       <- 1
  eta      <- 2
  lambda   <- as.matrix(lambda)
  y        <- as.matrix(y)
  x        <- x.init
  b        <- x
  Ax       <- A %*% x
  f        <- Inf
  f.prev   <- Inf
  iter     <- 0
  L.iter   <- 0
  status   <- STATUS_RUNNING
  duality.gap <- Inf
  infeasibility <- Inf

  if (verbose) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf("%5s  %9s  %9s  %9s\n", "Iter", "||r||_2", "Dual. gap", "Infeas.")
  }

  # Main loop
  while (TRUE) {

    # Compute the gradient
    r <- (A %*% b) - y
    g <- crossprod(A, r)

    # Stopping criteria

    duality.gap <- getDualityGapGroupSLOPE(group.id, n.group, b, lambda, g)
    infeasibility <- getInfeasibilityGroupSLOPE(group.id, n.group, lambda, g)

    if (verbose) {
      out_str <- sprintf("   %9.2e  %9.2e", duality.gap, infeasibility)
      printf("%5d  %9.2e%s\n", iter, f, out_str)
    }

    if ( (duality.gap < dual.gap.tol) & (infeasibility < infeas.tol) ) {
      status <- STATUS_OPTIMAL
    } else if ( (status == STATUS_RUNNING) & (iter >= max.iter) ) {
      status <- STATUS_ITERATIONS
    }

    if (status != STATUS_RUNNING) {
      if (verbose) {
        printf("Exiting with status %d -- %s\n", status, STATUS_MSG[[status]])
      }
      break
    }

    # Increment iteration count
    iter <- iter + 1

    # Keep copies of previous values
    f.prev  <- as.double(crossprod(r)) / 2
    Ax.prev <- Ax
    x.prev  <- x
    b.prev  <- b
    tt.prev <- tt

    # Lipschitz search
    while (TRUE) {
      x  <- proxGroupSortedL1(y = b - (1 / L) * g, group = group.id,
                              lambda = lambda / L, ...)
      d  <- x - b
      Ax <- A %*% x
      r  <- Ax - y
      f  <- as.double(crossprod(r)) / 2
      qq <- f.prev + as.double(crossprod(d, g)) + (L / 2) * as.double(crossprod(d))

      L.iter <- L.iter + 1

      if (qq >= f * (1 - 1e-12))
        break
      else
        L <- L * eta
    }

    # Update
    tt <- (1 + sqrt(1 + 4 * tt^2)) / 2
    b  <- x + ( (tt.prev - 1) / tt ) * (x - x.prev)
  }

  result <- recordResult(b = b, status = status, L = L,
                         iter = iter, L.iter = L.iter)
  return(result)
}

#' Alternating direction method of multipliers
#'
#' Compute the coefficient estimates for the Group SLOPE problem.
#'
#' \code{admmSolverGroupSLOPE} computes the coefficient estimates
#' for the Group SLOPE model. The employed optimization algorithm is
#' the alternating direction method of multipliers (ADMM).
#'
#' @param y the response vector
#' @param A the model matrix
#' @param group A vector describing the grouping structure. It should
#'   contain a group id for each predictor variable.
#' @param wt A vector of weights (per coefficient)
#' @param lambda A decreasing sequence of regularization parameters \eqn{\lambda}
#' @param rho Penalty parameter in the augmented Lagrangian (see Boyd et al., 2011)
#' @param max.iter Maximal number of iterations to carry out
#' @param verbose A \code{logical} specifying whether to print output or not
#' @param absolute.tol The absolute tolerance used in the stopping criteria for the primal and dual feasibility conditions (see Boyd et al., 2011, Sec. 3.3.1)
#' @param relative.tol The relative tolerance used in the stopping criteria for the primal and dual feasibility conditions (see Boyd et al., 2011, Sec. 3.3.1)
#' @param z.init An optional initial value for the iterative algorithm
#' @param u.init An optional initial value for the iterative algorithm
#' @param ... Options passed to \code{\link{prox_sorted_L1}}
#'
#' @return A list with the entries:
#'   \describe{
#'     \item{x}{Solution (n-by-1 matrix)}
#'     \item{status}{Convergence status: 1 if optimal, 2 if iteration limit reached}
#'     \item{iter}{Number of iterations of the ADMM method}
#'   }
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(runif(100, 0, 1), 10, 10)
#' grp <- c(0, 0, 1, 1, 2, 2, 2, 2, 2, 3)
#' wt  <- c(2, 2, 2, 2, 5, 5, 5, 5, 5, 1)
#' x   <- c(0, 0, 5, 1, 0, 0, 0, 1, 0, 3)
#' y   <- A %*% x
#' lam <- 0.1 * (10:7)
#' result <- admmSolverGroupSLOPE(y = y, A = A, group = grp, wt = wt,
#'                                lambda=lam, rho = 1, verbose = FALSE)
#' result$x
#' #           [,1]
#' #  [1,] 0.000000
#' #  [2,] 0.000000
#' #  [3,] 3.856002
#' #  [4,] 2.080742
#' #  [5,] 0.000000
#' #  [6,] 0.000000
#' #  [7,] 0.000000
#' #  [8,] 0.000000
#' #  [9,] 0.000000
#' # [10,] 3.512829
#'
#' @references S. Boyd, N. Parikh, E. Chu, B. Peleato, and J. Eckstein (2011) \emph{Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers.} Foundations and Trends in Machine Learning 3 (1).
#'
#' @export
admmSolverGroupSLOPE <- function(y, A, group, wt, lambda, rho = NULL,
                                 max.iter = 1e4, verbose = FALSE,
                                 absolute.tol = 1e-4, relative.tol = 1e-4,
                                 z.init = NULL, u.init = NULL, ...) {

  group.id <- getGroupID(group)
  n.group  <- length(group.id)
  p <- ncol(A)
  n.sample <- nrow(A)

  # Check function arguments
  n.lambda <- length(lambda)
  if ( (n.lambda > 1) & any(lambda[2:n.lambda] > lambda[1:(n.lambda - 1)]) ) {
    stop("Lambda must be non-increasing.")
  }
  if (n.lambda != n.group) {
    stop("Lambda must have exactly as many entries as there are groups.")
  }
  if (is.null(rho)) {
    stop("The argument rho needs to be explicitly specified.")
  }
  if (is.null(z.init)) z.init <- matrix(0, p, 1)
  if (is.null(u.init)) u.init <- matrix(0, p, 1)

  # Adjust matrix for prior weights
  Dinv <- diag(1 / wt)
  A    <- A %*% Dinv

  # Set constants
  STATUS_RUNNING    <- 0
  STATUS_OPTIMAL    <- 1
  STATUS_ITERATIONS <- 2
  STATUS_MSG        <- c('Optimal','Iteration limit reached')

  # Auxilliary function to record output information
  recordResult <- function(b, status, iter) {
    result           <- list()
    result$x         <- Dinv %*% b
    result$status    <- status
    result$iter      <- iter
    return(result)
  }

  # TODO: Run ADMM version of regular SLOPE if all groups are singletons

  # cache the vector (rho I + A^T A)^{-1} A^T y = M A^T y, and the matrix M
  if (n.sample >= p) {
    M <- solve(diag(rho, p) + crossprod(A))
    MAty <- tcrossprod(M, crossprod(y, A))
  } else {
    Woodbury_trick <- solve(diag(rho, n.sample) + tcrossprod(A))
    MAt.transposed <- crossprod(Woodbury_trick, A)
    MAty <- crossprod(MAt.transposed, y)
    M <- (diag(p) - crossprod(MAt.transposed, A)) / rho
  }
  lambda_rho <- lambda / rho
  # Initialize parameters and iterates
  z <- z.init
  u <- u.init
  x <- NULL
  z.prev <- NULL
  status <- STATUS_RUNNING
  dual.feasibility <- Inf
  primal.feasibility <- Inf
  iter <- 0

  if (verbose) {
    printf <- function(...) invisible(cat(sprintf(...)))
    printf("%5s  %11s  %11s\n", "Iter", "Prim. feas.", "Dual feas.")
  }

  # Main loop
  while (TRUE) {

    iter <- iter + 1

    # ADMM updates

    x <- MAty + crossprod(M, (rho * (z - u)))
    z.prev <- z
    z <- proxGroupSortedL1(y = as.vector(x + u), group = group.id,
                           lambda = lambda_rho, ...)
    u <- u + x - z

    # Stopping criteria

    primal.feasibility <- norm(as.matrix(z - x), type = "f")
    dual.feasibility <- norm(as.matrix(rho*(z - z.prev)), type = "f")

    primal.tol <- sqrt(p) * absolute.tol +
      relative.tol * max(norm(as.matrix(x), type = "f"), norm(as.matrix(z), type = "f"))
    dual.tol <- sqrt(p) * absolute.tol +
      relative.tol * rho * norm(as.matrix(u), type = "f")

    if (verbose) {
      printf("%5s  %11.2e  %11.2e\n", iter,
             dual.feasibility, primal.feasibility)
    }

    if ( (dual.feasibility < dual.tol) & (primal.feasibility < primal.tol) ) {
      status <- STATUS_OPTIMAL
    } else if ( (status == STATUS_RUNNING) & (iter >= max.iter) ) {
      status <- STATUS_ITERATIONS
    }

    if (status != STATUS_RUNNING) {
      if (verbose) {
        printf("Exiting with status %d -- %s\n", status, STATUS_MSG[[status]])
      }
      break
    }
  }

  result <- recordResult(b = z, status = status, iter = iter)
  return(result)
}
