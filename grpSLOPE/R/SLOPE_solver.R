# SLOPE: Sorted L1 Penalized Estimation (SLOPE)
# Copyright (C) 2015 Malgorzata Bogdan, Ewout van den Berg, Chiara Sabatti,
# Weijie Su, Emmanuel Candes, Evan Patterson
# Copyright (C) 2020 Alexej Gossmann (minor modifications to the original code)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#' Sorted L1 solver
#'
#' Solves the sorted L1 penalized regression problem: given a matrix \eqn{A},
#' a vector \eqn{b}, and a decreasing vector \eqn{\lambda}, find the vector
#' \eqn{x} minimizing
#' \deqn{\frac{1}{2}\Vert Ax - b \Vert_2^2 +
#'       \sum_{i=1}^p \lambda_i |x|_{(i)}.}
#' This optimization problem is convex and is solved using an accelerated
#' proximal gradient descent method.
#'
#' @param A an \eqn{n}-by-\eqn{p} matrix
#' @param b vector of length \eqn{n}
#' @param lambda vector of length \eqn{p}, sorted in decreasing order
#' @param initial initial guess for \eqn{x}
#' @param prox function that computes the sorted L1 prox
#' @param max_iter maximum number of iterations in the gradient descent
#' @param grad_iter number of iterations between gradient updates
#' @param opt_iter number of iterations between checks for optimality
#' @param tol_infeas tolerance for infeasibility
#' @param tol_rel_gap tolerance for relative gap between primal and dual
#'  problems
#'
#' @return An object of class \code{SLOPE_solver.result}. This object is a list
#'  containing at least the following components:
#'  \item{x}{solution vector \eqn{x}}
#'  \item{optimal}{logical: whether the solution is optimal}
#'  \item{iter}{number of iterations}
#'
#' @details This function has been adapted (with only cosmetic changes) from
#' the R package \code{SLOPE} version 0.1.3, due to this function being
#' deprecated and defunct in \code{SLOPE} versions which are newer than 0.1.3.
#'
#' @export
# Adapted from Adlas.R in original SLOPE distribution.
SLOPE_solver <- function(A, b, lambda, initial = NULL, prox = prox_sorted_L1,
                         max_iter = 10000, grad_iter = 20, opt_iter = 1,
                         tol_infeas = 1e-6, tol_rel_gap = 1e-6) {
  # Get problem dimension.
  n <- ncol(A)

  # Get initial lower bound on the Lipschitz constant.
  x <- with_seed(0, rnorm(n))
  x <- x / sqrt(sum(x ^ 2))
  x <- t(A) %*% (A %*% x)
  L <- sqrt(sum(x ^ 2))

  # Initialize parameters and iterates.
  x.init <- if (is.null(initial)) rep(0, n) else initial
  t      <- 1
  eta    <- 2
  x      <- x.init
  y      <- x
  Ax     <- A %*% x
  f.prev <- Inf
  iter   <- 0
  optimal <- FALSE

  # Main loop.
  repeat {
    # Compute the gradient at f(y).
    if ( (iter %% grad_iter) == 0) # Includes first iterations
      r <- (A %*% y) - b
    else
      r <- (Ax + ( (t.prev - 1) / t) * (Ax - Ax.prev)) - b
    g <- t(A) %*% r
    f <- as.double(crossprod(r)) / 2

    # Increment iteration count.
    iter <- iter + 1

    # Check optimality conditions.
    if ( (iter %% opt_iter) == 0) {
      # Compute 'dual', check infeasibility and gap.
      gs     <- sort(abs(g), decreasing = TRUE)
      ys     <- sort(abs(y), decreasing = TRUE)
      infeas <- max(max(cumsum(gs - lambda)), 0)

      # Compute primal and dual objective.
      obj_primal <-  f + as.double(crossprod(lambda, ys))
      obj_dual   <- (-f) - as.double(crossprod(r, b))

      # Check primal-dual gap.
      if ( (abs(obj_primal - obj_dual) / max(1, obj_primal) < tol_rel_gap) &&
            (infeas < tol_infeas * lambda[[1]]))
        optimal <- TRUE;
    }

    # Stopping criteria.
    if (optimal || (iter >= max_iter))
      break;

    # Store copies of previous values.
    Ax.prev <- Ax
    x.prev <- x
    f.prev <- f
    t.prev <- t

    # Lipschitz search.
    repeat {
      # Compute prox mapping.
      x <- prox(y - (1 / L) * g, lambda / L)
      d <- x - y

      Ax <- A %*% x
      r  <- Ax - b
      f  <- as.double(crossprod(r)) / 2
      q  <- (f.prev +
             as.double(crossprod(d, g)) +
             (L / 2) * as.double(crossprod(d)))

      if (q < f * (1 - 1e-12))
        L <- L * eta
      else
        break
    }

    # Update.
    t <- (1 + sqrt(1 + 4 * t ^ 2)) / 2
    y <- x + ( (t.prev - 1) / t) * (x - x.prev)
  }
  if (!optimal)
    warning("SLOPE solver reached iteration limit")

  # Package up the results.
  structure(list(x = as.vector(y),
                 optimal = optimal,
                 iter = iter,
                 infeas = infeas,
                 obj_primal = obj_primal,
                 obj_dual = obj_dual,
                 lipschitz = L),
            class = "SLOPE_solver.result")
}
