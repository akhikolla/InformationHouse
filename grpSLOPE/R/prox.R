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

#' Prox for sorted L1 norm
#'
#' Compute the prox for the sorted L1 norm. That is, given a vector \eqn{x}
#' and a decreasing vector \eqn{\lambda}, compute the unique value of \eqn{y}
#' minimizing
#' \deqn{\frac{1}{2} \Vert x - y \Vert_2^2 +
#'       \sum_{i=1}^n \lambda_i |x|_{(i)}.}
#' At present, two methods for computing the sorted L1 prox are
#' supported. By default, we use a fast custom C implementation. Since SLOPE
#' can be viewed as an isotonic regression problem, the prox can also be
#' computed using the \code{isotone} package. This option is provided
#' primarily for testing.
#'
#' @param x input vector
#' @param lambda vector of \eqn{\lambda}'s, sorted in decreasing order
#' @param method underlying prox implementation, either 'c' or 'isotone'
#'  (see Details)
#'
#' @details This function has been adapted (with only cosmetic changes) from
#' the R package \code{SLOPE} version 0.1.3, due to this function being
#' deprecated and defunct in \code{SLOPE} versions which are newer than 0.1.3.
#'
#' @export
prox_sorted_L1 <- function(x, lambda, method = c("c", "isotone")) {
  # Normalize input
  if (is.complex(x)) {
    sign_x <- complex(argument = Arg(x))
    x <- Mod(x)
  } else {
    sign_x <- sign(x)
    x <- abs(x)
  }

  # Sort input
  s <- sort(x, decreasing = TRUE, index.return = TRUE)

  # Compute prox
  impl <- switch(match.arg(method),
                 c = prox_sorted_L1_C,
                 isotone = prox_sorted_L1_isotone)
  result <- impl(s$x, lambda)

  # Restore original order and sign
  result[s$ix] <- result
  result * sign_x
}

prox_sorted_L1_isotone <- function(y, lambda) {
  loadNamespace("isotone")
  n <- length(y)

  # Solve the quadratic programming problem:
  #   min ||y - lambda - x||_2 s.t. x_1 >= x_2 >= ... >= x_n
  A_total_order <- cbind(2:n, 1:(n - 1))
  result <- isotone::activeSet(A_total_order,
                               y = y - lambda,
                               weights = rep(1, n))

  # Enforce non-negativity constraint.
  pmax(result$x, 0)
}
