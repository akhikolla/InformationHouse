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

#' Extract model coefficients
#' 
#' Extract the regression coefficients from a \code{grpSLOPE} object, either on the
#' scale of the normalized design matrix (i.e., columns centered and scaled to unit norm),
#' or on the original scale.
#'
#' If \code{scaled} is set to \code{TRUE}, then the coefficients are returned for the 
#' normalized version of the design matrix, which is the scale on which they were computed. 
#' If \code{scaled} is set to \code{FALSE}, then the coefficients are transformed to
#' correspond to the original (unaltered) design matrix.
#' In case that \code{scaled = FALSE}, an estimate for the intercept term is returned with
#' the other coefficients. In case that \code{scaled = TRUE}, the estimate of the intercept 
#' is always equal to zero, and is not explicitly provided.
#'
#' @param object A \code{grpSLOPE} object 
#' @param scaled Should the coefficients be returned for the normalized version of the design matrix?
#' @param ... Potentially further arguments passed to and from methods
#'
#' @return A named vector of regression coefficients where the names signify the group that each entry belongs to
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(letters[1:20]), each=5)
#' b   <- c(rep(1, 20), rep(0, 80))
#' y   <- A %*% b + rnorm(10) 
#' result <- grpSLOPE(X=A, y=y, group=grp, fdr=0.1)
#' head(coef(result), 8)
#' #       a_1       a_2       a_3       a_4       a_5       b_1       b_2       b_3 
#' #  7.942177  7.979269  8.667013  8.514861 10.026664  8.963364 10.037355 10.448692 
#' head(coef(result, scaled = FALSE), 8)
#' # (Intercept)         a_1         a_2         a_3         a_4         a_5         b_1         b_2 
#' #  -0.4418113   0.8886878   0.8372108   0.8422089   0.8629597   0.8615827   0.9323849   0.9333445 
#' 
#' @export
coef.grpSLOPE <- function(object, scaled = TRUE, ...) {
  if(is.null(object$beta)) {
    stop("beta is set to NULL. Maybe one of the group submatrices did not have full column rank? See documentation to grpSLOPE().")
  }

  coef.names <- rep(NA, length(object$beta))
  group.id <- getGroupID(object$group)
  group.length <- sapply(group.id, length)
  n.group <- length(group.id)
  for (i in 1:n.group) {
    coef.names[ group.id[[i]] ] <- paste0(names(group.id)[i], "_", 1:group.length[i])
  }

  if(scaled) {
    coefs <- object$beta
    names(coefs) <- coef.names
  } else {
    coefs <- c(object$original.scale$intercept, object$original.scale$beta)
    names(coefs) <- c("(Intercept)", coef.names)
  }

  return(coefs)
}

if (getRversion() < "3.3.0") {
  sigma <- function(object, ...) UseMethod("sigma", object)
}

#' Extract (estimated) noise level
#' 
#' Extract the noise level of the \code{grpSLOPE} model.
#'
#' This basically obtains \code{object$sigma}. For \code{R (>= 3.3.0)}
#' \code{sigma} is an S3 method with the default method coming from the
#' \code{stats} package.
#'
#' @param object A \code{grpSLOPE} object 
#' @param ... Potentially further arguments passed to and from methods
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(1:20), each = 5)
#' b   <- c(rep(1, 20), rep(0, 80))
#' y   <- A %*% b + rnorm(10) 
#' # model with unknown noise level
#' result <- grpSLOPE(X = A, y = y, group = grp, fdr = 0.1)
#' sigma(result)
#' # [1] 0.6505558
#' # model with known noise level
#' result <- grpSLOPE(X = A, y = y, group = grp, fdr = 0.1, sigma = 1)
#' sigma(result)
#' # [1] 1
#' 
#' @rawNamespace if (getRversion() >= "3.3.0") importFrom(stats,sigma)
#' @rawNamespace if (getRversion() < "3.3.0") export(sigma)
#' @name sigma
#' @export
sigma.grpSLOPE <- function(object, ...) {
  return(object$sigma)
}

#' Obtain predictions
#' 
#' Obtain predictions from a \code{grpSLOPE} model on new data
#'
#' Note that \code{newdata} must have the same shape, and must contain
#' the same predictor variables as columns in the same order as the
#' design matrix \code{X} that was used for the \code{grpSLOPE} model fit.
#'
#' @param object A \code{grpSLOPE} object 
#' @param newdata Predictor variables arranged in a matrix 
#' @param ... Potentially further arguments passed to and from methods
#'
#' @examples
#' set.seed(1)
#' A   <- matrix(rnorm(100^2), 100, 100)
#' grp <- rep(rep(1:20), each = 5)
#' b   <- c(rep(1, 20), rep(0, 80))
#' y   <- A %*% b + rnorm(10) 
#' result <- grpSLOPE(X = A, y = y, group = grp, fdr = 0.1)
#' newdata <- matrix(rnorm(800), 8, 100)
#' # group SLOPE predictions:
#' predict(result, newdata)
#' # [1] -5.283385 -6.313938 -3.173068  1.901488  9.796677 -0.144516 -0.611164 -5.167620
#' # true mean values:
#' as.vector(newdata %*% b)
#' # [1] -5.0937160 -6.5814111 -3.5776124  2.7877449 11.0668777  1.0253236 -0.4261076 -4.8622940
#' 
#' @export
predict.grpSLOPE <- function(object, newdata, ...) {
  intercept <- object$original.scale$intercept
  slopes <- object$original.scale$beta
  predictions <- intercept + newdata %*% slopes
  return(as.vector(predictions))
}
