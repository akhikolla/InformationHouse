#' Compute a common shared environment matrix
#'
#' Compute the common shared environment matrix for a set of related subjects.
#' The function is generic, and can accept a pedigree, or pedigreeList as the
#' first argument.
#'
#' When called with a pedigreeList, i.e., with multiple families, the routine
#' will create a block-diagonal-symmetric `bdsmatrix' object.  Since the [i,j]
#' value of the result is 0 for any two unrelated individuals i and j and a
#' `bdsmatix' utilizes sparse representation, the resulting object is often
#' orders of magnitude smaller than an ordinary matrix.  When called with a
#' single pedigree and ordinary matrix is returned.
#'
#' @param id either a pedigree object or pedigreeList object
#' @param rho The correlation between spouses
#' @param theta The partial path coefficient from parents to offspring
#' @param \dots Any number of optional arguments. Not used at the moment
#' @return a matrix of shared environment coefficients
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{pedigree}, \code{kinship},
#' @keywords datagen
#' @examples
#'
#' library(kinship2)
#' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
#'                     mom =c(0, 0, 0, 0, 0, 2, 2, 4, 0,  6,  8,  0, 10, 11),
#'                     dad =c(0, 0, 0, 0, 0, 1, 1, 3, 0,  5,  7,  0,  9, 12),
#'                     sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1,  2,  2,  1,  2,  2))
#'
#' tped <- with(test1, pedigree(id, dad, mom, sex))
#' extended.shared(tped)
#'
#' @export extended.shared
extended.shared <- function (id, rho=1, theta=1, ...)
{
  UseMethod("extended.shared")
}

#' @rdname extended.shared
#' @export
extended.shared.pedigreeList <- function (id, rho=1, theta=1, ...)
{
#    if (!require("Matrix")) {
#        stop("suggested package Matrix not installed")
#    }

    famlist <- unique(id$famid)
    nfam <- length(famlist)
    matlist <- vector("list", nfam)
    for (i in 1:length(famlist)) {
      family <- id[i]
      temp <- extended.shared(family)
      matlist[[i]] <- as(Matrix::forceSymmetric(temp), "dsCMatrix")
    }
    result <- Matrix::bdiag(matlist)

    if (any(duplicated(id$id))) {
      dimnames(result) <- list(NULL, paste(id$famid, id$id,
                                           sep = "/"))
    } else { dimnames(result) <- list(id$id, id$id) }
    result
}

#' @rdname extended.shared
#' @export
extended.shared.pedigree <- function (id, rho=1, theta=1, ...)
{
  n <- length(id$id)

  if (n == 1)
      return(matrix(1, 1, 1, dimnames = list(id$id, id$id)))
  if (any(duplicated(id$id)))
      stop("All id values must be unique")

  ## Create the I-B matrix
  IB <- diag(n)
  IB[cbind(seq(n),id$mindex)[id$mindex>0,]] <- -theta/2
  IB[cbind(seq(n),id$findex)[id$findex>0,]] <- -theta/2

  nonfounders <- (id$mindex+id$findex)>0

  ## Create the Gamma matrix
  ## and fix variances (standardization)
  GAMmat <- diag(n)
  diag(GAMmat)[nonfounders] <- sqrt(1-theta^2/2*(1+rho))

  ## Spouses
  spouses <- cbind(id$mindex, id$findex)[nonfounders,]
  U <- diag(n)
  U[spouses]       <- rho
  U[spouses[,2:1]] <- rho

  temp <- solve(IB) %*% GAMmat %*% U  %*%  t(GAMmat) %*% t(solve(IB))
  dimnames(temp) <- list(id$id, id$id)
  temp
}
