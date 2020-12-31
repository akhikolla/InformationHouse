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
#' @param \dots Any number of optional arguments. Not used at the moment
#' @return a matrix of shared environment coefficients
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @seealso \code{pedigree}, \code{kinship},
#' @keywords datagen
#' @examples
##'
##' library(kinship2)
##' test1 <- data.frame(id  =c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14),
##'                     mom =c(0, 0, 0, 0, 2, 2, 4, 4, 6,  2,  0,  0, 12, 13),
##'                     dad =c(0, 0, 0, 0, 1, 1, 3, 3, 3,  7,  0,  0, 11, 10),
##'                     sex =c(1, 2, 1, 2, 1, 2, 1, 2, 1,  1,  1,  2,  2,  2))
##' tped <- with(test1, pedigree(id, dad, mom, sex))
##' common.shared(tped)
#'
#' @importFrom methods as
#' @rdname common.shared
#' @export
common.shared <- function (id, ...)
{
  UseMethod("common.shared")
}


#' @rdname common.shared
#' @export
common.shared.pedigreeList <- function (id, ...)
{

    famlist <- unique(id$famid)
    nfam <- length(famlist)
    matlist <- vector("list", nfam)
    for (i in 1:length(famlist)) {
      family <- id[i]
      temp <- common.shared(family)
      matlist[[i]] <- as(Matrix::forceSymmetric(temp), "dsCMatrix")
    }
    result <- Matrix::bdiag(matlist)

    if (any(duplicated(id$id))) {
      dimnames(result) <- list(NULL, paste(id$famid, id$id,
                                           sep = "/"))
    } else { dimnames(result) <- list(id$id, id$id) }
    result
}

#' @rdname common.shared
#' @export
common.shared.pedigree <- function (id, ...)
{
  n <- length(id$id)

  if (n == 1)
    return(matrix(1, 1, 1, dimnames = list(id$id, id$id)))
  if (any(duplicated(id$id)))
    stop("All id values must be unique")

  temp <- matrix(rep(1, n*n), n)
  dimnames(temp) <- list(id$id, id$id)
  temp
}
