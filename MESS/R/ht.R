#' Show the head and tail of an object
#'
#' Show both the head and tail of an R object
#'
#' @param x The object to show
#' @param n The number of elements to list for the head
#' @param m The number of elements to list for the tail
#' @param returnList Logical. Should the result be returned as a list
#' @param ...  additional arguments passed to functions (not used at the moment)
#' @return NULL unless returnList is set to TRUE in which case a list is returned
#' @details This function does no error checking and it is up to the user to ensure that the input is indeed symmetric, positive-definite, and a matrix.
#'
#' @examples
#' ht(trees)
#' ht(diag(20))
#' ht(1:20)
#' ht(1:20, returnList=TRUE)
#'
#' @author Claus Ekstrom, \email{claus@@rprimer.dk}.
#' @keywords print
#' @export

ht <- function(x, n = 6L, m=n, returnList=FALSE, ...) {

    if (!returnList) {
        print(head(x, n))
        cat("...\n")
        print(tail(x, m))
        invisible(NULL)
    }
    else {
        list(head = head(x,n), tail = tail(x,m))
    }
}



