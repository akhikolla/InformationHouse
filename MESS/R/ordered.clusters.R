#' Check if unique elements of a vector appear in contiguous clusters
#' 
#' \code{ordered.clusters} determines if identical elements of a vector appear
#' in contiguous clusters, and returns \code{TRUE} if the do and \code{FALSE}
#' otherwise.
#' 
#' 
#' @param id a vector
#' @return The function returns TRUE if the elements appear in contiguous
#' clusters and FALSE otherwise
#' @author Claus Ekstrom \email{claus@@ekstroem.dk} with suggestions from Peter
#' Dalgaard.
#' @seealso \code{\link{duplicated}}
#' @keywords utilities
#' @examples
#' 
#' x <- c(1, 1, 1, 2, 2, 3, 4, 1, 5, 5, 5)
#' ordered.clusters(x)
#' ordered.clusters(sort(x))
#' ordered.clusters(x[order(x)])
#' 
#' @export ordered.clusters
ordered.clusters <- function(id) {
  d <- which(duplicated(id)) # NB: d > 1, poss. empty
  all(id[d]==id[d-1])
}
