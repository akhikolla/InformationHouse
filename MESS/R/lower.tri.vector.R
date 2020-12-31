#' Split Matrix by Clusters and Return Lower Triangular Parts as Vector
#' 
#' Split a matrix into block diagonal sub matrices according to clusters and
#' combine the lower triangular parts into a vector
#' 
#' 
#' @param x a square matrix
#' @param cluster numeric or factor. Is used to identify the sub-matrices of
#' \code{x} from which the lower triangular parts are extracted. Defaults to
#' the full matrix.
#' @param diag logical. Should the diagonal be included?
#' @return Returns a numeric vector containing the elements of the lower
#' triangular sub matrices.
#' @author Claus Ekstrom \email{claus@@ekstroem.dk}
#' @seealso \code{\link{lower.tri}}
#' @keywords manip
#' @examples
#' 
#' m <- matrix(1:64, ncol=8)
#' cluster <- c(1, 1, 1, 1, 2, 2, 3, 3)
#' lower.tri.vector(m, cluster)
#' 
#' @export lower.tri.vector
lower.tri.vector <- function(x, cluster=rep(1,nrow(x)), diag = FALSE) {

  if (!is.matrix(x) & nrow(x) != ncol(x))
    stop("Input must be a square matrix")

  # Simple check
  if (length(cluster) != nrow(x))
    stop("Length of cluster vector must match dimensions of matrix")

  # Remember to check new function on the ordering
  if (!ordered.clusters(cluster))
    stop("the cluster elements should be in contiguous blocks")
    
  unlist(lapply(unique(cluster), 
         function(id) { sel <- (id==cluster) ;
                        m <- x[sel,sel] ;
                        # as.vector(m[lower.tri(m, diag=diag)])
                        # Use this solution which works with sparse matrices
                        as.vector(m[which(lower.tri(m, diag=diag)==TRUE,arr.ind=TRUE)])
                      }
         )
  )
}
