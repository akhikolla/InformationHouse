#' Simulate randomized urn design
#'
#' Simulates a randomized treatment based on an urn model.
#'
#' The urn model can be described as follows: For k different treatments, the urn design is initiated with a number of balls in an urn corresponding to the start weight (the alpha argument), where each treatment has a specific colour. Whenever a patient arrives, a random ball is drawn from the urn and the colour decides the treatment for the patient. For each of the treatments that weren't chosen we add beta balls of the corresponding colour(s) to the urn to update the probabilities for the next patient.
#'
#' @param n the number of individuals to randomize
#' @param alpha a non-negative integer vector of weights for each treatment group. The length of the vector corresponds to the number of treatment groups.
#' @param beta a non-negative integer of weights added to the groups that were not given treatment
#' @param labels a vector of treatment labels. Must be the same length as the length of alpha.
#' @param data.frame A logical that determines if the function should return a vector of group indices (the default, if FALSE) or a data frame (if TRUE).
#' @param startid margin paramaters; vector of length 4 (see \code{\link[graphics]{par}})
#'
#' @return A vector with group indices. If the argument \code{data.frame=TRUE} is used then a data frame with three variables is returned: id, group, and treatment (the group label). 
#'
#' @examples
#' rud(5)
#' rud(5, alpha=c(1,1,10), beta=5)
#'
#' @export
rud <- function(n, alpha=c(1, 1), beta=1, labels=seq(1, length(alpha)), data.frame=FALSE, startid=1) {

    ncat <- length(alpha)  # Number of categories
    
    if (ncat<2)
        stop("alpha mst have a length of at least 2")
    
    if (length(labels) != length(alpha))
        stop("length of labels must be the same as length of groups")
    
    urn <- alpha
    group <- seq(1, ncat)
    res <- integer(n)
    
    for (i in 1:n) {
        res[i] <- sample(group, 1, prob=urn/(sum(urn)))
        urn[-res[i]] <- urn[-res[i]] + beta
    }
    
    if (data.frame)
        data.frame(id=(1:n)+(startid-1), group=res, treatment=labels[res])
    else
        res
}
