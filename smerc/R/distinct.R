#' Distinct elements of a list
#'
#' \code{distinct} takes a list of integer vectors and
#' returns the list indices that contain unique combinations
#' of elements.  This function is NOT robust against misuse,
#' so please use properly.
#'
#' Assume that \code{k} is the largest integer value in
#' \code{x}.  A vector of the largest \code{k} prime numbers
#' is obtained (call this \code{pri}).  The algorithm takes
#' the sum of the log of \code{pri[x[[i]]]} for each element
#' of \code{x}, and determines which sums are unique.  This
#' is why the elements of \code{x} must be integer vectors.
#' The prime aspect of the algorithm is critical, as it
#' ensures that a none of the values are multiples of the
#' others, ensuring uniqueness.
#'
#' Note: this algorithm has only been applied to data sets
#' where each element of \code{x[[i]]} appears only once,
#' though it should work for repeats also.
#'
#' @param x A list of integers
#' @param N The largest integer value across all elements of
#'   \code{x}.
#'
#' @return A vector with the distinct indices.
#' @export
#' @author Joshua French
#' @references Algorithm based on suggestion at
#'   \url{https://stackoverflow.com/a/29824978}.
#' @examples
#' x = list(1:3, 3:1, 1:4, 4:1, c(1, 2, 4, 6), c(6, 4, 1, 2))
#' x[distinct(x)]
distinct = function(x, N = max(unlist(x))) {
  pri = randtoolbox::get.primes(N)
  if (capabilities()[["long.double"]]) {
    sums = sapply(x, function(xi) sum(log(pri[xi])))
  } else {
    sums = sapply(x, function(xi) sum(sort(log(pri[xi]))))
  }
  which(!duplicated(sums))
}
