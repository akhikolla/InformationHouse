#' Perform equivalence testing on paired samples
#'
#' @param x numeric vector representing the (possibly surrogate) sample
#' @param y numeric vector representing the (possibly criterion) sample. Index
#'   paired with \code{x}
#' @param y_type classification of \code{y} for the purpose of analysis. Can be
#'   \code{"criterion"}, \code{"comparison"}, or \code{"both"}.
#' @param alpha the alpha level for the test
#' @param na.rm logical. Omit mean values for mean calculations?
#' @param scale character specifying whether the test should occur on an
#'   absolute or relative scale
#' @param absolute_region_width the region width for use when \code{scale =
#'   "absolute"}
#' @param relative_region_width the region width for use when \code{scale =
#'   "relative"}
#' @param ... further arguments passed to methods. Currently unused.
#'
#' @return a `paired_equivalence` object summarizing the test input and results
#' @note If a value is not specified for the region width that corresponds with
#'   \code{scale}, a default value will be assigned with a warning.
#' @references \href{https://www.ncbi.nlm.nih.gov/pubmed/29135817}{Dixon et al.}
#' @export
#'
#' @examples
#' set.seed(1544)
#' x <- data.frame(
#'   var1 = rnorm(500, 15, 4),
#'   var2 = rnorm(500, 23, 7.3)
#' )
#' y <- rnorm(500, 17.4, 9)
#'
#' test_result <- paired_equivalence_test(
#'   x, y, relative_region_width = 0.25
#' )
#'
#' lapply(test_result, head)
#'
paired_equivalence_test <- function(
  x, y, y_type = c("both", "criterion", "comparison"),
  alpha = 0.05, na.rm = TRUE,
  scale = c("relative", "absolute"),
  absolute_region_width = NULL,
  relative_region_width = NULL,
  ...
) {

    criterion_results <- NULL
    comparison_results <- NULL

    UseMethod("paired_equivalence_test", x)

}
