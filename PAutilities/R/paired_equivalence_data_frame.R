#' @rdname paired_equivalence_test
#' @export
paired_equivalence_test.data.frame <- function(
  x, y, y_type = c("both", "criterion", "comparison"),
  alpha = 0.05, na.rm = TRUE,
  scale = c("relative", "absolute"),
  absolute_region_width = NULL,
  relative_region_width = NULL,
  ...
) {

  object <- new_paired_equivalence(
    x, y, y_type, alpha, na.rm, scale,
    absolute_region_width, relative_region_width
  )

  column_wise <- sapply(
    object$x,
    paired_equivalence_test,
    y = object$y, y_type = object$y_type, alpha = object$alpha,
    na.rm = object$na.rm, scale = object$scale,
    absolute_region_width = object$absolute_region_width,
    relative_region_width = object$relative_region_width,
    simplify = FALSE
  )

  column_wise <- lapply(
    seq(column_wise), function(n) {
      result <- column_wise[[n]]$result
      result$variable <- names(column_wise)[n]
      result[ ,c("variable", setdiff(names(result), "variable"))]
    }
  )

  results <- do.call(rbind, column_wise)

  object$results <- results

  object

}
