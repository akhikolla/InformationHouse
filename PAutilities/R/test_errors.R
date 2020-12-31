#' Compare numeric variables in a data frame based on root-squared differences
#'
#' @param reference a data frame giving reference data
#' @param target a data frame giving target data
#' @param vars character vector of variable names to compare in each data frame
#' @param tolerance allowable difference between numeric values
#' @param return_logical logical. Should result be given as a logical vector
#'   (indicating TRUE/FALSE equality within tolerance) or a data frame of error
#'   summary values?
#'
#' @return If \code{return_logical = TRUE}, a named logical vector with one
#'   element per variable compared, indicating whether the maximum and
#'   root-mean-squared differences fall within the tolerance. If
#'   \code{return_logical = FALSE}, a data frame indicating the variables
#'   compared and the maximum and root-mean-squared differences.
#' @note It is assumed that \code{reference} and \code{target} have equal
#'   numbers of rows.
#' @export
#'
#' @examples
#' reference <- data.frame(
#' a = 1:100, b = 75:174
#' )
#'
#' target <- data.frame(
#'   a = 0.001 + (1:100),
#'   b = 76:175
#' )
#'
#' test_errors(reference, target, c("a", "b"))
#' test_errors(reference, target, c("a", "b"), return_logical = FALSE)
test_errors <- function(
  reference, target, vars,
  tolerance = 0.001005, return_logical = TRUE
) {

  stopifnot(
    inherits(reference, "data.frame"),
    inherits(target, "data.frame")
  )
  stopifnot(
    nrow(reference) == nrow(target),
    all(vars %in% names(reference)),
    all(vars %in% names(target))
  )
  stopifnot(
    all(sapply(reference[ ,vars], is.numeric)),
    all(sapply(target[ ,vars], is.numeric))
  )


  results <- sapply(
    vars,
    function(variable) {

      e <- reference[ ,variable] -
        target[ ,variable]
      se <- e^2

      if (anyNA(e)) {
        warning(paste(
          sum(is.na(e)), "missing error values",
          "will be ignored"
        ))
      }

      rmse <- sqrt(mean(se, na.rm = TRUE))
      max_e <- max(e, na.rm = TRUE)

      if (return_logical) {

        results <- all(
          rmse <= tolerance,
          max_e <= tolerance
        )

      } else {

        results <- data.frame(
          variable = variable,
          max_error = max_e,
          rmse = rmse,
          stringsAsFactors = FALSE,
          row.names = NULL
        )

      }

      results

    },
    simplify = FALSE
  )

  if (return_logical) results <- do.call(c, results)
  if (!return_logical) results <- do.call(rbind, results)

  results

}
