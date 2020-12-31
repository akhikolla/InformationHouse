#' Perform a spurious curve analysis
#'
#' Assess performance using the Transition Pairing Method when the spurious
#' pairing threshold is varied
#'
#' @param trans a \code{transition} object
#' @param predictions vector of predictions indicating transition \code{(1)} or
#'   non-transition \code{(2)}
#' @param references vector of criteria indicating transition \code{(1)} or
#'   non-transition \code{(2)}
#' @param thresholds the threshold settings to test
#'
#' @return an object with class \code{spurious_curve}
#' @export
#'
#' @examples
#' set.seed(100)
#' predictions <- (sample(1:100)%%2)
#' references  <- (sample(1:100)%%2)
#' \donttest{
#' trans <- get_transition_info(
#'   predictions, references, 7
#' )
#' head(spurious_curve(trans))
#' }
spurious_curve <- function(
  trans, predictions, references, thresholds = 1:20
) {

  ## Validate input

    if (missing(trans)) {

      stopifnot(!any(missing(predictions), missing(references)))

    } else {

      if (any(!missing(predictions), !missing(references))) {
        warning(paste(
          "Ignoring inputted values for `predictions`",
          "and `references` -- drawing from elements of `info`"
        ))
      }
      predictions <- trans$predictions
      references  <- trans$references

    }

  ## Run the tests
  as.list(thresholds) %>%
  lapply(function(x) {
    x <- try(summary(
      get_transition_info(predictions, references, x)
    ))
    if ("try-error" %in% class(x)) return(NULL)
    x
  }) %>%
  structure(., class = append(class(.), "spurious_curve"))

}

#' Plot a spurious curve
#'
#' @param x a \code{spurious_curve} object
#' @param ... further arguments (currently unused)
#'
#' @return a plot of the object
#' @export
#'
#' @seealso \code{\link{spurious_curve}}
#'
#' @examples
#' set.seed(100)
#' predictions <- (sample(1:100)%%2)
#' references  <- (sample(1:100)%%2)
#' \donttest{
#' trans <- get_transition_info(
#'   predictions, references, 7
#' )
#' result <- spurious_curve(trans)
#' plot(result)
#' }
plot.spurious_curve <- function(x, ...) {

  x <- lapply(x, as, Class = "data.frame")

  plot_vars <- c(
    "window_size", "recall", "precision",
    "rmse_prop", "aggregated_performance"
  )

  values <-
    lapply(x, function(y) y[ ,plot_vars]) %>%
    do.call(rbind, .) %>%
    as.list(.)

  par(mfrow = c(2,2))

  plot(
    values$window_size, values$recall, type = "b",
    ylab = "", xlab = ""
  )
  title(
    ylab = "Recall",
    xlab = "Spurious Pairing Threshold",
    line = 2, font.lab = 2
  )

  plot(
    values$window_size, values$precision, type = "b",
    ylab = "", xlab = ""
  )
  title(
    ylab = "Precision",
    xlab = "Spurious Pairing Threshold",
    line = 2, font.lab = 2
  )

  plot(
    values$window_size, values$rmse_prop, type = "b",
    ylab = "", xlab = ""
  )
  title(
    ylab = expression(bold("RMSE"["%"])),
    xlab = "Spurious Pairing Threshold",
    line = 2, font.lab = 2
  )

  plot(
    values$window_size, values$aggregated_performance,
    type = "b", ylab = "", xlab = ""
  )
  title(
    ylab = "Aggregated Performance",
    xlab = "Spurious Pairing Threshold",
    line = 2, font.lab = 2
  )

}
