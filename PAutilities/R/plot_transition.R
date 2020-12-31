#' Plot the transitions and matchings from a \code{transition} object
#'
#' @param x the object to plot
#' @param ... further methods passed to or from methods, currently unused
#'
#' @return A plot of the predicted and actual transitions in a \code{transition}
#'   object, as well as the matchings between them
#' @export
#'
#' @examples
#' predictions <- (sample(1:100)%%2)
#' references  <- (sample(1:100)%%2)
#' window_size <- 7
#' transitions <- get_transition_info(predictions, references, window_size)
#' plot(transitions)
plot.transition <- function(x, ...) {

  x$predictions <-
    as.character(x$predictions) %>%
    sapply(function(x) switch(
      x, "0" = 3, "1" = 2
    )) %>%
    unname(.)

  graphics::plot(
    seq(length(x$references)),
    x$references,
    pch = 16,
    ylab = "",
    yaxt = "n",
    ann = FALSE,
    xlab = "Reference",
    ylim = c(0.5,2.5)
  )

  graphics::mtext("Reference", 1, 3)

  graphics::axis(3)
  graphics::mtext("Prediction", line = 3)

  graphics::points(
    seq(length(x$predictions)),
    x$predictions
  )

  if (nrow(x$matchings) > 0) {
    sapply(
      seq(nrow(x$matchings)), function(i) {

        if (x$matchings$rejected[i]) {
          line_col <- "red"
        } else {
          line_col <- "blue"
        }

        graphics::lines(
          c(x$matchings$Reference_Index[i], x$matchings$Prediction_Index[i]),
          c(1,2),
          col = line_col
        )

      }
    )
  }
  invisible()

}
