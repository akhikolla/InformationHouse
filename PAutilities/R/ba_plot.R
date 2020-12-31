#' Create a Bland-Altman plot
#'
#' @param plotdata dataframe from which to build the plot
#' @param x_var character expression to evaluate for the x-axis
#' @param y_var character expression to evaluate for the y-axis
#' @param x_name axis label for the x-axis
#' @param y_name axis label for the y-axis
#' @param shape numeric. The point shape to display.
#' @param ... further arguments passed to \code{theme}
#'
#' @return a Bland-Altman plot
#' @export
#'
#' @references Bland, J. M., & Altman, D. G. (1986). Statistical methods for
#'   assessing agreement between two methods of clinical measurement. lancet,
#'   1(8476), 307-310.
#'
#' @examples
#' data(ex_data, package = "PAutilities")
#'
#' # Reduce the number of data points (for illustration purposes) by isolating
#' # the 150 largest cases
#'
#' illustration_threshold <-
#'     quantile(ex_data$Axis1, probs = 1 - (150 / nrow(ex_data)))
#' ex_data <- ex_data[ex_data$Axis1 > illustration_threshold, ]
#'
#' # Generate the plot
#' my_ba <- ba_plot(
#'     ex_data,
#'     "(Axis1 + Axis3) / 2",
#'     "Axis1 - Axis3",
#'     "mean(Axis1, Axis3)",
#'     "Axis1 - Axis3"
#' )
#'
#' my_ba
#'
#' \donttest{
#' # You can add to the plot as you would a normal ggplot object
#'     my_ba +
#'       ggplot2::geom_text(
#'       x = 2000, y = 9000, label = "A",
#'       size = 8, fontface = "bold", colour = "blue"
#'       )
#'
#' # With caution, you can change some automatic options (e.g. color of
#' # regression line) by overwriting in a new layer
#'
#'     my_ba + ggplot2::geom_smooth(method = "lm", se = FALSE, colour = "blue")
#'
#' }
ba_plot <- function(
  plotdata, x_var, y_var, x_name, y_name, shape = 16, ...
) {

  ggplot(plotdata, aes_string(x = x_var, y = y_var)) +
  geom_point(shape = shape) + theme_classic() +
  theme(axis.line = element_line(size = .5)) +
  scale_y_continuous(
    name = y_name
  ) +
  scale_x_continuous(
    name = x_name
  ) +
  geom_hline(
    yintercept =
      lazyeval::f_eval(
        ~mean(eval(parse(text = y_var)), na.rm = TRUE),
        data = plotdata
      ),
    size = 1.2
  ) +
  geom_smooth(
    method = 'lm', se = FALSE,
    colour = 'black'
  ) +
  geom_hline(
    aes(
      yintercept =
        lazyeval::f_eval(
          ~mean(eval(parse(text = y_var)), na.rm = TRUE) +
            (1.96*stats::sd(eval(parse(text = y_var)), na.rm = TRUE)),
          data = plotdata)
    ), size = 1.3, linetype = 'dashed'
  ) +
  geom_hline(
    aes(
      yintercept =
        lazyeval::f_eval(
          ~mean(eval(parse(text = y_var)), na.rm = TRUE) -
            (1.96*stats::sd(eval(parse(text = y_var)), na.rm = TRUE)),
          data = plotdata
        )
    ), size = 1.3, linetype = 'dashed'
  ) +
  theme(
    axis.title = element_text(size = 14, face = 'bold'),
    axis.text = element_text(size = 12)
  ) +
  theme(...)
}
