#' Plot \code{evgram} object
#' 
#' Plots \code{evgram} object produced by  
#' \code{\link[gear]{evgram}} using the 
#' \code{\link[graphics]{plot}} function.
#' 
#' The arguments after \code{...} are only considered
#' if a directional semivariogram is provided, i.e., 
#' when \code{x$ndir != 1}. 
#' 
#' \code{split} will split directional semivariograms into
#' different plots automatically using 
#' \code{\link[autoimage]{autosize}}.
#' 
#' \code{add_legend} and \code{args_legend} are only used
#' for a directional semivariogram when \code{split = FALSE}
#' 
#' @param x An \code{evgram} object produced by the \code{\link[gear]{evgram}} function.
#' @param ... Additional arguments to pass the \code{\link[graphics]{plot}} function to change aspects of the plot.
#' @param split A logical value indicating whether, for a directional 
#' semivariogram, the directional semivariograms should be displayed 
#' in a single or split panels.  Default is FALSE, for a single panel.
#' @param add_legend A logical value indicating whether a 
#' legend should be included for a plot of directional
#' semivariograms when \code{split = FALSE}.
#' @param args_legend A list with arguments for 
#' \code{\link[graphics]{legend}}.
#' @return NULL
#' @author Joshua French
#' @export
#' @method plot evgram
#' @seealso \code{\link[lattice]{xyplot}}, \code{\link[gear]{evgram}}
#' @examples 
#' data(co)
#' # omnidirectional example
#' v = evgram(Al ~ 1, co, ~ easting + northing)
#' plot(v)
#' plot(v, main = "semivariogram of Al")
#' 
#' # directional semivariograms overlaid
#' v2 = evgram(Al ~ 1, co, ~ easting + northing,
#'             angle = 22.5, ndir = 4)
#' plot(v2)
#' plot(v2, ylab = "semi-variance", pch = 2:5, type = "p")
#' plot(v2, lty = 2:5, type = "l", 
#'      args_legend = list(x = "bottomright",
#'                         legend = c("22.5", "67.5", "112.5", "157.5")))
#' # directional semivariograms split
#' plot(v2, split = TRUE)
#' plot(v2, split = TRUE, col = 2, pch = 3, type = "b",
#'      main = c("(a)", "(b)", "(c)", "(d)"))
plot.evgram = function(x, ..., split = FALSE,
                       add_legend = TRUE,
                       args_legend = list()) {
  arg_check_split(split)
  arg_check_add_legend(add_legend)
  if (!is.list(args_legend)) {
    stop("args_legend must be a list")
  }
  
  # setup arguments for plotting
  largs = list(...)
  if (is.null(largs$ylim)) {
    largs$ylim = c(0, max(x$semivariogram$semivariance))
  }
  if (is.null(largs$xlim)) {
    largs$xlim = c(0, max(x$semivariogram$distance))
  }
  largs$y = x$semivariogram$semivariance
  largs$x = x$semivariogram$distance
  if (is.null(largs$xlab)) {
    largs$xlab = "distance"
  }
  if (is.null(largs$ylab)) {
    largs$ylab = "semivariance"
  }  
  
  # different plotting approaches depending on options
  if (x$ndir == 1)   {
    do.call(graphics::plot, largs)
  } else if (split) {
    old_par = graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par))
    graphics::par(mfrow = autoimage::autosize(x$ndir))
    split_semivariogram = split(x$semivariogram, x$semivariogram$angle)
    temp_largs = largs
    if (is.null(largs$main)) largs$main = names(split_semivariogram)
    if (length(largs$main) != x$ndir) {
      stop("length of main != x$ndir")
    }
    
    for (i in seq_along(split_semivariogram)) {
      temp_largs$y = split_semivariogram[[i]]$semivariance
      temp_largs$x = split_semivariogram[[i]]$distance
      temp_largs$main = largs$main[i]
      do.call(graphics::plot, temp_largs)
    }
  } else {
    # blank initial plot with correct labeling
    temp_largs = largs
    temp_largs$type = "n"

    # setup various options for overlaying semivariograms
    if (is.null(largs$pch)) largs$pch = seq_len(x$ndir)
    if (length(largs$pch) != x$ndir) {
      stop("length of pch != x$ndir")
    }
    if (is.null(largs$col)) largs$col = seq_len(x$ndir)
    if (length(largs$col) != x$ndir) {
      stop("length of col != x$ndir")
    }
    if (is.null(largs$lty)) largs$lty = seq_len(x$ndir)
    if (length(largs$lty) != x$ndir) {
      stop("length of lty != x$ndir")
    }
    my_type = largs$type
    if (is.null(my_type)) {
      my_type = "b"
    }
    my_pch = largs$pch
    largs$pch = NULL
    my_col = largs$col
    largs$col = NULL
    my_lty = largs$lty
    largs$lty = NULL

    # fix plotting issues depending on type
    if (my_type == "l") {
      my_pch = NULL
    }
    if (my_type == "p") {
      my_lty = NULL
    }
    
    split_semivariogram = split(x$semivariogram, x$semivariogram$angle)
    do.call(graphics::plot, temp_largs) # main plot
    # overlay semivariograms
    for (i in seq_along(split_semivariogram)) {
      graphics::points(semivariance ~ distance,
                       data = split_semivariogram[[i]],
                       col = my_col[i], pch = my_pch[i],
                       lty = my_lty[i], type = my_type)
    }
    if (add_legend) {
      args_legend$pch = my_pch
      args_legend$lty = my_lty
      args_legend$col = my_col
      if (is.null(args_legend$x)) {
        args_legend$x = "topleft"
      }
      if (is.null(args_legend$legend)) {
        args_legend$legend = names(split_semivariogram)
      }
      do.call(graphics::legend, args_legend)
    }
  }
}
