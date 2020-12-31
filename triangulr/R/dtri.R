# Probability Density Function

#' @include Triangular.R
#' @rdname Triangular
#' @export
dtri <- function(x,
                 min = 0,
                 max = 1,
                 mode = 0.5,
                 log = FALSE) {
  if (!is.numeric(x) || !is.numeric(min) || !is.numeric(max) ||
      !is.numeric(mode)) {
    cnd_signal(tri_error_numeric("x", x, min, max, mode))
  }

  if (!is.logical(log) || length(log) > 1L) {
    cnd_signal(tri_error_logical(log))
  }

  if (length(min) == 1L &&
      length(max) == 1L && length(mode) == 1L) {
    DTriC(x, min, max, mode, log)
  } else {
    tryCatch({
      params <- vec_recycle_common(min, max, mode, .size = length(x))
    }, error = function(c) {
      cnd_signal(tri_error_recycle("x", x, min, max, mode))
    })
    DTriC2(x, params[[1L]], params[[2L]], params[[3L]], log)
  }
}
