# Moment Generating Function

#' @include Triangular.R
#' @rdname Triangular
#' @export
mgtri <- function(t,
                  min = 0,
                  max = 1,
                  mode = 0.5) {
  if (!is.numeric(t) || !is.numeric(min) || !is.numeric(max) ||
      !is.numeric(mode)) {
    cnd_signal(tri_error_numeric("t", t, min, max, mode))
  }

  if (length(min) == 1L &&
      length(max) == 1L && length(mode) == 1L) {
    MGTriC(t, min, max, mode)
  } else {
    tryCatch({
      params <- vec_recycle_common(min, max, mode, .size = length(t))
    }, error = function(c) {
      cnd_signal(tri_error_recycle("t", t, min, max, mode))
    })
    MGTriC2(t, params[[1L]], params[[2L]], params[[3L]])
  }
}
