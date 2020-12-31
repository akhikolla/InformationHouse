# Quantile Function

#' @include Triangular.R
#' @rdname Triangular
#' @export
qtri <- function(p,
                 min = 0,
                 max = 1,
                 mode = 0.5,
                 lower_tail = TRUE,
                 log_p = FALSE) {
  if (!is.numeric(p) || !is.numeric(min) || !is.numeric(max) ||
      !is.numeric(mode)) {
    cnd_signal(tri_error_numeric("p", p, min, max, mode))
  }

  if (!is.logical(lower_tail) ||
      length(lower_tail) > 1L || !is.logical(log_p)
      || length(log_p) > 1L) {
    cnd_signal(tri_error_logical2(lower_tail, log_p))
  }

  if (length(min) == 1L &&
      length(max) == 1L && length(mode) == 1L) {
    QTriC(p, min, max, mode, lower_tail, log_p)
  } else {
    tryCatch({
      params <- vec_recycle_common(min, max, mode, .size = length(p))
    }, error = function(c) {
      cnd_signal(tri_error_recycle("p", p, min, max, mode))
    })
    QTriC2(p, params[[1L]], params[[2L]], params[[3L]], lower_tail, log_p)
  }
}
