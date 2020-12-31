# Cumulative Distribution Function

#' @include Triangular.R
#' @rdname Triangular
#' @export
ptri <- function(q,
                 min = 0,
                 max = 1,
                 mode = 0.5,
                 lower_tail = TRUE,
                 log_p = FALSE) {
  if (!is.numeric(q) || !is.numeric(min) || !is.numeric(max) ||
      !is.numeric(mode)) {
    cnd_signal(tri_error_numeric("q", q, min, max, mode))
  }

  if (!is.logical(lower_tail) ||
      length(lower_tail) > 1L || !is.logical(log_p)
      || length(log_p) > 1L) {
    cnd_signal(tri_error_logical2(lower_tail, log_p))
  }

  if (length(min) == 1L &&
      length(max) == 1L && length(mode) == 1L) {
    PTriC(q, min, max, mode, lower_tail, log_p)
  } else {
    tryCatch({
      params <- vec_recycle_common(min, max, mode, .size = length(q))
    }, error = function(c) {
      cnd_signal(tri_error_recycle("q", q, min, max, mode))
    })
    PTriC2(q, params[[1L]], params[[2L]], params[[3L]], lower_tail, log_p)
  }
}
