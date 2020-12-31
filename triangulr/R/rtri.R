# Random Variate Generator

#' @include Triangular.R
#' @rdname Triangular
#' @export
rtri <- function(n,
                 min = 0,
                 max = 1,
                 mode = 0.5,
                 dqrng = FALSE) {
  if (!is.numeric(n) || !is.numeric(min) || !is.numeric(max) ||
      !is.numeric(mode)) {
    cnd_signal(tri_error_numeric("n", n, min, max, mode))
  }

  if (length(n) != 1L || n < 1L) {
    cnd_signal(tri_error_n(n))
  }

  if (length(min) == 1L &&
      length(max) == 1L && length(mode) == 1L) {
    RTriC(n, min, max, mode, dqrng)
  } else {
    n <- as.integer(n)
    tryCatch({
      params <- vec_recycle_common(min, max, mode, .size = n)
    }, error = function(c) {
      cnd_signal(tri_error_recycle("n", n, min, max, mode, value = TRUE))
    })
    RTriC2(n, params[[1L]], params[[2L]], params[[3L]], dqrng)
  }
}
