tri_error_numeric <- function(f_nm, f, min, max, mode) {
  title <- paste0("Arguments ",
                  f_nm,
                  ", min, max, and mode must have numeric values.")
  msg <- tri_error_bullets(
    title,
    tri_error_bullet("Type", typeof, f, f_nm),
    tri_error_bullet("Type", typeof, min, "min"),
    tri_error_bullet("Type", typeof, max, "max"),
    tri_error_bullet("Type", typeof, mode, "mode")
  )
  tri_error(msg)
}

tri_error_logical <- function(log) {
  msg <- tri_error_bullets(
    "Argument log must have a single logical value.",
    tri_error_bullet("Size", length, log, "log"),
    tri_error_bullet("Type", typeof, log, "log")
  )
  tri_error(msg)
}

tri_error_logical2 <- function(lower_tail, log_p) {
  msg <- tri_error_bullets(
    "Arguments lower_tail and log_p must each have a single logical value.",
    tri_error_bullet("Size", length, lower_tail, "lower_tail"),
    tri_error_bullet("Size", length, log_p, "log_p"),
    tri_error_bullet("Type", typeof, lower_tail, "lower_tail"),
    tri_error_bullet("Type", typeof, log_p, "log_p")
  )
  tri_error(msg)
}

tri_error_recycle <-
  function(f_nm, f, min, max, mode, value = FALSE) {
    title <- paste0("Arguments ",
                    f_nm,
                    ", min, max, and mode must have compatible sizes.")
    msg <- tri_error_bullets(
      title,
      if (value) {
        tri_error_bullet("Value", c, f, f_nm)
      } else {
        tri_error_bullet("Size", length, f, f_nm)
      },
      tri_error_bullet("Size", length, min, "min"),
      tri_error_bullet("Size", length, max, "max"),
      tri_error_bullet("Size", length, mode, "mode"),
      i = "Only min, max, and mode values of size one are recycled."
    )
    tri_error(msg)
  }

tri_error_n <- function(n) {
  msg <- tri_error_bullets(
    "Argument n must have a non-zero positive numeric value.",
    tri_error_bullet("Size", length, n, "n"),
    tri_error_bullet("Positive", all, n > 0, "n")
  )
  tri_error(msg)
}

tri_error_bullet <- function(prefix, f, arg, arg_nm) {
  paste0(prefix, " ", f(arg), ": Argument ", arg_nm, ".")
}

#' @importFrom rlang format_error_bullets
#' @importFrom vctrs vec_c
tri_error_bullets <- function(title, ...) {
  error_bullets <-
    format_error_bullets(vec_c(..., .name_spec = "{outer}"))
  paste0(title, "\n", error_bullets)
}

tri_error_class <- function(class) {
  c(paste0("tri_error_", class), "tri_error")
}

#' @importFrom rlang as_name error_cnd
tri_error <- function(x) {
  call <- sys.call(-1)
  fn_name <- as_name(call[[1]])
  class <- tri_error_class(gsub("^error_", "", fn_name))
  error_cnd(class, message = x)
}
