#' Printing and timing utility for managing processes
#'
#' @param part character scalar, either \code{Start} or \code{End}.
#' @param ... character strings to print. Default messages will print if no
#'   arguments are provided.
#' @param timer a proc_time object. Required for \code{manage_procedure} only if
#'   using the default message for \code{part = "End"} default message.
#' @param verbose logical. Print to console?
#'
#' @return For \code{part = "Start"}, a proc_time object (i.e., a timer passable
#'   to an eventual \code{part = "End"} command); for \code{part = "End"},
#'   invisible
#' @export
#'
#' @examples
#'
#' manage_procedure("Start", "String will be printed\n")
#' timer <- manage_procedure(
#' "Start", "Printing a string is optional", verbose = FALSE
#' )
#'
#' ## Default starting message
#' manage_procedure("Start")
#'
#' ## Default ending message
#' manage_procedure("End", timer = timer)
#'
#' ## Other examples
#' get_duration(timer)
#' manage_procedure("End", "Custom ending message")
#'
manage_procedure <- function(
  part = c("Start", "End"),
  ..., timer = NULL, verbose = TRUE
) {

  part <- match.arg(part)
  default <- !length(list(...))

  ## Default starting message
  if (all(default, verbose, part == "Start")) {
    cat("\nBeginning new process...\n")
  }

  ## Default ending message
  if (all(default, verbose, part == "End")) {
    if(is.null(timer)) {

      warning(paste(
        "\n`timer` argument must be provided",
        "in order to use\n  default message with",
        "part = \"End\""
      ))

    } else {

      cat(
        "\n...Process successful. Elapsed time",
        get_duration(timer), "minutes.\n"
      )

    }
  }

  ## Manage the printing
  if (verbose) cat(...)

  ## Manage the returning
  if (part == "Start") return(proc.time()) else invisible()

}

#' @rdname manage_procedure
#'
#' @export
get_duration <- function (timer) {

  format(
    (proc.time() - timer)[3]/60,
    digits = 1,
    nsmall = 1
  )

}
