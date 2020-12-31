#' A table function to use with magrittr pipes
#'
#' Accepts a data frame as input and computes a contingency table for direct use in combination with the magrittr package.
#'
#' categorize is a wrapper to xtabs or table such that a data frame can be given as the first argument.
#'
#' @param .data A data frame
#' @param ... A formula (as in xtabs) or one or more objects which can be interpreted as factors (including character strings), or a list (or data frame) whose components can be so interpreted.
#' @return A table (possibly as an xtabs class if a model formula was used)
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' if (requireNamespace("magrittr", quietly = TRUE)) {
#'     library(magrittr)
#'
#'     esoph %>% categorize(alcgp, agegp)
#'     esoph %>% categorize(~ alcgp + agegp)
#' }
#'
#' @export
categorize <- function(.data, ...) {
    dots <- eval(substitute(alist(...)));
    if ("call" %in% sapply(dots, class))
        xtabs(..., data=.data)
    else {
        with(.data, do.call("table", dots))
    }
}
