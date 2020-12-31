#' Expand table or matrix to data frame
#'
#' Expands a contingency table to a data frame where each observation in the table becomes a single observation in the data frame with corresponding information for each for each combination of the table dimensions.
#'
#' @param x A table or matrix
#' @return A data frame with the table or matrix expanded
#' @author Claus Ekstrom \email{claus@@rprimer.dk}
#' @keywords manip
#' @examples
#'
#' expand_table(diag(3))
#' m <- matrix(c(2, 1, 3, 0, 0, 2), 3)
#' expand_table(m)
#' result <- expand_table(UCBAdmissions)
#' head(result)
#'
#' # Combine into table again
#' xtabs(~Admit + Gender + Dept, data=result)
#'
#' @export
expand_table <- function(x) {

    if (!any(c("matrix", "table") %in% class(x)))
        stop("needs matrix or table as input")

    ndim <- length(dim(x))

    if ("matrix" %in% class(x))
        x <- as.data.frame(as.table(x))

    if ("table" %in% class(x))
        x <- as.data.frame(x)

#    x[rep(seq.int(1,nrow(x)), x$Freq), 1:2]
    x[rep(seq.int(1,nrow(x)), x$Freq), seq.int(ndim)]
}
