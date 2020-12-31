#' Write a data frame in XML format
#'
#' Writes the data frame to a file in the XML format.
#'
#' @param data the data frame object to save
#' @param file the file name to be written to.
#' @param collapse logical. Should the output file be collapsed to make it fill less? (Defaults to TRUE)
#' @return None
#' @details This function does not require the \pkg{XML} package to be installed to function properly.
#'
#' @examples
#'
#' \dontrun{
#' data(trees)
#' write.xml(trees, file="mydata.xml")
#' }
#'
#' @author Claus Ekstrom, \email{claus@@rprimer.dk} based on previous work by Duncan Temple Lang.
#' @keywords file
#' @export write.xml

write.xml <- function(data, file=NULL, collapse=TRUE) {

  if(is.null(file))
    stop("filename not specified")

  if (!is.data.frame(data))
    stop("data must be a data frame")

  con <- file(file, "w")

  vnames <- names(data)
  pre    <- paste0(ifelse(collapse, "", "    "), "<", vnames, ">")
  post   <- paste0("</", vnames, ">")
  colchar <- ifelse(collapse, "", "\n")
  rowend  <- ifelse(collapse, "</row>\n", "  </row>")

  ## Write header
  writeLines('<?xml version="1.0"?>\n<document>', con=con)
  sapply(1:nrow(data), function(rowi)
      ## Start
      writeLines(c("  <row>", paste0(pre, sapply(data[rowi,], as.character), post, collapse=colchar), rowend), con=con, sep=colchar)
         )
  writeLines("</document>", con=con)
  close(con)

  invisible(0)
}
