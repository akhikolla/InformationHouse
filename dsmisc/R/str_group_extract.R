#' Extract Regular Expression Groups
#'
#' @param string string to extract from
#' @param pattern pattern with groups to match
#' @param group groups to extract
#' @param nas return NA values (TRUE) or filter them out (FALSE)
#'
#' @return string vector or string matrix
#'
#' @export
#'
#' @examples
#' 
#' strings <- paste(LETTERS, seq_along(LETTERS), sep = "_")
#' str_group_extract(strings, "([\\w])_(\\d+)")
#' str_group_extract(strings, "([\\w])_(\\d+)", 1)
#' str_group_extract(strings, "([\\w])_(\\d+)", 2)
#' 
str_group_extract <- 
  function(string, pattern, group = NULL, nas = TRUE){
    
    # handle group option
    if ( is.null(group) ){
      m <- stringr::str_match(string, pattern)[,1]
    } else {
      m <- stringr::str_match(string, pattern)[, (group + 1) ]
    }
    
    # handle 
    if ( nas == FALSE ) {
      m <- m[ !is.na(m) ]
    }
    
    # return
    m
  }

