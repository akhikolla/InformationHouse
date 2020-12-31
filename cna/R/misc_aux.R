
# ==== General auxiliary functions ====

# m_any, m_all, m_is.null: 
# shortcuts to vapply(..., any, logical(1)), etc.
m_any <- function(x, na.rm = FALSE, USE.NAMES = FALSE){
  vapply(x, any, logical(1), na.rm = na.rm, USE.NAMES = USE.NAMES)
}
m_all <- function(x, na.rm = FALSE, USE.NAMES = FALSE){
  vapply(x, all, logical(1), na.rm = na.rm, USE.NAMES = USE.NAMES)
}
m_is.null <- function(x, USE.NAMES = FALSE){
  vapply(x, is.null, logical(1), USE.NAMES = USE.NAMES)
}

# function isConstant
# returns TRUE if all values in x are equal, FALSE otherwise
isConstant <- function(x) if (length(x)) all(x==x[1]) else TRUE


# function isASCII
# Returns TRUE, if char vector x has only ASCII characters, otherwise FALSE
isASCII <- function(x) !anyNA(iconv(x, "", "ASCII"))