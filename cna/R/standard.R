
# stdCond: standardize conditions with Rcpp functions
#   x   character Vector with 1-sided formulas. No blanks allowed!
stdCond <- function(x){
  l <- strsplit(noblanks(x), "+", fixed = TRUE)
  if (length(l) == 0L) return(character(0))
  l <- lapply(l, unique.default)
  u <- unlist(l, use.names = FALSE, recursive = FALSE)
  out <- 
    C_relist_Char(C_mconcat(strsplit(u, "*", fixed = TRUE), "*", sorted = TRUE),
                  lengths(l))
  C_mconcat(out, "+", sorted = TRUE)
}
# matchCond: find formulas in x with no equivalent formula in y
#  x,y  character Vector with 1-sided formulas. No blanks allowed!
matchCond <- function(x, y){
  stdCond(x) %in% unique(stdCond(y))
}
