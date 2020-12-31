
# function allCombs: 
# ==================
# Generates a data frame of all possible combinations of length(x) multi-valued factors
# the factors are labelled LETTERS[seq_along(x)] and x is taken as the number of levels of the factors
# Arguments:
#   x   integer vector with values >0
# Value: A data frame
allCombs <- function(x){    
  stopifnot(x>0, x%%1 == 0)
  out <- do.call(expand.grid, c(lapply(x, seq_len), KEEP.OUT.ATTRS = FALSE))
  names(out) <- make.unique(rep(LETTERS, length.out = ncol(out)))
  out
}


# function makeFuzzy:
# ==================
# Subsets a data frame x based on a condition cond
# Arguments:
#   x          data frame or configTable
#   fuzzvalues values to add/substract at random from elements of x
# Value: resulting configTable of type "fs"
makeFuzzy <- function(x, fuzzvalues = c(0, .05, .1), ...){
  if (is.null(colnames(x))){
    colnames(x) <- make.unique(LETTERS[seq_len(ncol(x))])
  }
  xfuzzy <- ct2df(x)
  ulx <- unlist(xfuzzy, use.names = FALSE)
  xfuzzy[, ] <- ifelse(as.logical(ulx),
                       1 - sample(fuzzvalues, length(ulx), replace = TRUE),
                       sample(fuzzvalues, length(ulx), replace = TRUE))
  configTable(xfuzzy, type = "fs", rm.dup.factors = FALSE, rm.const.factors = FALSE, ...)
}

