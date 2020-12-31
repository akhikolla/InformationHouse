any_equal <- function(target, current, tolerance = sqrt(.Machine$double.eps)) {
  any(sapply(target, function(targ) abs(current - targ) < tolerance ))
}
all_equal <- function(target, current, tolerance = sqrt(.Machine$double.eps)) {
  all(sapply(target, function(targ) abs(current - targ) < tolerance ))
}


#' Convert a codeset to a contingency table
#'
#' @param x codeset 
#'
#' @return contingency table as a 2x2 matrix
#' 
#' @export
as.contingency.table <- function(x) {
  y <- NULL
  
  if ( nrow(x) > 2 && ncol(x) == 2 ) {
    tp = length(which(x[,1] == 1 & x[,2] == 1))
    tn = length(which(x[,1] == 0 & x[,2] == 0))
    fp = length(which(x[,1] == 0 & x[,2] == 1))
    fn = length(which(x[,1] == 1 & x[,2] == 0))
    
    y <- matrix(
      c(tp, fp, fn, tn), 
      ncol = 2, nrow = 2
    )
  }
  else {
    y <- x
  }
  
  dimnames(y) = list(c("first_positive", "first_negative"), c("second_positive", "second_negative"))
  class(y) <- c("contingency.table", "rating.set", class(y))
  attr(y, "agreement") <- list(
    "true_positive" = y[1, 1],
    "false_negative" = y[1, 2],
    "false_positive" = y[2, 1],
    "true_negative" = y[2, 2]
  )
  attr(y, "baserate") <- baserateCT(y)
  attr(y, "kappa") <- kappaCT(y)
  
  return(y)
}

#' Convert codeset to contingency table
#'
#' @param x matrix contingency table (2x2)
#'
#' @return 2-column matrix representation of the contingency table
#' 
#' @export
as.code.set <- function(x) {
  y <- NULL
  if (all(dim(x) == c(2,2))) {
    y <- matrix(
      c(rep(1, x[1,1]*2), rep(1:0, x[1,2]), rep(0:1, x[2,1]), rep(0, x[2,2]*2)),
      byrow = TRUE,
      ncol = 2
    )
  } 
  else {
    y <- x
  }
  
  class(y) <- c("code.set", "rating.set", class(y))
  attr(y, "baserate") <- baserateSet(y)
  attr(y, "kappa") <- kappaSet(y)
  colnames(y) <- c("first_rater", "second_rater")
  
  return(y)
}

#' Helper function to return special values on a rating set
#' 
#' @param x Set or Contingency.Table
#'
#' @param i Value to search for
#'
#' @export
"$.rating.set" <- function (x, i) {
  if (i %in% names(attributes(x))) {
    return(attr(x, i))
  }
  else {
    if (is(x, "code.set")) {
      if (i %in% colnames(x)) {
        return(x[, which(colnames(x) == i)])
      }
    }
    else if (is(x, "contingency.table")) {
      if (i %in% colnames(x)) {
        return(sum(x[, which(colnames(x) == i)]))
      }
      else if (i %in% rownames(x)) {
        return(sum(x[which(rownames(x) == i), ]))
      }
    }
  }
  
  return(NULL)
}

#' @export
.DollarNames.rating.set <- function(x, pattern) {
  vals <- .codegen_dollar(x, pattern)
  
  cols <- colnames(x)
  if (!is.null(cols)) {
    vals <- c(cols, vals)
  }

  rows <- rownames(x)
  if (!is.null(rows)) {
    vals <- c(rows, vals)
  }
  
  vals
}

.codegen_dollar <- function(x, pattern = "") {
  attrs <- names(attributes(x))
  vals <- attrs[attrs %in% c("baserate", "kappa", "agreement")]
  return(vals)
}

#' @export
print.rating.set <- function(x, ...) {
  y <- unclass(x)

  for(a in additional_attributes) {
    attr(y, a) <- NULL
  }
  
  print(y, ...)
}