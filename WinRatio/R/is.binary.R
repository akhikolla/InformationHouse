is.binary <- function(x){
  x <- stats::na.omit(x)
  if (is.numeric(x)){
    all(x %in% c(0, 1))
  } else {FALSE}
}