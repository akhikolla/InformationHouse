limit <- function(x) {
  y <- ifelse((x%%2 == 0),x/2,(x-1)/2)
  return(y)
}
