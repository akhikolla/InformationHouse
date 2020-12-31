
R_xlen_t <- function(x) {
  if (max(x, na.rm = TRUE) <= 2147483647) {
    x <- as.integer(x)
  }
  x
}
