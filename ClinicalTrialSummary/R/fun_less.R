fun_less <- function(x, ...) {
  nx <- as.numeric(x)
  res <- ifelse(nx < 1e-04, paste("<0.0001"), x)
  return(res)
}

