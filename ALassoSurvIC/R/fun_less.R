fun_less <- function(x, ...) {
  nx <- as.numeric(x)
  ax <- c()
  for (i in 1:length(x)) {
    ax[i] <- ifelse(nx[i] < 1e-04, paste("<0.0001"), x[i])
  }
  return(ax)
}

