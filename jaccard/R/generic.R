getp <- function(lr,lr0,ties.method = "average") {
  # Adapted from J.D.Storey, see qvalue::empPvals
  # Get resampled p-values, pulling across variables (e.g., genes)
  # lr: observed statistics
  # lr0: null statistics (i.e. from resampled residuals)

  m = length(lr)
  m0 = length(lr0)
  v = c(rep(TRUE,m),rep(FALSE,m0))
  v = v[rev(order(c(lr,lr0)))]
  u = 1:length(v)
  w = 1:m
  p = ((u[v==TRUE]-w))/(length(lr0))
  #to account for the extreme observed values
  #p = ((u[v==TRUE]-w)+1)/(length(lr0)+2)
  p = p[rank(-lr, ties.method = ties.method)]
  p = pmax(p, 1/m0)
  return(p)
}
