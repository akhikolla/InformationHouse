#' Calculate the self exciting point process residuals
#' @param jtms A numerical vector containing the event times of the SEPP in ascending order
#' @param Tau A numerical scalar giving the censoring time, which should be greater than or equal to the event time
#' @param Inu A function that gives the integral of the baseline event rate function \eqn{\nu(t)} from 0 to the argument of the function
#' @param Ig A function that gives the integral of the excitation function \eqn{g(t)} from 0 to the argument of the function
#' @return A numerical vector containing the SEPP residuals \eqn{\Lambda(t_i)} where \eqn{\Lambda} is the cumulative intensity process. 
sepp.resid <- function(jtms,Tau,Inu,Ig){
  ttms <- c(jtms,Tau)
  res <- sapply(ttms,Inu)
  i <- 2;
  N <- length(ttms);
  while(i<=N){
    res[i] <- res[i]+sum(Ig(ttms[i]-ttms[1:(i-1)]))
    i <- i+1;
  }
  return(res)
}
