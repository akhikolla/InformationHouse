##' @title Full Significance Bayesian Testing
##' @description Performs Full Significance Bayesian Testing (FSBT) for univariate sharp
##' null hypothesis based on a posterior sample. The marginal posterior density is obtained by kernel density estimation from \code{sim.sample}.
##' @usage kdeFSBT(H0, sim.sample)
##'
##' @param H0 a scalar value under the null hypothesis.
##' @param sim.sample a sample from the marginal posterior distribution.
##'
##' @return double
##'
##' @examples
##'
##' x <-  rnorm(1000, 0, 1)
##' kdeFSBT(-1, x)
##'
##' @references
##' Pereira, C. A. d. B., Stern, J. M. and Wechsler, S. (2008) Can a significance test be genuinely Bayesian? \emph{Bayesian Analysis} \bold{3}, 79-100.
##'
##' @export
kdeFSBT <- function(H0, sim.sample){
  if(missing(H0) || missing(sim.sample) || !is.numeric(H0))
    stop("\'H0\' and \'sample\' must be a scalar and a vector respectively.")

   EV <- NA

   if(!(H0%in%range(sim.sample))) {
     warning("\'H0\' is outside the range of \'sim.sample\'")
     EV = 0.0
   } else {
     try({
       den <- splinefun(density(sim.sample))
       d.H0 <- den(H0)
       cent_dens <- median(sim.sample)

       if(cent_dens < H0){
         param.H0.other <- uniroot(function(x)
           den(x)-d.H0, interval=c(min(sim.sample), cent_dens))$root
         EV <- integrate(function(x) den(x), min(sim.sample), param.H0.other)$value +
           integrate(function(x) den(x), H0, max(sim.sample))$value
       } else {
         param.H0.other <- uniroot(function(x)
           den(x)-d.H0, interval=c(cent_dens, max(sim.sample)))$root
         EV <- integrate(function(x) den(x), param.H0.other, max(sim.sample))$value +
           integrate(function(x) den(x), min(sim.sample), H0)$value
       }
     })
   }
   return(EV)
}