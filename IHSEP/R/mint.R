#' Sequence convolution

#' conv.seq calculates the convolution of two sequences
#'
#' @param x,y numeric vectors
#' @param real logical value, indicating whether the result should be real valued
#' @return a numeric vector of length equal to \code{length(x)+length(y)}, giving the convolution of \code{x} and \code{y}
#' @seealso \code{\link{convolve}}
#' @examples
#' conv.seq(1:5,1:4)
#' convolve(1:5,4:1,type="open")
conv.seq <- function(x,y,real=TRUE){
  m <- length(x); n <- length(y);
  x <- c(x,rep.int(0,n-1)); y <- c(y,rep.int(0,m-1));
  res <- fft(fft(x)*fft(y),inverse=TRUE)/(m+n-1);
  if(real)Re(res)else res;
}
#' Mean Intensity of the Self-Exciting Point Process With an Exponential Excitation Function

#' \code{h.fn.exp} calculates the mean intensity function \eqn{h(t)} which solves the integral equation \deqn{h(t)=\nu(t)+\int_0^t g(t-s)h(s)ds, t\geq 0}, where the excitation function is exponential: \eqn{g(t)= \gamma_1 e^{-\gamma_2 t}}.

#' @param x numerical scalar, at which the mean intensity \eqn{h} is to be evaluated
#' @param nu a function, which gives the baseline event rate
#' @param g.p a numeric vector of two elements giving the two parameters \eqn{\gamma_1,\gamma_2} of the exponential excitation function
#' @return a numric scalar which gives the value of the function \eqn{h} at \code{x}.
#' @seealso \code{\link{h.fn}}
#' @examples
#' nu <- function(x)200+100*cos(pi*x);
#' x <- 1:500/100;
#' y <- sapply(x,h.fn.exp,nu=nu,g.p=c(2,1));
#' h <- splinefun(x,y);
#' g <- function(x)2*exp(-x)
#' round(nu(x)+sapply(x,function(x)integrate(function(u)g(x-u)*h(u),0,x)$value) - y,5)
h.fn.exp <- function(x,nu,g.p=c(4,8)){
  nu(x)+integrate(function(y)nu(y)*g.p[1]*exp(-(g.p[2]-g.p[1])*(x-y)),
                  0,x)$value
}

#' Mean Intensity Function of the Self-Exciting Point Process

#' \code{h.fn} calculate the values of the mean intensity function of the self-exciting process with given baseline event rate and excitation function at a (fairly large) number of points. Values of the function at other points can be obtained by interpolation (e.g. spline interpolation).

#' @param nu a (vectorized) function specifying the baseline invent rate of the SEPP
#' @param g a (vectorized) function specifying the excitation function of the SEPP
#' @param N an integer giving the number of equal sized intervals to partition the domain into during the calculation The larger this value is, the more accurately the solution approxmates the truth, and the more time requred to evaluate.
#' @param to a numeric scalar, the end point of the estimation domain
#' @param abs.tol a numeric scalar specifying the absolute tolerance of error
#' @param maxit an integer specifying the maximal number of iterations allowed
#' @return a list with elelents, \code{x}: the vector of the points where \eqn{h} is evaluated; \code{y}: the vector of the corresponding \eqn{h} values; \code{nit}: the number of iterations used; \code{G.err}: the approximation error in \eqn{G}. 
#' @seealso \code{\link{h.fn1}}
#' @examples
#' nu <- function(x)(200+100*cos(pi*x))*(x>=0);
## g <- function(x)2*exp(-x)*(x>=0)
#' h.l <- h.fn1(nu=nu,g=g,from=0,to=5);
#' h1 <- splinefun(h.l$x,h.l$y);
#' round(nu(x)+sapply(x,function(x)integrate(function(u)g(x-u)*h(u),0,x)$value) - h1(x),5)

h.fn <- function(nu,g,N=2^12,to=1,abs.tol=1e-10,maxit=100){
  xx <- (1:N -1)/N;
  inc <- h <- nu(xx*to)*to;
  g.v <- g(xx*to)*to;
  nit <- 0;
  while(max(abs(inc <- conv.seq(inc,g.v)[1:N]/N))>(abs.tol*to) &&
        nit<=maxit){
    h <- h+inc; nit <- nit+1;
  }
  if(nit>maxit)warning("maxit reached!\n")
  list(x=xx*to,y=h/to,nit=nit,G.err=max(abs(inc))/to)
}

