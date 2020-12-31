# 1<j<d, d- dimensionality of W vector
#
#### New version of h_kr ####
#' Function h_kr
#'
#' Function \eqn{h_{k,r}: R \to R} is given by \deqn{h_{k,r} (x) = \sum_{j=0}^k (-1)^j {{k}\choose{j}} (x-rj)_+^{H-1/\alpha}, \ \ \ x\in R}
#' @param x real number
#' @param l a value by which we shift x. Is used for computing function f_.+l and is passed to integrate function.
#' @inheritParams increment
#' @inheritParams path
#' @examples
#' #### Plot h_kr ####
#' s<-seq(0,10, by=0.01)
#' h_val<-sapply(s,h_kr, k=5, r=1, H=0.3, alpha=1)
#' plot(s,h_val)
#' @references \insertRef{MOP18}{rlfsm}
#' @export
#'
#'

h_kr<- function(k,r,x,H,alpha,l=0)
{
    S<-0
    for (j in 0:k)
    {
	      S<-S + ifelse(x+l-r*j>0,(-1)^j * choose(k,j) *(x+l-r*j)^(H-1/alpha),0)
    }
    S
}



#### Norm of an arbitrary function ####

#' Alpha-norm of an arbitrary function
#' @param fun a function to compute a norm
#' @param ... a set of parameters to pass to integrate
#' @inheritParams path
#' @details
#' fun must accept a vector of values for evaluation. See ?integrate for further details.
#' Most problems with this function appear because of rather high precision. Try to tune rel.tol parameter first.
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#' Norm_alpha(h_kr,alpha=1.8,k=2,r=1,H=0.8,l=4)
#' @export
#'
Norm_alpha<-function(fun,alpha,...) {

    if (sum(formalArgs(fun)=='alpha')) {
	      int<-function(ff,alpha,...){(abs(ff(alpha=alpha,...)))^alpha}
    } else {
	      int<-function(ff,alpha,...){(abs(ff(...)))^alpha}
    }
    # !!!!!!!!!!!!!!!!!!
    res<-integrate(f=int, lower=-Inf, upper=Inf, ff=fun, alpha=alpha, subdivisions=20000, ... ) #, subdivisions=20000, rel.tol = .Machine$double.eps^0.05)
    list(result=(res$value)^(1/alpha),abs.error=res$abs.error)
}


