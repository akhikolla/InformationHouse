# !!!!!! Figure out where is it from
# Write docs and assign normal names in order to make it
# obvious which theorem the functions are from


#'
#' Function a_p.
#'
#' Computes the corresponding function value from Mazur et al. 2018.
#' @param p power, real number from (-1,1)
#' @references \insertRef{MOP18}{rlfsm}
#' @export
a_p <- function(p){
  if(p<1 & p>0){

    integrand<-function(y) {(1-exp(1i*y))*(abs(y))^(-1-p)}
    elliptic::myintegrate(integrand,lower=-Inf, upper=Inf,subdivisions=3000, rel.tol=.Machine$double.eps^.1,stop.on.error = TRUE, keep.xy = FALSE, aux = NULL)

  } else {

    if(p<0 & p>-1) sqrt(2*pi)*gamma(-p/2)/(2^(p+1/2))/gamma((p+1)/2) else stop("p must be in the range (-1,0) & (0,1)")

  }
  }

#### U_gh ####
#' alpha-norm of u*g + v*h.
#' @inheritParams U_ghuv
#' @param ... additional parameters to pass to Norm_alpha
#' @examples
#' g<-function(x) exp(-x^2)
#' h<-function(x) exp(-abs(x))
#' U_gh(g=g, h=h, u=4, v=3, alpha=1.7)
#' @export
U_gh<-function(g,h,u,v,...){

  ffun<-function(s,g,h,u,v) u*g(s)+v*h(s)
  l<-Norm_alpha(fun=ffun, g=g, h=h, u=u, v=v, ...)
  l
}


#### U_g ####
#' alpha norm of u*g
#'
#' @inheritParams U_ghuv
#' @param ... additional parameters to pass to Norm_alpha
#' @examples
#' g<-function(x) exp(-x^2)
#' g<-function(x) exp(-abs(x))
#' U_g(g=g,u=4,alpha=1.7)
#' @export
U_g<-function(g,u,...){

  ffun<-function(s, g, u) u*g(s)
  l<-Norm_alpha(fun=ffun, g=g, u=u, ...)
  l
}


#### U_ghuv ####
# We need to pass parameters to integrate !!!!!!
#' A dependence structure of 2 random variables.
#'
#' It is used when random variables do not have finite second moments, and thus, the covariance matrix is not defined.
#' For \eqn{X= \int_{\R} g_s dL_s} and \eqn{Y= \int_{\R} h_s dL_s} with \eqn{\| g \|_{\alpha}, \| h\|_{\alpha}< \infty}. Then the measure of dependence is given by \eqn{U_{g,h}: \R^2 to \R} via
#' \deqn{U_{g,h} (u,v)=\exp(- \sigma^{\alpha}{\| ug +vh \|_{\alpha}}^{\alpha} ) - \exp(- \sigma^{\alpha} ({\| ug \|_{\alpha}}^{\alpha} + {\| vh \|_{\alpha}}^{\alpha}))}
#' @inheritParams path
#' @param g,h functions \eqn{g,h: \R to \R} with finite alpha-norm (see \code{\link{Norm_alpha}}).
#' @param v,u real numbers
#' @param ... additional parameters to pass to U_gh and U_g
#' @examples
#' g<-function(x) exp(-x^2)
#' h<-function(x) exp(-abs(x))
#' U_ghuv(alpha=1.5, sigma=1, g=g, h=h, u=10, v=15,
#' rel.tol = .Machine$double.eps^0.25, abs.tol=1e-11)
#' @export
U_ghuv<-function(alpha, sigma, g, h, u, v, ...){

  R1<-U_gh(g=g, h=h, u=u, v=v, alpha=alpha, ...)
  R2<-U_g(g=g,u=u,alpha=alpha, ...)
  R3<-U_g(g=h,u=v,alpha=alpha, ...)

  r1<-R1$result ; er1<-R1$abs.error
  r2<-R2$result ; er2<-R2$abs.error
  r3<-R3$result ; er3<-R3$abs.error

  r<-exp(-(sigma*r1)^alpha) -
    exp(-sigma^alpha*((r2)^alpha+(r3)^alpha))

  ap_error<- alpha*sigma^alpha*exp(-(sigma*r1)^alpha)*r1^(alpha-1)*er1+
                    alpha*sigma^alpha*exp(-(sigma^alpha*(r2^alpha+r3^alpha)))*r2^(alpha-1)*er2+
                    alpha*sigma^alpha*exp(-(sigma^alpha*(r2^alpha+r3^alpha)))*r3^(alpha-1)*er3
  list(result=r,abs_aprox_error=ap_error)
}



#### Theta ####
# How to compute an error of an integral of an integral ? !!!!!!!!!!!!!!!!!!!!!!
#' Function of the form
#' \deqn{\theta(g,h)_{p} = a_p^{-2} \int_{\R^2}  |xy|^{-1-p}U_{g,h}(x,y) dxdy}
#' @export
#' @inheritParams U_ghuv
#' @inheritParams path
#' @inheritParams a_p
#' @references \insertRef{MOP18}{rlfsm}
theta<-function(p,alpha,sigma,g,h){

    ap<-a_p(p)

    f_xy<-function(x,y,p,g,h,...){
      (abs(x*y))^(-1-p)*U_ghuv(g=g, h=h, u=x, v=y,...)
    }
    # f_xy(x=4,y=1,p=0.1,g,h,alpha=1,sigma=1)

    f_xy_vec<-Vectorize(FUN=f_xy, vectorize.args = 'x')
    # Fvec(x=c(1,7,0,4),y=1,p=0.1,g,h,alpha=1,sigma=1)

    int_over_x<-function(y,p,g,h,...){
       integrate(f_xy_vec, lower=-Inf, upper=Inf, y=y, p=p, g=g, h=h,
                 rel.tol=.Machine$double.eps^0.1, subdivisions = 300, ...)$value
    }
    # int_over_x(y=1, p=0.5, g=g, h=h, alpha=1.5, sigma=1)

    int_over_x_vec<-Vectorize(FUN=int_over_x, vectorize.args = 'y')

    int_over_y<-function(p,g,h,...){
      integrate(int_over_x_vec, lower=-Inf, upper=Inf, p=p, g=g, h=h,
                rel.tol=.Machine$double.eps^0.1, subdivisions = 300, ...)$value
    }
    # int_over_y(p=0.5,g=g,h=h,alpha=1.5, sigma=1)

    ap^(-2)*int_over_y(p=p, g=g, h=h, alpha=alpha, sigma=sigma)
}
# theta(p=0.5, alpha=1.5, sigma=1, g=g, h=h)


#### cov_theta_hh ####
# r_j=j ; j={1,2}
# !!!! Don't forget to test cov_theta_hh !!!!
# cov_theta_hh is cov(W¹,W¹)
#' @inheritParams path
#' @inheritParams U_ghuv
cov_theta_hh<-function(p,j,j_p,H,alpha,sigma,l_max,k=2){

    summ<-0
    # the first h in 2.19(i)
    h_0<-function(x) h_kr(k=k,H=H,alpha=alpha,l=0,r=j, x=x)
    # h_0(9)


    for(l in -l_max:l_max){
        # the second h in 2.19(i)
        h_l<-function(x) h_kr(x=x,k=k,H=H,alpha=alpha,l=l,r=j_p)
        # h_l(19)
        h_0_vec<-Vectorize(FUN=h_0, vectorize.args = 'x')
        h_l_vec<-Vectorize(FUN=h_l, vectorize.args = 'x')
        summ<-summ+theta(p=p, alpha=alpha, sigma=sigma, g=h_0_vec, h=h_l_vec)
    }
    summ
}
# cov_theta_hh(p=0.4,j=1,j_p=1,H=0.8,alpha=1.8,sigma=1,l_max=1,k=2)


#### cov(W²,W²)####
### And here we have a problem ###
# (2.19ii)
# j1<-1
# j2<-1
# t<-c(1,2)
# l_max<-1000; H<-0.95; alpha<-1.9; sigma<-1

# summ<-0
# summ_err<-0
# lsum<-vector()
# h_1<-function(x) h_kr(x=x,k=k,H=H,alpha=alpha,l=0,r=1)

# for(l in 900:l_max){

  # h_l1<-function(x) h_kr(x=x,k=k,H=H,alpha=alpha,l=l,r=1)
  # Part1<-U_ghuv(alpha=alpha, sigma=sigma, g=h_1, h=h_l1, u=t[j1], v=t[j2],
    #             subdivisions=2000, rel.tol = .Machine$double.eps^0.05, abs.tol=1e-11)

  # h_l2<-function(x) -h_l1(x)
  # Part2<-U_ghuv(alpha=alpha, sigma=sigma, g=h_1, h=h_l2, u=t[j1], v=t[j2],
    #             subdivisions=2000, rel.tol = .Machine$double.eps^0.05, abs.tol=1e-11)

  # summ<-summ+Part1$result+Part2$result
  # summ_err<-summ_err+Part1$abs_aprox_error+Part2$abs_aprox_error
  # lsum<-rbind(lsum,c(Part1,Part2))
  # list(result=summ,abs_aprox_error=summ_err)
# }


# l<-21
# h_l1<-function(x) h_kr(x=x,k=k,H=H,alpha=alpha,l=l,r=1)
# U_ghuv(alpha=alpha, sigma=sigma, g=h_1, h=h_l1, u=t[j1], v=t[j2])

