# Here we have a problem with import
#' Parameter estimation procedure in continuous case.
#'
#' Parameter freq is preserved to allow for investigation of the inference procedure in high frequency case.
#' @import doParallel
#' @import foreach
#' @import ggplot2
#' @import stats
#' @import methods
#' @import plyr
#' @import grid
#' @import Rcpp
#' @importFrom stabledist rstable
#' @importFrom reshape2 melt
#' @importFrom Rdpack reprompt
#' @inheritParams path
#' @param t1,t2 real number such that  t2 > t1 > 0
#' @param p power
#' @param k increment order
#' @param path sample path of lfsm on which the inference is to be performed
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' alpha<-0.8; H<-0.8; sigma<-0.3
#' p<-0.3; k=3; t1=1; t2=2
#'
#' lfsm<-path(N=N,m=m,M=M,alpha=alpha,H=H,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#' ContinEstim(t1,t2,p,k,path=lfsm,freq='L')
#' @export
ContinEstim<-function(t1,t2,p,k,path,freq){

    H_est  <- H_hat(p,k,path)
    if(H_est<=0 | H_est>=1) return('H estimate is not in (0,1)')

    alpha_est <- alpha_hat(t1,t2,k,path,H_est,freq)
    if(alpha_est<=0 | alpha_est>=2) return('alpha estimate is not in (0,2)')

    sigma_est<-tryCatch(
        sigma_hat(t1,k,path,alpha_est,H_est,freq),
        error=function(c) 'Something went wrong while computing sigma_hat')

    list(alpha=alpha_est,H=H_est,sigma=sigma_est)

}



# !!!! In theorem 4.1 frequency is always low
#' Low frequency estimation procedure for lfsm.
#'
#' General estimation procedure for low frequency case when 1/alpha is not a natural number.
#' @inheritParams ContinEstim
#' @inheritParams path
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' sigma<-0.3
#' p<-0.3; k=3; t1=1; t2=2
#'
#' #### Continuous case
#' lfsm<-path(N=N,m=m,M=M,alpha=1.8,H=0.8,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#'
#' GenLowEstim(t1=t1,t2=t2,p=p,path=lfsm,freq="L")
#'
#' #### H-1/alpha<0 case
#' lfsm<-path(N=N,m=m,M=M,alpha=0.8,H=0.8,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#'
#' GenLowEstim(t1=t1,t2=t2,p=p,path=lfsm,freq="L")
#'
#' #### The procedure works also for high frequency case
#' lfsm<-path(N=N,m=m,M=M,alpha=1.8,H=0.8,
#'            sigma=sigma,freq='H',disable_X=FALSE,seed=3)$lfsm
#'
#' GenLowEstim(t1=t1,t2=t2,p=p,path=lfsm,freq="H")
#' @export
GenLowEstim<-function(t1,t2,p,path,freq='L'){

    H_0<-H_hat(p=-p,k=1,path=path)

    alpha_0<-alpha_hat(t1=t1,t2=t2,k=1,path=path,H=H_0,freq=freq)
    if(alpha_0<=0) return('alpha_0 estimate is less or equal to 0') else {

    k_new<-2+floor(alpha_0^(-1))

    if(2*k_new>length(path)) return('2k estimate is larger than the path length. Can\'t use R_hl')

    H_est<-H_hat(p=-p,k=k_new,path=path)
        if(is.nan(H_est)) return('H cannot be estimated')
    alpha_est<-alpha_hat(t1=t1,t2=t2,k=k_new,path=path,H=H_est,freq=freq)
        if(is.nan(alpha_est)) return(list(alpha=alpha_est,H=H_est,sigma=NaN))

    if(alpha_est<=0 | alpha_est>=2) return('alpha estimate is not in (0,2)') else {

        if(H_est<=0 | H_est>=1) return('H estimate is not in (0,1)') else {

                sigma_est<-tryCatch(
                    sigma_hat(t1=t1,k=k_new,path=path,alpha=alpha_est,H=H_est,freq=freq),
                    error=function(c) 'H*alpha-1 estimate must be too close to -1 to compute sigma estimate')

                ifelse(is.character(sigma_est),
                       return(sigma_est),return(list(alpha=alpha_est,H=H_est,sigma=sigma_est)))

        }
    }
    }
}


#' Inverse alpha estimator
#'
#' A function from a general estimation procedure which is defined as m^p_-p'_k /m^p'_-p_k, originally proposed in [13].
#' @references \insertRef{MOP18}{rlfsm}
#' @inheritParams path
#' @inheritParams GenHighEstim
#' @export
phi_of_alpha<-function(p,p_prime,alpha){

    numerator<-((2/alpha)^(p-p_prime))*((a_p(-p))^p_prime)*((gamma(p_prime/alpha))^p)

    denominator<-((a_p(-p_prime))^p)*((gamma(p/alpha))^p_prime)

    numerator/denominator
}


f_p<-function(x,p) (abs(x))^p

#' High frequency estimation procedure for lfsm.
#'
#' General estimation procedure for high frequency case when 1/alpha is not a natural number.
#' "Unnecessary" parameter freq is preserved to allow for investigation of the inference procedure in low frequency case
#' @inheritParams ContinEstim
#' @inheritParams path
#' @param p_prime power
#' @param low_bound positive real number
#' @param up_bound positive real number
#' @details
#' In this algorithm the preliminary estimate of alpha is found via using \code{\link{uniroot}} function. The latter is
#' given the lower and the upper bounds for alpha via low_bound and up_bound parameters. It is not possible to
#' pass 0 as the lower bound because there are numerical limitations on the alpha estimate, caused by the
#' length of the sample path and by numerical errors. p and p_prime must belong to the interval (0,1/2) (in the notation kept in rlfsm package)
#' The two powers cannot be equal.
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' sigma<-0.3
#' p<-0.2; p_prime<-0.4
#'
#' #### Continuous case
#' lfsm<-path(N=N,m=m,M=M,alpha=1.8,H=0.8,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#'
#' GenHighEstim(p=p,p_prime=p_prime,path=lfsm,freq="H")
#'
#' #### H-1/alpha<0 case
#' lfsm<-path(N=N,m=m,M=M,alpha=0.8,H=0.8,
#'            sigma=sigma,freq='H',disable_X=FALSE,seed=3)$lfsm
#'
#' GenHighEstim(p=p,p_prime=p_prime,path=lfsm,freq="H")
#'
#' @export
GenHighEstim<-function(p,p_prime,path,freq,low_bound=0.01,up_bound=4){

    # p from (0, 1/2)
    # p =/= p'

    # Find a preliminary estimate for H
    H_0<-H_hat(p=-p,k=1,path=path)
    if(H_0<=0 | H_0>=1) return('H estimate is not in (0,1)')

    m_numer_est<-sf(path=path,f=f_p,k=1,r=1,H=H_0,freq=freq, p=-p_prime)
    m_denom_est<-sf(path=path,f=f_p,k=1,r=1,H=H_0,freq=freq, p=-p)

    # Find a preliminary estimate for alpha
    Funct<-function(x) (m_numer_est^p)/(m_denom_est^p_prime) - phi_of_alpha(p=p,p_prime=p_prime,alpha=x)

    alpha_0<-try(uniroot(Funct, interval= c(low_bound,up_bound))$root, silent=T)
    if(is(alpha_0,"try-error")) return("alpha_0 doesnt belong to the interval")

    # Find a preliminary estimate for k
    k_est<-2+floor(alpha_0^(-1))

    ### final estimates

    # Find an estimate for H
    H_est<-H_hat(p=-p,k=k_est,path=path)
    if(H_est<=0 | H_est>=1) return('H estimate is not in (0,1)')

    m_numer_est<-sf(path=path,f=f_p,k=k_est,r=1,H=H_0,freq=freq, p=-p_prime)
    m_denom_est<-sf(path=path,f=f_p,k=k_est,r=1,H=H_0,freq=freq, p=-p)

    # Find an estimate for alpha
    Funct<-function(x) (m_numer_est^p)/(m_denom_est^p_prime) - phi_of_alpha(p=p,p_prime=p_prime,alpha=x)

    alpha_est<-try(uniroot(Funct, interval= c(low_bound,up_bound))$root, silent=T)
    if(is(alpha_est,"try-error")) return('alpha_est doesn\'t belong to the interval')

    sigma_est<-try(((alpha_est*a_p(-p)*m_denom_est/2/gamma(p/alpha_est))^(-1/p))/(Norm_alpha(h_kr,alpha=alpha_est,k=k_est,r=1,H=H_est,l=0)$result), silent=T)
    if(is(sigma_est,"try-error")) return('Problems with computing sigma estimate')

    list(alpha=alpha_est, H=H_est, sigma=sigma_est)

}

