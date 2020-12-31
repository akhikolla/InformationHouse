
#### Higher order increments ####

#' Higher order increments
#'
#' Difference of the kth order. Defined as following:
#' \deqn{\Delta_{i,k}^{n,r} X:= \sum_{j=0}^k (-1)^j{{k}\choose{j}}X_{(i-rj)/n}, i\geq rk.}
#' Index i here is a coordinate in terms of point_num. Although R uses vector indexes
#' that start from 1, increment has i varying from 0 to N, so that a
#' vector has a length N+1. It is done in order to comply with the notation of the paper.
#' This function doesn't allow for choosing frequency n. The frequency is determined by the
#' path supplied, thus n equals to either the length of the path in high frequency setting
#' or 1 in low frequency setting. increment() gives increments at certain point passed as i,
#' which is a vector here. increments() computes high order increments for the whole sample
#' path. The first function evaluates the formula above, while the second one uses structure
#' diff(diff(...)) because the formula is slower at higher k.
#'
#' @param r difference step, a natural number
#' @param i index of the point at which the increment is to be computed, a natural number.
#' @param k order of the increment, a natural number
#' @param path sample path for which a kth order increment is computed
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#'
#' m<-45; M<-60; N<-2^10-M
#' alpha<-0.8; H<-0.8; sigma<-0.3
#'
#' lfsm<-path(N=N,m=m,M=M,alpha=alpha,H=H,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#'
#' tryCatch(
#' increment(r=1,i=length(lfsm),k=length(lfsm)+100,path=lfsm),
#' error=function(c) 'An error occures when k is larger then the length of the sample path')
#'
#' increment(r=3,i=50,k=3,path=lfsm)
#'
#'
# The functions produce the same result on the whole path
#' path=c(1,4,3,6,8,5,3,5,8,5,1,8,6)
#'
#' r=2; k=3
#' n <- length(path) - 1
#' DeltaX = increment(seq(r*k, n), path = path, k = k, r = r)
#' DeltaX == increments(k=k,r=r,path)
#' @export
increment<-function(r,i,k,path) {

    if (sum(i<r*k)) stop("i must be greater or equal to r*k") else
    {
	    S<-0
	    for (j in 0:k)
	    {
		    S<-S + (-1)^j * choose(k,j) * path[i+1-r*j]
	    }

	S
    }
}


#' @rdname increment
#' @export
increments<-function(k,r,path) {

        for (j in 1:k)
        {
            n<-length(path)
            path<-(path[(r+1):n] - path[1:(n-r)])
        }

    path
}

#' Statistic V
#'
#' Statistic of the form
#' \deqn{V_{\textnormal{high}}(f; k,r)_n:=\frac{1}{n}\sum_{i=rk}^n f\left( n^H \Delta_{i,k}^{n,r} X \right),  }
#' \deqn{V_{\textnormal{low}}(f; k,r)_n :=\frac{1}{n}\sum_{i=rk}^n f\left( \Delta_{i,k}^{r} X \right)}
#'
#' @param path sample path for which the statistic is to be calculated.
#' @param f function applied to high order increments.
#' @param k order of the increments.
#' @param r step of high order increments.
#' @param H Hurst parameter.
#' @param freq frequency.
#' @param ... parameters to pass to function f
#' @details Hurst parameter is required only in high frequency case. In the low frequency, there is no need to assign H a value because it will not be evaluated.
#' @references \insertRef{MOP18}{rlfsm}
#' @seealso \code{\link{phi}} computes V statistic with f(.)=cos(t.)
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' freq='L'
#' r=1; k=2; p=0.4
#' S<-(1:20)*100
#'
#' path_lfsm<-function(...){
#'
#'     List<-path(...)
#'     List$lfsm
#'
#' }
#'
#' Pths<-lapply(X=S,FUN=path_lfsm,
#'              m=m, M=M, alpha=alpha, sigma=sigma, H=H,
#'              freq=freq, disable_X = FALSE,
#'              levy_increments = NULL, seed = NULL)
#'
#' f_phi<-function(t,x) cos(t*x)
#' f_pow<-function(x,p) (abs(x))^p
#'
#' V_cos<-sapply(Pths,FUN=sf,f=f_phi,k=k,r=r,H=H,freq=freq,t=1)
#' ex<-exp(-(abs(sigma*Norm_alpha(h_kr,alpha=alpha,k=k,r=r,H=H,l=0)$result)^alpha))
#'
#'  # Illustration of the law of large numbers for phi:
#' plot(y=V_cos, x=S, ylim = c(0,max(V_cos)+0.1))
#' abline(h=ex, col='brown')
#'
#' # Illustration of the law of large numbers for power functions:
#' Mpk<-m_pk(k=k, p=p, alpha=alpha, H=H, sigma=sigma)
#'
#' sf_mod<-function(Xpath,...) {
#'
#'     Path<-unlist(Xpath)
#'     sf(path=Path,...)
#' }
#'
#' V_pow<-sapply(Pths,FUN=sf_mod,f=f_pow,k=k,r=r,H=H,freq=freq,p=p)
#' plot(y=V_pow, x=S, ylim = c(0,max(V_pow)+0.1))
#' abline(h=Mpk, col='brown')
#'
#############################################################################
#'
#' @export
#### V high / low ####

sf<-function(path,f,k,r,H,freq,...){

    n<-length(path)-1 # -1 because the scalling factor is 1/(n-1)
    v1<-increments(k=k,r=r,path=path)

    if(freq=='L') v2<-sapply(v1,FUN=f,...) else{
        if(freq=='H') v2<-sapply((n^H)*v1,FUN=f,...) else{
            stop('Parameter freq could take either "H" or "L" values')
        }
    }
    sum(v2)/(length(v2))
}

# Deterministic version from Theorem 4.5
#' m(-p,k)
#'
#' defined as \eqn{m_{p,k} := E[|\Delta_{k,k} X|^p]} for positive powers.  When p is negative (-p is positive) the equality does not hold.
#'
#' The following identity is used for computations: \deqn{m_{-p,k} =  \frac{(\sigma \|h_k\|_{\alpha})^{-p}}{a_{-p}} \int_{\R} \exp(-|y|^{\alpha}) |y|^{-1+p} dy = \frac{2(\sigma \|h_k\|_{\alpha})^{-p}}{\alpha a_{-p}} \Gamma(p/\alpha)}
#' @inheritParams path
#' @inheritParams ContinEstim
#' @param p a positive number
#' @references \insertRef{MOP18}{rlfsm}
#' @export
m_pk<-function(k,p,alpha,H,sigma){

    h_k<-Norm_alpha(h_kr,alpha=alpha,k=k,r=1,H=H,l=0)$result
    2*((sigma*h_k)^(-p))*gamma(p/alpha)/alpha/a_p(-p)

}


#### R_high/low function ####
#' R high /low
#'
#' Defined as
#' \deqn{R_{\textnormal{high}} (p,k)_n := \frac{\sum_{i=2k}^n \left| \Delta_{i,k}^{n,2} X \right|^p}
#' {\sum_{i=k}^n \left| \Delta_{i,k}^{n,1} X \right|^p}, \qquad}
#' \deqn{R_{\textnormal{low}} (p,k)_n := \frac{\sum_{i=2k}^n \left| \Delta_{i,k}^{2} X \right|^p}
#' {\sum_{i=k}^n \left| \Delta_{i,k}^{1} X \right|^p}}
#'
#' The computation procedure for high- and low frequency cases is the same, since there is no way to control frequency given a sample path.
#' @inheritParams ContinEstim
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' alpha<-0.8; H<-0.8; sigma<-0.3
#' p<-0.3; k=3
#'
#' lfsm<-path(N=N,m=m,M=M,alpha=alpha,H=H,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#' R_hl(p=p,k=k,path=lfsm)
#' @export
R_hl<-function(p,k,path){

    n=length(path)-1
    if (n<2*k) stop('X has too few points for the chosen k')

    S_up<-sum((abs(increments(r=2, k=k, path=path)))^p)
    S_low<-sum((abs(increments(r=1, k=k, path=path)))^p)
    S_up/S_low
}



#### Statistical estimator of H in high/low frequency setting ####
#' Statistical estimator of H in high/low frequency setting
#'
#' The statistic is defined as
#' \deqn{\widehat{H}_{\textnormal{high}} (p,k)_n:= \frac 1p \log_2  R_{\textnormal{high}} (p,k)_n,
#' \qquad \widehat{H}_{\textnormal{low}} (p,k)_n:= \frac 1p \log_2  R_{\textnormal{low}} (p,k)_n}
#' @inheritParams ContinEstim
#' @references \insertRef{MOP18}{rlfsm}
#' @export
H_hat<-function(p,k,path){

    1/p*log2(R_hl(p,k,path))

}


#### phi ####
#' Phi
#'
#' Defined as
#' \deqn{\varphi_{\textnormal{high}}(t; H,k)_n := V_{\textnormal{high}}(\psi_t; k)_n \qquad \textnormal{and}
#' \qquad \varphi_{\textnormal{low}}(t; k)_n := V_{\textnormal{low}}(\psi_t; k)_n},
#' where \eqn{\psi_t(x):=cos(tx)}
#' @details Hurst parameter is required only in high frequency case. In the low frequency, there is no need to assign H a value because it will not be evaluated.
#'
#' @inheritParams path
#' @inheritParams ContinEstim
#' @param t positive real number
#' @references \insertRef{MOP18}{rlfsm}
#' @export
phi<-function(t,k,path,H,freq){

    N<-length(path)-1

    if(freq=="L") {
        FreqFactor<-1
    } else {

        if(freq=="H"){
            FreqFactor<-N^H
        } else {
            stop("parameter freq can take only two values- 'H' and 'L' ")
        }

    }

    v1<-increments(k=k,r=1,path=path)
    sum(cos(t*FreqFactor*v1))/(N-k)
}

# not for export
phi_low<-function(t,k,path){

    N<-length(path)-1
    v1<-increments(k=k,r=1,path=path)
    sum(cos(t*v1)/(N-k))

}

## alpha_hat estimation
# |phi| is used
#' Statistical estimator for alpha
#'
#' Defined for the two frequencies as
#' \deqn{\widehat \alpha_{high} := \frac{\log | \log \varphi_{high} (t_2; \widehat H_{high} (p,k)_n, k)_n|  -  \log | \log \varphi_{high} (t_1; \widehat H_{high} (p,k)_n, k)_n|}{\log t_2 - \log t_1}}
#' \deqn{\widehat \alpha_{low} := \frac{\log | \log \varphi_{low} (t_2;k)_n|  -  \log | \log \varphi_{low} (t_1; k)_n|}{\log t_2 - \log t_1}}
#' @inheritParams path
#' @inheritParams ContinEstim
#' @details The function triggers function \code{\link{phi}}, thus Hurst parameter is required only in high frequency case. In the low frequency, there is no need to assign H a value because it will not be evaluated.
#' @references \insertRef{MOP18}{rlfsm}
#' @examples
#'
#' m<-45; M<-60; N<-2^14-M
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' freq='H'
#' r=1; k=2; p=0.4; t1=1; t2=2
#'
#' # Estimating alpha in the high frequency case
#' # using preliminary estimation of H
#' lfsm<-path(N=N,m=m,M=M,alpha=alpha,H=H,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=3)$lfsm
#'
#' H_est<-H_hat(p=p,k=k,path=lfsm)
#' H_est
#' alpha_est<-alpha_hat(t1=t1,t2=t2,k=k,path=lfsm,H=H_est,freq=freq)
#' alpha_est
#'
#' @export
alpha_hat<-function(t1,t2,k,path,H,freq){

    V1<-abs(phi(t1,k,path,H,freq))
    V2<-abs(phi(t2,k,path,H,freq))
    (log(abs(log(V2)))-log(abs(log(V1))))/(log(t2)-log(t1))

}


## sigma_hat estimation
# H and alpha should be taken from estimates
#' Statistical estimator for sigma
#' @inheritParams path
#' @inheritParams ContinEstim
#' @examples
#' m<-45; M<-60; N<-2^14-M
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' freq='H'
#' r=1; k=2; p=0.4; t1=1; t2=2
#'
#' # Reproducing the work of ContinEstim
#' # in high frequency case
#' lfsm<-path(N=N,m=m,M=M,alpha=alpha,H=H,
#'            sigma=sigma,freq='L',disable_X=FALSE,seed=1)$lfsm
#'
#' H_est<-H_hat(p=p,k=k,path=lfsm)
#' H_est
#' alpha_est<-alpha_hat(t1=t1,t2=t2,k=k,path=lfsm,H=H_est,freq=freq)
#' alpha_est
#'
#' sigma_est<-tryCatch(
#'                     sigma_hat(t1=t1,k=k,path=lfsm,
#'                     alpha=alpha_est,H=H_est,freq=freq),
#'                     error=function(c) 'Impossible to compute sigma_est')
#' sigma_est
#' @export
sigma_hat<-function(t1,k,path,alpha,H,freq){

    V1<-abs(phi(t=t1,k=k,path=path,H=H,freq=freq))
    abs(((-log(V1))^(1/alpha))/(t1*Norm_alpha(h_kr,alpha=alpha,k=k,r=1,H=H,l=0)$result))
    # abs is not from MOP18

}


#### Generate sample paths ####
# LofF- a list of functions to load on each node
LofF<-unclass(lsf.str())[1:length(lsf.str())]

#' Generator of a set of lfsm paths.
#'
#' It is essentially a wrapper for \code{\link{path}} generator, which exploits the latest to create a matrix with paths in its columns.
#' @param N_var number of lfsm paths to generate
#' @param ... arguments to pass to path
#' @param parallel a TRUE/FALSE flag which determines if the paths will be created in parallel or sequentially
#' @param seed_list a numerical vector of seeds to pass to \code{\link{path}}
#' @seealso \code{\link{path}}
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' freq='L'
#' r=1; k=2; p=0.4
#'
#' Y<-paths(N_var=10,parallel=TRUE,N=N,m=m,M=M,
#'          alpha=alpha,H=H,sigma=sigma,freq='L',
#'          disable_X=FALSE,levy_increments=NULL)
#'
#' Hs<-apply(Y,MARGIN=2,H_hat,p=p,k=k)
#' hist(Hs)
#'
#' @export
#'
paths<-function(N_var,parallel,seed_list=rep(x=NULL, times = N_var),...){
    ind_par<-NULL # avoids NOTEs when being built
    if(length(seed_list)!=N_var & !is.null(seed_list)) {stop('Seed_list must be equal to N_var')} else
    {
        # Switcher between parallel and consecutive executions
        `%switch_do%` <- ifelse(parallel, `%dopar%`, `%do%`)
        r<-matrix()
        r<-foreach::foreach (ind_par = 1:N_var, .combine = cbind, .packages='stabledist', .export = LofF) %switch_do% {
	        X<-path(seed=seed_list[ind_par],...)
	        X$lfsm
        }
        r
    }
}


#### Apply statistics to sample paths ####
#' Retrieve statistics(bias, variance) of estimators based on a set of paths
#' @param paths real-valued matrix representing sample paths of the stochastic process being studied
#' @param true_val true value of the estimated parameter
#' @param Est estimator (i.e. H_hat)
#' @param ... parameters to pass to Est
#' @examples
#' m<-45; M<-60; N<-2^10-M
#' alpha<-1.8; H<-0.8; sigma<-0.3
#' freq='L';t1=1; t2=2
#' r=1; k=2; p=0.4
#'
#' Y<-paths(N_var=10,parallel=TRUE,N=N,m=m,M=M,
#'          alpha=alpha,H=H,sigma=sigma,freq='L',
#'          disable_X=FALSE,levy_increments=NULL)
#'
#' Retrieve_stats(paths=Y,true_val=sigma,Est=sigma_hat,t1=t1,k=2,alpha=alpha,H=H,freq="L")
#' @export
Retrieve_stats<-function(paths,true_val,Est,...){

    r<-apply(paths,MARGIN=2,Est,...)

    v<-var(r) # n-1 variance
    bias<-function(r,true_val){mean(r)-true_val}
    b<-bias(r,true_val)
    list(estimates=r,variance=v,bias=b)
}
