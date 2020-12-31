
#### MCestimLFSM function ####
#' Numerical properties of statistical estimators operating on the linear fractional stable motion.
#'
#' The function is useful, for instance, when one needs to compute standard deviation of \eqn{\widehat \alpha_{high}}
#' estimator given a fixed set of parameters.
#'
#' MCestimLFSM  performs Monte-Carlo experiments to compute parameters according to procedure Inference.
#' More specifically, for each element of s it generates Nmc lfsm sample paths with length equal to s[i], performs the statistical
#' inference on each, obtaining the estimates, and then returns their different statistics. It is vital that the estimator
#' returns a list of named parameters (one or several of 'sigma', 'alpha' and 'H'). MCestimLFSM uses the names to lookup the true
#' parameter value and compute its bias.
#'
#' For sample path generation MCestimLFSM uses a light-weight version of path, path_fast. In order to be applied,
#' function Inference must accept argument 'path' as a sample path.
#' @param Nmc Number of Monte Carlo repetitions
#' @param Inference statistical function to apply to sample paths
#' @param s sequence of path lengths
#' @inheritParams path
#' @param fr frequency. Either "H" or "L"
#' @param ... parameters to pass to Inference
#' @return It returns a list containing the following components:
#'     \item{data}{a data frame, values of the estimates depending on path length s}
#'     \item{data_nor}{a data frame, normalized values of the estimates depending on path length s}
#'     \item{means, biases, sds}{data frames: means, biases and standard deviations of the estimators depending on s}
#'     \item{Inference}{a function used to obtain estimates}
#'     \item{alpha, H, sigma}{the parameters for which MCestimLFSM performs path generation}
#'     \item{freq}{frequency, either 'L'  for low- or 'H' for high frequency}
#' @export
#' @examples
#' #### Set of global parameters ####
#' m<-25; M<-60
#' p<-.4; p_prime<-.2; k<-2
#' t1<-1; t2<-2
#' NmonteC<-5e1
#' S<-c(1e2,3e2)
#' alpha<-1.8; H<-0.8; sigma<-0.3
#'
#'
#' # How to plot empirical density
#' \donttest{
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,ContinEstim,
#'                      t1=t1,t2=t2,p=p,k=k)
#' l_plot<-Plot_dens(par_vec=c('sigma','alpha','H'),
#'                   MC_data=theor_3_1_H_clt, Nnorm=1e7)
#'
#' }
#'
#' # For MCestimLFSM() it is vital that the estimator returns a list of named parameters
#'
#' H_hat_f <- function(p,k,path) {hh<-H_hat(p,k,path); list(H=hh)}
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,H_hat_f,
#'                      p=p,k=k)
#'
#'
#' # The estimator can return one, two or three of the parameters.
#'
#' est_1 <- function(path) list(H=1)
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,est_1)
#'
#' est_2 <- function(path) list(H=0.8, alpha=1.5)
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,est_2)
#'
#' est_3 <- function(path) list(sigma=5, H=0.8, alpha=1.5)
#' theor_3_1_H_clt<-MCestimLFSM(s=S,fr='H',Nmc=NmonteC,
#'                      m=m,M=M,alpha=alpha,H=H,
#'                      sigma=sigma,est_3)
MCestimLFSM <- function (Nmc, s, m, M, alpha, H, sigma, fr, Inference, ...){

    i <- integer(0)
    ind <- integer(0)
    means <- data.frame()
    sds <- data.frame()
    data_whole <- data.frame()
    data_nor <- data.frame()
    s_whole <- integer(0) # factor vector with s for data_whole
    s_means<-vector(mode='numeric')
    s_sds<-vector(mode='numeric')

    for (i in 1:length(s)) {

        data <- matrix()
        data_nor_i <- data.frame()

        # path simulation + statistics extraction
        data <- foreach(ind = 1:Nmc, .combine = rbind, .packages = "stabledist", .export = LofF, .inorder = FALSE) %dopar% {

            sample_path <- path_fast(N = s[i], m = m, M = M, alpha = alpha, H = H, sigma = sigma, freq = fr)

            # Inference must return a number, or NA or NaN for every parameter, or throw an error
            # c(1,2,'abcd') will result in halting of execution
            if (is.null(unlist(formals(Inference)["freq"]))) {
                LL <- tryCatch(Inference(path = sample_path, ...),
                               error=function(c) 'Inference produced error(s)')
            }

            else {
                LL <- tryCatch(Inference(path = sample_path, freq = fr, ...),
                               error=function(c) 'Inference produced error(s)')
            }

            if (!is.character(LL)) # find a character in one element at least!
            as.data.frame(LL)
        }

        if (!is.null(data) & !is.null(names(data))) {

            data <- as.matrix(data) #!!!! this double-conversion seems to work best
            vect_s<-rep(s[i], length.out=nrow(data))
            means_i <- aaply(data, .margins=2, .fun=mean, na.rm = TRUE)
            sds_i   <- aaply(data, .margins=2, .fun=sd, na.rm = TRUE)

            data_nor_i <- sweep(data, MARGIN=2, STATS=means_i, FUN = "-", check.margin = TRUE)
            data_nor_i <- sweep(data_nor_i, MARGIN=2, STATS=sds_i, FUN = "/", check.margin = TRUE)

            # what we will aggreggate

            if(!is.null(means_i)) {means <- rbind(means, means_i); s_means<-c(s_means,s[i])}
            if(!is.null(sds_i)) {sds <- rbind(sds, sds_i); s_sds<-c(s_sds,s[i])}
            data_nor<- rbind(data_nor, data_nor_i)

            data_whole <- rbind(data_whole,data)
            s_whole <- c(s_whole,vect_s)
        }

        else {
            # there used to be stop function instead of warning
            warning("The inference function either hasn't produced any data, or the estimates have no names")
        }
    }

    nms <- names(data_whole)

    ########

    # here we create a data frame with the real parameters
    # to be able to subtract them from data in the future
    real_pars<-data_whole[1,]
    for(ind_names in seq_along(nms)) real_pars[ind_names]<-get(nms[ind_names])

    # whole data + lengths
    data_whole$s <- s_whole

    # data_nor is essentually normalized whole_data
    data_nor$s <- s_whole
    names(means) <- nms; names(sds) <- nms # names

    # biases
    biases <- sweep(as.matrix(means), MARGIN=2, STATS=as.matrix(real_pars), FUN = "-", check.margin = TRUE)
    means$s <- s_means;  sds$s <- s_sds
    biases<-as.data.frame(biases); biases$s <- s_means

    list(data=data_whole, data_nor=data_nor, means=means, sds=sds, biases=biases, Inference=Inference, params=list(sigma=sigma, alpha=alpha, H=H), freq=fr)
}
