#' Path array generator
#'
#' The function takes a list of parameters (alpha, H) and uses \code{\link{expand.grid}} to obtain all possible combinations of them.
#' Based on each combination, the function simulates an lfsm sample path. It is meant to be used in conjunction with function \code{\link{Plot_list_paths}}.
#' @param l a list of parameters to expand
#' @inheritParams path
#' @return The returned value is a data frame containing paths and the corresponding values of alpha, H and frequency.
#' @examples
#' l=list(H=c(0.2,0.8),alpha=c(1,1.8), freq="H")
#' arr<-Path_array(N=300,m=30,M=100,l=l,sigma=0.3)
#' str(arr)
#' head(arr)
#' @export
Path_array<-function(N,m,M,l,sigma){
  # rename as path_per_par?
  table<-as.matrix(expand.grid(l))

  # this function takes elements of a list and applies rlfsm to it
  gen_path_from_el<-function(el,N,m,M,sigma) {

    a<-as.numeric(el[["alpha"]])
    h<-as.numeric(el[["H"]])
    f<-el[["freq"]]
    P<-path_fast(N,m,M,alpha=a,H=h,freq=f,sigma)#disable_X=FALSE,...)
    ProcP<-cbind(n=(1:(N+1)),X=P,alpha=a,H=h,freq=f) # N+1 because we add 0 in the beginning of lfsm and Levy motion
    ProcP
  }

  res<-adply(table,1,gen_path_from_el,N,m,M,sigma)

  # Here we coerce those factors that should be numeric
  # X is not coerced correctly
  res$n<-as.numeric(as.character(res$n)) ; res$X<-as.numeric(as.character(res$X))
  res$X1 <- NULL
  res
}




#### Image rendering functions ####
#' Rendering of path lattice
#'
#' @param arr a data frame produced by \code{\link{Path_array}}.
#' @examples
#' l=list(H=c(0.2,0.8),alpha=c(1,1.8), freq="H")
#' arr<-Path_array(N=300,m=30,M=100,l=l,sigma=0.3)
#' p<-Plot_list_paths(arr)
#' p
#' @export
Plot_list_paths<-function(arr){

  # avoiding NOTEs during CRAN checks
  n<-vector(mode = "numeric", length = 0)
  X<-vector(mode = "numeric", length = 0)
  beta_ind<-vector(mode = "logical", length = 0)

  # Draws either a transparent line over jumps or a solid line  in cont. case.
  arr$beta_ind<-ifelse(as.numeric(as.character(arr$H))-
                     1/as.numeric(as.character(arr$alpha))>0,1,0.25)


  pl <- ggplot(arr, aes(x=n, y=X)) + geom_point(size = 0.25, colour = "brown")
  pl<-pl + geom_line(aes(alpha=beta_ind), colour ="brown", size = 0.8)
  pl<-pl + facet_wrap(H ~ alpha,scales = "free", labeller = label_both)
  pl<-pl + theme(strip.background = element_blank(), strip.placement = "outside")
  pl<-pl + theme(legend.position = "none")
  pl
}

