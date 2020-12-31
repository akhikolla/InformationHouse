##' MuChPoint fitting procedure
##'
##' Produce a block-wise estimation of a symmetric matrix.
##'
##' @param Y symmetric matrix of observations.
##' @param Lmax a positive integer less than number of columns (and number of rows).
##' By default, \code{nrow(Y)/2}.
##' @param N a positive integer vector less than number of columns (and number of rows).
##' N is used when the break-points are known.
##' By default, \code{NULL}.
##' @param cores a positive integer giving the number of cores used. If you use windows,
##'  the parallelization is impossible.
##' By default, 1.
##' @param verbose logical. To display the progression bars. By default TRUE.
##'
##' @references Article: BRAULT V., OUADAH S., SANSONNET L. and LEVY-LEDUC C. Nonparametric
##' homogeneity tests and multiple change-point estimation for analyzing large Hi-C data matrices.
##' Journal of Multivariate Analysis, 2017
##'
##' @rdname MuChPoint-proc
##'
##' @examples
##' require(MuChPoint)
##' mu=c(rep(c(rep(1,25),rep(0,25)),3))%*%t(rep(c(rep(0,25),rep(1,25)),3))
##' Y=matrix(rnorm(150^2,0,5),150)+mu+t(mu)
##' Y=as.matrix(Matrix::forceSymmetric(Y))
##' res=MuChPoint(Y)
##' plot(res,Y,L=5,shiny=FALSE)
##' plot(res,Y,L=1:5,shiny=FALSE,ask=FALSE)
##'
##'
##' @export MuChPoint

MuChPoint=function(Y,Lmax=nrow(Y)/2,N=NULL,cores=1,verbose=TRUE){
  ### Verification et param de bases
  ## Verification Y
  if (!is.matrix(Y)){
    Y=as.matrix(Y)
  }
  if (!is.numeric(Y)){
    stop("Y must be a numerical matrix")
  }

  n=nrow(Y)# Taille matrice

  if (n!=ncol(Y)){
    stop("Y must have the same number of lines and columns")
  }

  ## Verfication cores
  if (floor(cores)!=cores){
    warning("\n cores must be an positive integer, enforcing cores to '1'")
    cores <- 1
  }
  if (cores!=1){
    if (Sys.info()[['sysname']] == "Windows") {
      warning("\n Windows does not support fork, enforcing cores to '1'.")
      cores <- 1
    }else{
      if (cores>parallel::detectCores()){
        warning(paste("\n Available cores is ",as.character(parallel::detectCores()),
                      ", enforcing cores to '",as.character(parallel::detectCores()),"'",
                      sep=""))
        cores=parallel::detectCores()
      }
      if (cores<1){
        warning("\n cores must be an positive integer, enforcing cores to '1'")
        cores <- 1
      }
      if (cores>n){
        cores=n
      }
    }
  }

  ### Verification Verbose
  if (!is.logical(verbose)){
    warning("\n verbose must be logical, enforcing verbose to 'FALSE'")
    verbose=FALSE
  }

  ### Verification Lmax
  if (is.null(Lmax)){
    if(is.null(N)){
      stop("\n N or Lmax is required")
    }
  }else{
    if (!is.numeric(Lmax)){
      if(is.null(N)){
        stop("\n Lmax must be an integer between 1 and n")
      }else{
        warning("\n Lmax is not used because it must be an integer")
        Lmax=NULL
      }
    }else{
      if ((floor(Lmax)!=Lmax)||(Lmax<1)||(Lmax>=n)){
        if(is.null(N)){
          stop("\n Lmax must be an integer between 1 and n")
        }else{
          warning("\n Lmax is not used because it must be an integer between 1 and n")
          Lmax=NULL
        }
      }
      if ((!is.null(Lmax))&&(!is.null(N))){
        if (Lmax!=nrow(Y)/2){
          warning("\n N is not used because Lmax is given")
          N=NULL
        }
      }
    }
  }


  ### Calcul des rangs
  if (verbose){
    cat("Computation of Rij\n")
    pb=utils::txtProgressBar(min=0,max=n,style=3)
    CompuRij=function(i){
      utils::setTxtProgressBar(pb, i)
      rank(Y[i,])
    }
    if (cores==1){
      Rij=matrix(unlist(lapply(1:n,CompuRij)),n,byrow=TRUE)
    }else{
      Rij=matrix(unlist(parallel::mclapply(1:n,CompuRij,mc.cores=cores)),n,byrow=TRUE)
    }
  }else{
    if (cores==1){
      Rij=matrix(unlist(lapply(1:n,function(i){rank(Y[i,])})),n,byrow=TRUE)
    }else{
      Rij=matrix(unlist(parallel::mclapply(1:n,function(i){rank(Y[i,])},mc.cores=cores)),n,byrow=TRUE)
    }
  }

  ### Calcul si N donne
  if (!is.null(N)){

    ## Verification N
    if (any(floor(N)!=N)){
      stop("N must be an increasing positive integer vector less than n")
    }

    L=length(N) # Nb ruptures

    if (L>1){
      if (any((N[2:L]-N[1:(L-1)])<=0)){
        N=sort(N)
        warning("\n N has been sorted by increasing order")
      }
    }
    if ((N[1]<1)||(N[L]>=n)){
      stop("N must be an increasing positive integer vector between 1 and n-1")
    }

    selection=FALSE

    Nbis=c(0,N,n) ## ajout des valeurs 0 et n
    if (verbose){
      cat("\nComputation of Sn\n")
      pb=utils::txtProgressBar(min=0,max=L+1,style=3)
      if (cores==1){
        ### Calcul de la statistique
        S=0
        for (l in (1:(L+1))){ ## Decalage de 0:L
          S=S+(Nbis[l+1]-Nbis[l])*
            sum(((apply(Rij[,(Nbis[l]+1):(Nbis[l+1])],1,sum))/(Nbis[l+1]-Nbis[l]) ##Rbar
                 -(n+1)/2)^2)
          utils::setTxtProgressBar(pb,l)
        }
      }else{
        ### Calcul de la statistique
        if (cores>(L+1)){
          cores=L+1
        }
        S=Reduce("+",parallel::mclapply(1:(L+1),function(l){(Nbis[l+1]-Nbis[l])*
            sum(((apply(Rij[,(Nbis[l]+1):(Nbis[l+1])],1,sum))/(Nbis[l+1]-Nbis[l]) ##Rbar
                 -(n+1)/2)^2)},mc.cores=cores))
        utils::setTxtProgressBar(pb,l)
      }
    }else{
      if (cores==1){
        ### Calcul de la statistique
        S=0
        for (l in (1:(L+1))){ ## Decalage de 0:L
          S=S+(Nbis[l+1]-Nbis[l])*
            sum(((apply(Rij[,(Nbis[l]+1):(Nbis[l+1])],1,sum))/(Nbis[l+1]-Nbis[l]) ##Rbar
                 -(n+1)/2)^2)
        }
      }else{
        ### Calcul de la statistique
        if (cores>(L+1)){
          cores=L+1
        }
        S=Reduce("+",parallel::mclapply(1:(L+1),function(l){(Nbis[l+1]-Nbis[l])*
            sum(((apply(Rij[,(Nbis[l]+1):(Nbis[l+1])],1,sum))/(Nbis[l+1]-Nbis[l]) ##Rbar
                 -(n+1)/2)^2)},mc.cores=cores))
      }
    }
    S=S*4/n^2
    bt=matrix(0,0,0)
  }else{
    ## Recherche des ruptures si N non donne
    ### Calcul de C
    # Cn1n2=matrix(0,n,n)
    # if (verbose){
    #   cat("\nComputation of Cn[n1,n2]\n")
    #   pbbis=utils::txtProgressBar(min=0,max=n*(n+1)/2,style=3) ### Decalage pour n1
    #   C=function(n1,n2){
    #     if (n2>(n1+1)){
    #       utils::setTxtProgressBar(pbbis,n1*(2*n-n1+1)/2+n2)
    #       (n2-n1)*sum(((apply(Rij[,(n1+1):n2],1,sum))/(n2-n1) ##Rbar
    #                    -(n+1)/2)^2)
    #     }else{
    #       0
    #     }
    #   }
    # }else{
    #   C=function(n1,n2){
    #     if (n2>(n1+1)){
    #       (n2-n1)*sum(((apply(Rij[,(n1+1):n2],1,sum))/(n2-n1) ##Rbar
    #                    -(n+1)/2)^2)
    #     }else{
    #       0
    #     }
    #   }
    # }
    # if (cores==1){
    #   Cn1n2[1:(n-1),2:n]=matrix(unlist(lapply(0:(n-2),function(n1){
    #     unlist(lapply(2:n,function(n2){C(n1,n2)}))})),nrow=n-1,byrow=TRUE)
    # }else{
    #   Cn1n2[1:(n-1),2:n]=matrix(unlist(parallel::mclapply(0:(n-2),function(n1){
    #     unlist(lapply(2:n,function(n2){C(n1,n2)}))},mc.cores=cores)),
    #     nrow=n-1,byrow=TRUE)
    # }
    # diag(Cn1n2)=colSums(((Rij[,1:n])/-(n+1)/2)^2) ## Cas ou n1+1=n2
    ### Version Cpp
    if (verbose){
      cat("\nComputation of Cn[n1,n2]\n")
    }
    Cn1n2=Compute_Cn1n2(Rij)
    ### Calcul de I_n0(L) (optimise/A paralleliser)
    IL=matrix(0,Lmax+1,n)
    ind=IL
    IL[1,]=Cn1n2[1,]
    ind[1,]=0
    if (verbose){
      cat("\nComputation of Sn and N\n")
      pb=utils::txtProgressBar(min=0,max=2*(Lmax+1),style=3)
      if (Lmax>0){
        for (L in 2:(Lmax+1)){
          utils::setTxtProgressBar(pb,L)
          ind[L,L:n]=sapply(L:n,function(n0){which.max(IL[L-1,L:(n0-1)]+Cn1n2[(L+1):n0,n0])})+L-1
          IL[L,L:n]=sapply(L:n,function(n0){IL[L-1,ind[L,n0]]+Cn1n2[ind[L,n0]+1,n0]})
        }
      }
      bt=matrix(0,Lmax+1,Lmax+1)
      for (L in 1:(Lmax+1)){
        utils::setTxtProgressBar(pb,Lmax+1+L)
        bt[L,L]=ind[L,n]
        if (L>1){
          for (k in (L-1):1){
            bt[L,k]=ind[k,bt[L,k+1]]
          }
        }
      }
    }else{
      if (Lmax>0){
        for (L in 2:(Lmax+1)){
          ind[L,L:n]=sapply(L:n,function(n0){which.max(IL[L-1,L:(n0-1)]+Cn1n2[(L+1):n0,n0])})+L-1
          IL[L,L:n]=sapply(L:n,function(n0){IL[L-1,ind[L,n0]]+Cn1n2[ind[L,n0]+1,n0]})
        }
      }
      bt=matrix(0,Lmax,Lmax)
      bt[1,1]=ind[1,n]
      if (Lmax>1){
        for (L in 2:Lmax){
          bt[L,L]=ind[L,n]
          for (k in (L-1):1){
            bt[L,k]=ind[k,bt[L,k+1]]
          }
        }
      }
    }
    S=IL[-1,n]*4/n^2
    bt=bt[-1,-1]
    N=bt[Lmax,]
  }
  methods::new(Class="MuChPoint", S=S,N=N,bt=bt)
}


##' @useDynLib MuChPoint
##' @importFrom Rcpp sourceCpp
NULL
