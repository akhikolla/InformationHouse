
dettwobytwo <- function(D) {
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}

# The determinant of a 3 by 3 matrix
detthreebythree <- function(D){
  D[1,1]*(D[2,2]*D[3,3]-D[2,3]*D[3,2])-D[1,2]*(D[2,1]*D[3,3]-D[2,3]*D[3,1])+D[1,3]*(D[2,1]*D[2,3]-D[2,2]*D[3,1])
}

# The log of the BGe/BDe score, but simplified as much as possible
# see arXiv:1402.6863 

dettwobytwo <- function(D) {
  D[1,1]*D[2,2]-D[1,2]*D[2,1]
}

# The determinant of a 3 by 3 matrix
detthreebythree <- function(D){
  D[1,1]*(D[2,2]*D[3,3]-D[2,3]*D[3,2])-D[1,2]*(D[2,1]*D[3,3]-D[2,3]*D[3,1])+D[1,3]*(D[2,1]*D[2,3]-D[2,2]*D[3,1])
}

# The log of the BGe/BDe score, but simplified as much as possible
# see arXiv:1402.6863 
DAGcorescore<-function(j,parentnodes,n,param) {
  if(param$DBN){
    if(param$stationary) {
      internalparents <- parentnodes[which(parentnodes<=param$nsmall)]
      corescore <- DAGcorescore(j,parentnodes,param$n+param$nsmall,param$otherslices)+
      DAGcorescore(j,internalparents,param$n,param$firstslice) 
    } else {
      corescore<-0
      for(i in 1:length(param$paramsets)) {
        corescore<-corescore+DAGcorescore(j,parentnodes,param$n+param$nsmall,
                                          param$paramsets[[i]])
      }
    }
  } else if(param$type=="bge"){
    TN<-param$TN
    awpN<-param$awpN
    scoreconstvec<-param$scoreconstvec
    
    lp<-length(parentnodes) #number of parents
    awpNd2<-(awpN-n+lp+1)/2
    A<-TN[j,j]
    switch(as.character(lp),
           "0"={# just a single term if no parents
             corescore <- scoreconstvec[lp+1] -awpNd2*log(A)
           },
           
           "1"={# no need for matrices
             D<-TN[parentnodes,parentnodes]
             logdetD<-log(D)
             B<-TN[j,parentnodes]
             logdetpart2<-log(A-B^2/D)
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           },
           
           "2"={# can do matrix determinant and inverse explicitly
             # but this is numerically unstable for large matrices!
             # so we use the same approach as for 3 parents
             D<-TN[parentnodes,parentnodes]
             detD<-dettwobytwo(D)
             logdetD<-log(detD)
             B<-TN[j,parentnodes]
             #logdetpart2<-log(A-(D[2,2]*B[1]^2+D[1,1]*B[2]^2-2*D[1,2]*B[1]*B[2])/detD) #also using symmetry of D
             logdetpart2<-log(dettwobytwo(D-(B)%*%t(B)/A))+log(A)-logdetD
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           },
           
           {# otherwise we use cholesky decomposition to perform both
             D<-as.matrix(TN[parentnodes,parentnodes])
             choltemp<-chol(D)
             logdetD<-2*log(prod(choltemp[(lp+1)*c(0:(lp-1))+1]))
             B<-TN[j,parentnodes]
             logdetpart2<-log(A-sum(backsolve(choltemp,B,transpose=TRUE)^2))
             corescore <- scoreconstvec[lp+1]-awpNd2*logdetpart2 - logdetD/2
           })
  } else if (param$type=="bde"){
    
    lp<-length(parentnodes) # number of parents
    noparams<-2^lp # number of binary states of the parents
    chi<-param$chi
    scoreconstvec<-param$scoreconstvec
    
    switch(as.character(lp),
           "0"={# no parents
             N1<-sum(param$d1[,j])
             N0<-sum(param$d0[,j])
             NT<-N0+N1
             corescore <- scoreconstvec[lp+1] + lgamma(N0+chi/(2*noparams)) + lgamma(N1+chi/(2*noparams)) - lgamma(NT+chi/noparams)
           },
           "1"={# one parent
             corescore<-scoreconstvec[lp+1]  
             summys<-param$data[,parentnodes]
             for(i in 1:noparams-1){
               totest<-which(summys==i)
               N1<-sum(param$d1[totest,j])
               N0<-sum(param$d0[totest,j])
               NT<-N0+N1
               corescore <- corescore + lgamma(N0+chi/(2*noparams)) + lgamma(N1+chi/(2*noparams)) - lgamma(NT+chi/noparams)
             }
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - param$logedgepmat[parentnodes, j]
             }
           },     
           { # more parents
             summys<-colSums(2^(c(0:(lp-1)))*t(param$data[,parentnodes]))
             
             N1s<-collectC(summys,param$d1[,j],noparams)
             N0s<-collectC(summys,param$d0[,j],noparams)
             
             NTs<-N1s+N0s       
             corescore <- scoreconstvec[lp+1] + sum(lgamma(N0s+chi/(2*noparams))) + 
               sum(lgamma(N1s+chi/(2*noparams))) - sum(lgamma(NTs+chi/noparams))
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
             }
           })
    
  } else if (param$type=="bdecat"){
    
    lp<-length(parentnodes) # number of parents
    chi<-param$chi
    
    corescore <- param$scoreconstvec[lp+1] # starting value
    
    Cj <- param$Cvec[j] # number of levels of j
    
    switch(as.character(lp),
           "0"={# no parents
             Cp <- 1 # effectively 1 parent level
             summys <- rep(0, nrow(param$data))
           },
           "1"={# one parent
             Cp <- param$Cvec[parentnodes] # number of parent levels
             summys <- param$data[,parentnodes]
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - param$logedgepmat[parentnodes, j]
             }
           },     
           { # more parents
             Cp <- prod(param$Cvec[parentnodes])
             # use mixed radix mapping to unique parent states
             summys<-colSums(cumprod(c(1,param$Cvec[parentnodes[-lp]]))*t(param$data[,parentnodes]))
             if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
               corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
             }
           })
    
    if(!is.null(param$weightvector)){
      Ns <- collectCcatwt(summys, param$data[,j], param$weightvector, Cp, Cj)
    } else{
      Ns <- collectCcat(summys, param$data[,j], Cp, Cj)
    }
    
    NTs <- rowSums(Ns)
    
    corescore <- corescore + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
    
  } else if(param$type=="usr"){
    corescore <- usrDAGcorescore(j,parentnodes,n,param)
  } 
  
  return(corescore)
}






