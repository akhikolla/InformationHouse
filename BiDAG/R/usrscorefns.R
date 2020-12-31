### This function returns the objects needed to evaluate the user defined score

# edgepf is the factor to penalise each edge
# edgepmat is a matrix of edge penalisation factors
# chi is the general number of pseudo count
# delta is a scaling factor for zero pseudo counts when the parents are not all on
# eta is a scaling factor for zero pseudo counts when the node has no parents 
usrscoreparameters <- function(initparam, usrpar=list(pctesttype="usrCItest", edgepf=2, edgepmat=NULL, 
                                                      chi=0.5, delta=NULL, eta=NULL)){
  
  if(is.null(usrpar$chi)) {usrpar$chi <- 0.5}
  if(is.null(usrpar$edgepf)) {usrpar$edgepf <- 2}
  
  initparam$chi <- usrpar$chi #1
  initparam$pf <- usrpar$edgepf
  
  if(is.null(usrpar$delta)) {usrpar$delta <- 100*initparam$chi}
  if(is.null(usrpar$eta)) {usrpar$eta <- 10*initparam$chi}
  
  initparam$delta <- usrpar$delta
  initparam$eta <- usrpar$eta
  
  if (is.null(initparam$weightvector)) {
    initparam$N <- nrow(initparam$data)
    initparam$d1 <- initparam$data
    initparam$d0 <- (1-initparam$data)
  } else {
    initparam$N <- sum(initparam$weightvector)
    initparam$d1 <- initparam$data*initparam$weightvector
    initparam$d0 <- (1-initparam$data)*initparam$weightvector
  }
  
  maxparents <- ncol(initparam$data) - 1
  initparam$scoreconstvec <- rep(0, maxparents+1)
  
  if (is.null(usrpar$edgepmat)) {
    initparam$logedgepmat <- NULL
  } else {
    initparam$logedgepmat <- log(usrpar$edgepmat)
  }
  
  initparam$scoreconstvec <- lgamma(initparam$chi/2)+lgamma((1+initparam$delta)*initparam$chi/4)-3*lgamma(initparam$chi/4)-lgamma(initparam$delta*initparam$chi/4) - c(0:maxparents)*log(initparam$pf)
  initparam$scoreconstvec[1] <- lgamma((1+initparam$eta)*initparam$chi/2)-lgamma(initparam$chi/2)-lgamma(initparam$eta*initparam$chi/2) # simpler result with no parents
  
  initparam
  
}

### This function evaluates the log score of a node given its parents
usrDAGcorescore <- function (j,parentnodes,n,param) {
  
  lp<-length(parentnodes) # number of parents
  chi<-param$chi
  scoreconstvec<-param$scoreconstvec
  
  switch(as.character(lp),
         "0"={# no parents
           N1<-sum(param$d1[,j])
           N0<-sum(param$d0[,j])
           NT<-N0+N1
           corescore <- scoreconstvec[lp+1] + lgamma(N0+param$eta*chi/2) + lgamma(N1+chi/2) - lgamma(NT+(1+param$eta)*chi/2)
         },
         "1"={# one parent
           corescore<-scoreconstvec[lp+1]  
           summys<-param$data[,parentnodes]
           for(i in 0:1){
             totest<-which(summys==i)
             N1<-sum(param$d1[totest,j])
             N0<-sum(param$d0[totest,j])
             NT<-N0+N1
             if(i==0){
               corescore <- corescore + lgamma(N0+param$delta*chi/4) + lgamma(N1+chi/4) - lgamma(NT+(1+param$delta)*chi/4)
             } else {
               corescore <- corescore + lgamma(N0+chi/4) + lgamma(N1+chi/4) - lgamma(NT+chi/2)
             }
           }
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - param$logedgepmat[parentnodes, j]
           }
         },     
         { # more parents
           summys<-1*(rowSums(param$data[,parentnodes])==lp)
           N1s<-collectC(summys,param$d1[,j],2)
           N0s<-collectC(summys,param$d0[,j],2)
           NTs<-N1s+N0s       
           corescore <- scoreconstvec[lp+1] + sum(lgamma(N0s+c(param$delta,1)*chi/4)) + 
             sum(lgamma(N1s+chi/4)) - sum(lgamma(NTs+c(1+param$delta,2)*chi/4))
           if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
             corescore <- corescore - sum(param$logedgepmat[parentnodes, j])
           }
         })
  
  corescore 
  
}



### This function prints out the user defined score objects

usrprint.scoreparameters <- function (x,...){
  
  cat("CBN mimic score is being used","\n")
  cat("Prior pseudo counts:", x$chi,"\n")
  cat("Edge penalization factor:", x$pf,"\n")
  cat("Prior pseudo count scaling with 1 parent:", x$delta,"\n")
  cat("Prior pseudo count scaling with no parent:", x$eta,"\n")
  if(is.null(x$weightvector)) {
    cat("Data is not weighted\n")
  } else {
    "Data is weighted according to the weight vector"
  }
  cat("Score constant vector:", x$scoreconstvec,"\n")
}


### This function defines the CI tests for the starting skeleton

# actually not needed here as we revert to the bde score
# but just to test in case!

usrdefinestartspace <- function(alpha,param,cpdag,n){
  if(cpdag){
    pc.skel<-pc(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                verbose = FALSE)
    
  } else {
    pc.skel<-pcalg::skeleton(suffStat = list(d1=param$d1,d0=param$d0,data=param$data),
                             indepTest = weightedbinCItest, alpha = alpha, labels = colnames(param$data),
                             verbose = FALSE)
  }
  pc.skel
}
