DAGcattablescore<-function(j,parentnodes,n,param,parenttable,tablemaps,numparents,numberofparentsvec){
  
  lp <- length(parentnodes) # number of parents
  noparams <- 2^lp # number of binary states of the parents
  
  corescores <- rep(NA,noparams)
  
  chi <- param$chi
  scoreconstvec <- param$scoreconstvec
  
  Cj <- param$Cvec[j] # number of levels of j
  
  Nslist <- vector("list", noparams)
  
  if(lp==0){ # no parents
    Cp <- 1 # effectively 1 parent level
    summys <- rep(0, nrow(param$data))
    if(!is.null(param$weightvector)){
      Ns <- collectCcatwt(summys, param$data[,j], param$weightvector, Cp, Cj)
    } else{
      Ns <- collectCcat(summys, param$data[,j], Cp, Cj)
    }
    NTs <- rowSums(Ns)
    corescores[noparams] <- scoreconstvec[lp+1] + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
  } else {
    
    if(lp==1){
      Cp <- param$Cvec[parentnodes] # number of parent levels
      summys <- param$data[,parentnodes]
    } else {
      Cp <- prod(param$Cvec[parentnodes])
      # use mixed radix mapping to unique parent states
      summys <- colSums(cumprod(c(1,param$Cvec[parentnodes[-lp]]))*t(param$data[,parentnodes]))
    }
    
    if(!is.null(param$weightvector)){
      Ns <- collectCcatwt(summys, param$data[,j], param$weightvector, Cp, Cj)
    } else{
      Ns <- collectCcat(summys, param$data[,j], Cp, Cj)
    }
    
    Nslist[[noparams]] <- Ns
    
    NTs <- rowSums(Ns)
    
    corescores[noparams] <- scoreconstvec[lp+1] + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
    
    if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
      corescores[noparams] <- corescores[noparams] - sum(param$logedgepmat[parentnodes, j])
    }
    
    for (jj in (noparams-1):1){ # use poset to combine sets
      
      lplocal <- numberofparentsvec[jj] # size of parent set
      
      if(lplocal>0) {
        localparents <- parentnodes[parenttable[jj, 1:lplocal]]
        Cp <- prod(param$Cvec[localparents])
      } else {
        Cp <- 1
      }
      
      missingparentindex <- parenttable[tablemaps$backwards[noparams-tablemaps$forward[jj]+1],1] # get first element of complement of parent set
      higherlayer <- tablemaps$backwards[tablemaps$forward[jj]+2^(missingparentindex-1)] # get row of poset element with this missing parent included
      missingparent <- which(parenttable[higherlayer,]==missingparentindex) # which component it is in the higher layer
      
      Nstemp <- Nslist[[higherlayer]] # map to the previous lists with missing parent added
      
      # we take the Ns calculated previously for the case with this element included
      # since we know which power it takes in the mapping to summys, we can marginalise it out
      if(missingparent>1){
        size1 <- prod(param$Cvec[parentnodes[parenttable[higherlayer,1:(missingparent-1)]]])
      } else {
        size1 <- 1
      }
      if(missingparent<(lplocal+1)){
        size2 <- prod(param$Cvec[parentnodes[parenttable[higherlayer,missingparent:lplocal+1]]])
      } else {
        size2 <- 1
      }
      
      Cmissing <- param$Cvec[parentnodes[parenttable[higherlayer,missingparent]]]
      
      elementstocombine <- as.vector(t(t(matrix(c(1:size1),nrow=size1,ncol=size2))+Cmissing*(c(1:size2)-1)*size1)) 
      # collect the elements we want to combine to remove the missing parent from the previous tables
      
      Ns <- Nstemp[elementstocombine,]+Nstemp[elementstocombine+size1,]
      if(Cmissing>2){
        for(ii in 3:Cmissing-1)
          Ns <- Ns + Nstemp[elementstocombine+ii*size1,]
      }
      
      Nslist[[jj]] <- Ns
      
      if(lplocal==0){ # the combining turns the matrix into a vector at the end
        NTs <- sum(Ns)
      } else {
        NTs <- rowSums(Ns)
      }
      
      corescores[jj] <- scoreconstvec[lplocal+1] + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
      
      if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
        if(lplocal>0) {
          corescores[jj] <- corescores[jj] - sum(param$logedgepmat[localparents, j])
        }
      }
    }
    
  }
  
  return(corescores)
}


DAGcattablescoreplus1<-function(j,parentnodes,additionalparent,n,param,parenttable,tablemaps,numparents,numberofparentsvec){
  
  lp <- length(parentnodes) # number of parents
  noparams <- 2^lp # number of binary states of the parents
  
  corescores <- rep(NA,noparams)
  
  allparents <- c(parentnodes,additionalparent) # combine the sets, but put the additional one last!
  
  lpadd <- lp+1 # including the additional parent
  
  chi <- param$chi
  scoreconstvec <- param$scoreconstvec
  
  Cj <- param$Cvec[j] # number of levels of j
  
  Nslist <- vector("list", noparams)
  
  if(lpadd==1){
    Cp <- param$Cvec[allparents] # number of parent levels
    summys <- param$data[,allparents]
  } else {
    Cp <- prod(param$Cvec[allparents])
    summys <- colSums(cumprod(c(1,param$Cvec[allparents[-lpadd]]))*t(param$data[,allparents]))
  }
  
  if(!is.null(param$weightvector)){
    Ns <- collectCcatwt(summys, param$data[,j], param$weightvector, Cp, Cj)
  } else{
    Ns <- collectCcat(summys, param$data[,j], Cp, Cj)
  }
  
  Nslist[[noparams]] <- Ns
  
  NTs <- rowSums(Ns)
  
  corescores[noparams] <- scoreconstvec[lpadd+1] + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
  
  if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
    corescores[noparams] <- corescores[noparams] - sum(param$logedgepmat[allparents, j])
  }
  
  if(lpadd>1){ # otherwise there are no further terms to compute!  
    
    for (jj in (noparams-1):1){ # use poset to combine sets
      
      lplocal<-numberofparentsvec[jj]+1 # size of parent set
      
      if (lplocal>1) {
        localparents <- c(parentnodes[parenttable[jj, 1:(lplocal-1)]], additionalparent)
        Cp <- prod(param$Cvec[localparents])
      } else {
        localparents <- additionalparent
        Cp <- param$Cvec[localparents]
      }
      
      missingparentindex<-parenttable[tablemaps$backwards[noparams-tablemaps$forward[jj]+1],1] # get first element of complement of parent set
      higherlayer<-tablemaps$backwards[tablemaps$forward[jj]+2^(missingparentindex-1)] # get row of poset element with this missing parent included
      missingparent<-which(parenttable[higherlayer,]==missingparentindex) # which component it is in the higher layer
      
      Nstemp <- Nslist[[higherlayer]] # map to the previous lists with missing parent added
      
      # we take the Ns calculated previously for the case with this element included
      # since we know which power it takes in the mapping to summys, we can marginalise it out
      if(missingparent>1){
        size1 <- prod(param$Cvec[parentnodes[parenttable[higherlayer,1:(missingparent-1)]]])
      } else {
        size1 <- 1
      }
      if(missingparent<lplocal){
        size2 <- prod(param$Cvec[c(parentnodes[parenttable[higherlayer,(missingparent+1):lplocal]],additionalparent)])
      } else {
        size2 <- param$Cvec[additionalparent]
      }
      
      Cmissing <- param$Cvec[parentnodes[parenttable[higherlayer,missingparent]]]
      
      elementstocombine <- as.vector(t(t(matrix(c(1:size1),nrow=size1,ncol=size2))+Cmissing*(c(1:size2)-1)*size1)) 
      # collect the elements we want to combine to remove the missing parent from the previous tables
      
      Ns <- Nstemp[elementstocombine,]+Nstemp[elementstocombine+size1,]
      if(Cmissing>2){
        for(ii in 3:Cmissing-1)
          Ns <- Ns + Nstemp[elementstocombine+ii*size1,]
      }
      
      Nslist[[jj]] <- Ns
      
      NTs <- rowSums(Ns)
      
      #lplocal+1 because we have 1 additional parent and indexing in scoreconstvec started with 0
      corescores[jj] <- scoreconstvec[lplocal+1] + sum(lgamma(Ns+chi/(Cp*Cj))) - sum(lgamma(NTs + chi/Cp)) + Cp*lgamma(chi/Cp) - (Cp*Cj)*lgamma(chi/(Cp*Cj))
      
      if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
        corescores[jj] <- corescores[jj] - sum(param$logedgepmat[localparents, j])
      }
    }
  }
  
  return(corescores)
}

