DAGbinarytablescore<-function(j,parentnodes,n,param,parenttable,tablemaps,numparents,numberofparentsvec){
  lp<-length(parentnodes) # number of parents
  noparams<-2^lp # number of binary states of the parents
  
  corescores<-rep(NA,noparams)
  
  chi<-param$chi
  scoreconstvec<-param$scoreconstvec
  
  N1slist<-vector("list", noparams)
  N0slist<-vector("list", noparams)
  
  if(lp==0){ # no parents
    N1<-sum(param$d1[,j])
    N0<-sum(param$d0[,j])
    NT<-N0+N1
    corescores[noparams] <- scoreconstvec[lp+1] + lgamma(N0+chi/(2*noparams)) + lgamma(N1+chi/(2*noparams)) - lgamma(NT+chi/noparams)
  } else {
    
    if(lp==1){
      summys<-param$data[,parentnodes]
    } else {
      summys<-colSums(2^(c(0:(lp-1)))*t(param$data[,parentnodes]))
    }
    
    N1s<-collectC(summys,param$d1[,j],noparams)
    N0s<-collectC(summys,param$d0[,j],noparams)
    
    N1slist[[noparams]]<-N1s
    N0slist[[noparams]]<-N0s
    
    NTs<-N1s+N0s
    
    corescores[noparams] <- scoreconstvec[lp+1] + sum(lgamma(N0s+chi/(2*noparams))) + sum(lgamma(N1s+chi/(2*noparams))) - sum(lgamma(NTs+chi/noparams))
    
    if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
      corescores[noparams] <- corescores[noparams] - sum(param$logedgepmat[parentnodes, j])
    }
    
    for (jj in (noparams-1):1){ # use poset to combine sets
      
      lplocal<-numberofparentsvec[jj] # size of parent set
      noparamslocal<-2^lplocal # number of parameters of this size, for scoring later
      
      missingparentindex<-parenttable[tablemaps$backwards[noparams-tablemaps$forward[jj]+1],1] # get first element of complement of parent set
      higherlayer<-tablemaps$backwards[tablemaps$forward[jj]+2^(missingparentindex-1)] # get row of poset element with this missing parent included
      missingparent<-which(parenttable[higherlayer,]==missingparentindex) # which component it is in the higher layer
      
      
      N1stemp<-N1slist[[higherlayer]] # map to the previous lists with missing parent added
      N0stemp<-N0slist[[higherlayer]] 
      
      # we take the N1s and N0s calculated previously for the case with this element included
      # since we know which power of 2 it takes in the mapping to summys, we can marginalise it out
      size1<-2^(missingparent-1)
      size2<-2^(lplocal-missingparent+1)
      
      elementstocombine<-as.vector(t(t(matrix(c(1:size1),nrow=size1,ncol=size2))+2*(c(1:size2)-1)*size1)) 
      # collect the elements we want to combine to remove the missing parent from the previous tables
      
      N1s<-N1stemp[elementstocombine]+N1stemp[elementstocombine+size1]
      N0s<-N0stemp[elementstocombine]+N0stemp[elementstocombine+size1]
      
      
      N1slist[[jj]]<-N1s
      N0slist[[jj]]<-N0s
      
      NTs<-N1s+N0s
      
      corescores[jj] <- scoreconstvec[lplocal+1] + sum(lgamma(N0s+chi/(2*noparamslocal))) + sum(lgamma(N1s+chi/(2*noparamslocal))) - sum(lgamma(NTs+chi/noparamslocal))
      
      if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
        if(lplocal>0) {
          localparents <- parentnodes[parenttable[jj, 1:lplocal]]
          if(length(localparents)>0) {
            corescores[jj] <- corescores[jj] - sum(param$logedgepmat[localparents, j])
          }
        }
      }
    }
    
  }
  
  return(corescores)
}


DAGbinarytablescoreplus1<-function(j,parentnodes,additionalparent,n,param,parenttable,tablemaps,numparents,numberofparentsvec){
  
      lp<-length(parentnodes) # number of parents
      noparams<-2^lp # number of binary states of the parents
      
      allparents<-c(parentnodes,additionalparent) # combine the sets, but put the additional one last!
      
      lpadd<-lp+1 # including the additional parent
      noparamsadd<-2*noparams 
      chi<-param$chi
      scoreconstvec<-param$scoreconstvec
      
      corescores<-rep(NA,noparams)
      
      N1slist<-vector("list", noparams)
      N0slist<-vector("list", noparams)
      
      if(lpadd==1){
        summys<-param$data[,allparents]
      } else {
        summys<-colSums(2^(c(0:(lpadd-1)))*t(param$data[,allparents]))
      }
      
      N1s<-collectC(summys,param$d1[,j],noparamsadd)
      N0s<-collectC(summys,param$d0[,j],noparamsadd)
      
      N1slist[[noparams]]<-N1s
      N0slist[[noparams]]<-N0s
      
      NTs<-N1s+N0s
      
      corescores[noparams] <- scoreconstvec[lpadd+1] + sum(lgamma(N0s+chi/(2*noparamsadd))) + sum(lgamma(N1s+chi/(2*noparamsadd))) - sum(lgamma(NTs+chi/noparamsadd))
      
      if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
        corescores[noparams] <- corescores[noparams] - sum(param$logedgepmat[allparents, j])
      }
      
      if(lpadd>1){ # otherwise there are no further terms to compute!  
        
        for (jj in (noparams-1):1){ # use poset to combine sets
          
          lplocal<-numberofparentsvec[jj]+1 # size of parent set
          noparamslocal<-2^lplocal # number of parameters of this size, for scoring later
          
          missingparentindex<-parenttable[tablemaps$backwards[noparams-tablemaps$forward[jj]+1],1] # get first element of complement of parent set
          higherlayer<-tablemaps$backwards[tablemaps$forward[jj]+2^(missingparentindex-1)] # get row of poset element with this missing parent included
          missingparent<-which(parenttable[higherlayer,]==missingparentindex) # which component it is in the higher layer
          
          N1stemp<-N1slist[[higherlayer]] # map to the previous lists with missing parent added
          N0stemp<-N0slist[[higherlayer]] 
          
          # we take the N1s and N0s calculated previously for the case with this element included
          # since we know which power of 2 it takes in the mapping to summys, we can marginalise it out
          size1<-2^(missingparent-1)
          size2<-2^(lplocal-missingparent+1)
          
          elementstocombine<-as.vector(t(t(matrix(c(1:size1),nrow=size1,ncol=size2))+2*(c(1:size2)-1)*size1)) # collect the elements we want to combine to remove the missing parent from the previous tables
          
          N1s<-N1stemp[elementstocombine]+N1stemp[elementstocombine+size1]
          N0s<-N0stemp[elementstocombine]+N0stemp[elementstocombine+size1]
          
          N1slist[[jj]]<-N1s
          N0slist[[jj]]<-N0s
          
          NTs<-N1s+N0s
          
          #lplocal+1 because we have 1 additional parent and indexing in scoreconstvec started with 0
          corescores[jj] <- scoreconstvec[lplocal+1] + sum(lgamma(N0s+chi/(2*noparamslocal))) + sum(lgamma(N1s+chi/(2*noparamslocal))) - sum(lgamma(NTs+chi/noparamslocal))
          
          if (!is.null(param$logedgepmat)) { # if there is an additional edge penalisation
            if (lplocal>1) {
              localparents <- c(parentnodes[parenttable[jj, 1:(lplocal-1)]], additionalparent)
            } else {
              localparents<-additionalparent
            }
            corescores[jj] <- corescores[jj] - sum(param$logedgepmat[localparents, j])
          }
        }
      
      } 
      return(corescores)
}
