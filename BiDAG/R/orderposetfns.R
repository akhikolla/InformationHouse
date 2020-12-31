#builds a table with banned scores for sampling version
poset.scores<-function(posetparenttable,scoretable,numberofparentsvec,rowmaps,n,plus1lists=NULL,
                       numparents,updatenodes=c(1:n)){
  orderscore<-list(length=n) #first level of list: the only level - list of matrices - 3 dimensions
  revnumberofparentsvec<-lapply(numberofparentsvec,rev)
  if (is.null(plus1lists)){
    for (j in updatenodes){
      len<-numparents[j] #maximum number of parents for node j
      binomcoefs<-choose(len,c(0:len))
      nrows<-nrow(posetparenttable[[j]])
      P_local<-vector("numeric",length=nrows)
      P_local[nrows] <-scoretable[[j]][1,1]
      maxoverall<-max(scoretable[[j]][,1])
      P_local[1]<-log(sum(exp(scoretable[[j]][,1]-maxoverall)))+maxoverall
      cutoff<-1

      if(nrows>2){
        for(level in 1:(len-1)){
          cutoff<-cutoff+binomcoefs[level]

          for (i in (nrows-1):cutoff)  {
            #so we go through all rows where non-zero entries more than number of banned parents
            # find the parents in the poset graph
            posetparentnodes <- posetparenttable[[j]][i,c(1:revnumberofparentsvec[[j]][i])]
            maxparents<-max(P_local[posetparentnodes])
            parentsum<-log(sum(exp(P_local[posetparentnodes]-maxparents)))+maxparents-log(len-revnumberofparentsvec[[j]][i]-level+1)
            conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1],1]
            maxoverall<-max(parentsum,conjugatescore)
            P_local[i]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
          }
        }
      }
      orderscore[[j]]<-as.matrix(P_local)
    }

    return(orderscore)
  } else {
    for (j in updatenodes) {
      len<-numparents[j] #maximum number of parents for node j
      binomcoefs<-choose(len,c(0:len))
      ll<-length(plus1lists$parents[[j]])+1
      nrows<-nrow(posetparenttable[[j]])
      P_local <- matrix(nrow=nrows,ncol=ll)

      for (li in 1:ll){
        P_local[nrows,li] <-scoretable[[j]][[li]][1,1]
        maxoverall<-max(scoretable[[j]][[li]][,1])
        P_local[1,li]<-log(sum(exp(scoretable[[j]][[li]][,1]-maxoverall)))+maxoverall
        cutoff<-1
        if(nrows>2){
          for(level in 1:(len-1)){
            cutoff<-cutoff+binomcoefs[level]
            for (i in (nrows-1):cutoff)  {
              #so we go through all rows where non-zero entries more than number of banned parents
              # find the parents in the poset graph
              posetparentnodes <- posetparenttable[[j]][i,c(1:revnumberofparentsvec[[j]][i])]
              maxparents<-max(P_local[posetparentnodes,li])
              parentsum<-log(sum(exp(P_local[posetparentnodes,li]-maxparents)))+maxparents-log(len-revnumberofparentsvec[[j]][i]-level+1)
              conjugatescore<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1],1]
              maxoverall<-max(parentsum,conjugatescore)
              P_local[i,li]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
            }
          }
        }

      }
      orderscore[[j]]<-P_local
    }
    return(orderscore)
  }

}

#builds a table with banned scores for max version
posetscoremax<-function(posetparenttable,scoretable,numberofparentsvec,rowmaps,n,plus1lists=NULL,
                        updatenodes=c(1:n)) {
    listy<-list()
    revnumberofparentsvec<-lapply(numberofparentsvec,rev)
    if (is.null(plus1lists)) {
      maxmatrix<-list()
      maxrows<-list()
          for (j in updatenodes) {
              nrows<-nrow(posetparenttable[[j]])
              P_local <- numeric(nrows)
              maxrow<-numeric(nrows)
              P_local[nrows]<-scoretable[[j]][1,1]
              maxrow[nrows]<-1

              #last element, when all nodes are banned 1st row in a score table
              #for all other we have
  if(nrows>1) {
    for (i in (nrows-1):1)  { # find the parents in the poset graph
      posetparentnodes <- posetparenttable[[j]][i,c(1:revnumberofparentsvec[[j]][i])]
      # take the maximum of the parent scores, and the conjugate score
      prevmax<- max(P_local[posetparentnodes])
      candmax<-scoretable[[j]][rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1],1]
      if (prevmax>candmax) {
        P_local[i]<-prevmax
        maxrow[i]<-maxrow[posetparentnodes[which.max(P_local[posetparentnodes])]]
          } else {
            P_local[i]<-candmax
            maxrow[i]<-rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1]
          }
      }
  }

  maxmatrix[[j]]<-as.matrix(P_local)
  maxrows[[j]]<- maxrow

}
listy$maxmatrix<-maxmatrix
listy$maxrows<-maxrows

return(listy)
} else {

  maxmatrix<-list()
  maxrows<-list()
  for (j in updatenodes)
  {
    ll<-length(plus1lists$parents[[j]])+1
    nrows<-nrow(posetparenttable[[j]])
    P_local <- matrix(nrow=nrows,ncol=ll)
    maxrow<-matrix(nrow=nrows,ncol=ll)
    for (li in 1:ll)
    {
      P_local[nrows,li]<-scoretable[[j]][[li]][1,1] #in plus1 lists last row means that all nodes banned but plus1
      maxrow[nrows,li]<-1
      #last element, when all nodes are banned 1st row in a score table
      #for all oter we have
      if(nrows>1) {
        for (i in (nrows-1):1) {
          # find the parents in the poset graph
          posetparentnodes <- posetparenttable[[j]][i,c(1:revnumberofparentsvec[[j]][i])]
          # take the maximum of the parent scores, and the conjugate score
          prevmax<- max(P_local[posetparentnodes,li])
          candmax<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1],1]
          if (prevmax>candmax) {
            P_local[i,li]<-prevmax
            maxrow[i,li]<-maxrow[posetparentnodes[which.max(P_local[posetparentnodes,li])],li]
          } else {
            P_local[i,li]<-candmax
            maxrow[i,li]<-rowmaps[[j]]$backwards[nrows-rowmaps[[j]]$forward[i]+1]

          }
        }
      }
    }
    maxmatrix[[j]]<-P_local
    maxrows[[j]]<-maxrow

  }
  listy$maxmatrix<-maxmatrix
  listy$maxrow<-maxrows

  return(listy)

}

}



parentsmapping<-function(parenttable,numberofparentsvec,n,updatenodes=c(1:n)) {
  maps<-list()
  mapi<-list()
  
  for (i in updatenodes) {
    nrows<-nrow(parenttable[[i]])
    P_local <- numeric(nrows)
    P_localinv <- numeric(nrows)
    
    P_local[1]<-1
    P_localinv[1]<-1
    if (nrows>1){
      for (j in 2:nrows)  {
        parentnodes <- parenttable[[i]][j,c(1:numberofparentsvec[[i]][j])]
        #numberofparentsvec stores number of non zero entries in i-th row in a parnttable,
        P_local[j]<-sum(2^parentnodes)/2+1
        #so extracting those non-zero entries we get row index
        P_localinv[P_local[j]]<-j # the inverse mapping
      }
    }
    mapi$forward<-P_local
    mapi$backwards<-P_localinv
    maps[[i]]<- mapi
  }
  
  return(maps)
}

poset<-function(parenttable,numberofparentsvec,rowmaps,n,updatenodes=c(1:n)){
  posetparenttables<-list(length=n)
  for (i in updatenodes) {
    
    nrows<-nrow(parenttable[[i]])
    ncols<-ncol(parenttable[[i]])
    posetparenttables[[i]]<-matrix(NA,nrow=nrows,ncol=ncols)
    offsets<-rep(1,nrows)
    
    if(nrows>1) {
      for(j in nrows:2){
        
        parentnodes<- parenttable[[i]][j,c(1:numberofparentsvec[[i]][j])]
        children<-rowmaps[[i]]$backwards[rowmaps[[i]]$forward[j]-2^parentnodes/2]
        posetparenttables[[i]][cbind(children,offsets[children])]<-j
        offsets[children]<-offsets[children]+1
      }
    }
    
  }
  return(posetparenttables)
}
