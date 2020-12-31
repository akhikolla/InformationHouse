# this function scores by propagating through the partition poset graph
posetpartitionscoreparents<-function(numberofparents,neededposetparenttable,
                                     neededparentsvec,numberofparentsvec,rowmaps,
                                     needednodebannedrow,scoretable,n,updatenodes=c(1:n)) {
  
  scorematrices<-list()

  for (j in updatenodes)  {
    np<-numberofparents[j]
    if(np==1) {
      scorematrices[[j]]<-as.matrix(scoretable[[j]][2,1])
    } else if (np>1){
      binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
      nrows<-nrow(neededposetparenttable[[j]])
      P_local <- numeric(nrows)
      nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph
      # at the top of the graph we only allow one possible parent set

      for (i in nrows:(nrows-np+1)) {
        k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
        conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
        P_local[i]<-conjugatescore
      }

      # for each other level of the poset graph we need to add the parents divided by the difference in levels

      cutoff<-1

      for(level in 0:(np-2)) {

        for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
          k <- needednodebannedrow[[j]][i] # the banned nodes row
          parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
          maxparents<-max(P_local[parentnodes])
          # take the sum of the parent scores and divide by the relevant factor
          parentsum<-log(sum(exp(P_local[parentnodes]-maxparents)))+maxparents -
            log(np-rev(numberofparentsvec[[j]])[k]-level+1)
          # take the conjugate score
          conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
          # find max and combine
          maxoverall<-max(parentsum,conjugatescore)
          P_local[i]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
        }

        cutoff<-cutoff+binomcoefs[level+1]
      }
      scorematrices[[j]]<-as.matrix(P_local)
    }
  }
  return(scorematrices)
}
#builds a table where entries correspond to partition scores of required and allowed nodes
plus1needed.partition<-function(numparents,parenttable,neededposetparenttable,neededparentsvec,
                                 numberofparentsvec,rowmapsneeded,needednodebannedrow,
                                scoretable,plus1lists,n,updatenodes=c(1:n)) {
  
  revnumberofparentsvec<-lapply(numberofparentsvec,rev)
  rowmaps<-rowmapsneeded
  if (is.null(plus1lists)) {
    scorematrices<-list()
    for (j in updatenodes){
      np<-numparents[j]
      if(np==1) {
        scorematrices[[j]]<-as.matrix(scoretable[[j]][2,1])
      } else if (np>1) {
        binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
        nrows<-nrow(neededposetparenttable[[j]])
        P_local <- numeric(nrows)
        nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph
        # at the top of the graph we only allow one possible parent set
        for (i in nrows:(nrows-np+1)) {
          k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
          conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
          P_local[i]<-conjugatescore
        }
        # for each other level of the poset graph we need to add the parents divided by the difference in levels
        cutoff<-1
        for(level in 0:(np-2)) {
          for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
            k <- needednodebannedrow[[j]][i] # the banned nodes row
            parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
            maxparents<-max(P_local[parentnodes])
            # take the sum of the parent scores and divide by the relevant factor
            parentsum<-log(sum(exp(P_local[parentnodes]-maxparents)))+maxparents - log(np-revnumberofparentsvec[[j]][k]-level+1)
            # take the conjugate score
            conjugatescore<-scoretable[[j]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
            # find max and combine
            maxoverall<-max(parentsum,conjugatescore)
            P_local[i]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
          }
          cutoff<-cutoff+binomcoefs[level+1]
        }
        scorematrices[[j]]<-as.matrix(P_local)
      }
    }
    return(scorematrices)
  } else{
    scorematrices.needed<-list()
    for (j in updatenodes){
      np<-numparents[j] #number of possible parents for node j
      binomcoefs<-choose(np,c(np:1))*(2^c(np:1)-1)
      ll<-length(plus1lists$parents[[j]])+1
      nrows<-nrow(neededposetparenttable[[j]])
      nrowspar<-nrow(parenttable[[j]])
      if(np>0) {P_local <- matrix(nrow=nrows,ncol=ll)}
      for (li in 1:ll){
        if (np==0){
          scorematrices.needed[[j]]<-NULL
          break #we don't have allowed table, all scores already in the scoretable plus1 lists
        }
        else if(np==1) {
          #we need just 1 additional parent set which is not in the scoretable plus1 node + the only parent in scoretable
          P_local[1,li]<-as.matrix(scoretable[[j]][[li]][2,1])
          scorematrices.needed[[j]]<-P_local
        } else if (np>1){
          nrowold<-length(rowmaps[[j]]$forward) # size of the other poset graph

            for (i in nrows:(nrows-np+1)) {
              k <- needednodebannedrow[[j]][i] # the banned nodes row, there should be (n-1) banned nodes
              conjugatescore<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
              P_local[i,li]<-conjugatescore
          }

          # for each other level of the poset graph we need to add the parents divided by the difference in levels
          cutoff<-1
          for(level in 0:(np-2)){
            for (i in (nrows-np):cutoff)  { # find the parents in the poset graph
              k <- needednodebannedrow[[j]][i] # the banned nodes row
              parentnodes <- neededposetparenttable[[j]][i,c(1:neededparentsvec[[j]][i])]
              maxparents<-max(P_local[parentnodes,li])
              # take the sum of the parent scores and divide by the relevant factor
              parentsum<-log(sum(exp(P_local[parentnodes,li]-maxparents)))+maxparents-log(np-revnumberofparentsvec[[j]][k]-level+1)
              # take the conjugate score
              conjugatescore<-scoretable[[j]][[li]][rowmaps[[j]]$backwards[nrowold-rowmaps[[j]]$forward[k]+1],1]
              # find max and combine
              maxoverall<-max(parentsum,conjugatescore)
              P_local[i,li]<- log(exp(parentsum-maxoverall)+exp(conjugatescore-maxoverall)) + maxoverall
            }

            cutoff<-cutoff+binomcoefs[level+1]
          }
          scorematrices.needed[[j]]<-as.matrix(P_local)
        }

      }
    }
    return(scorematrices.needed)
  }

}

parentlistnonempty<-function(elements,n){
  
  matrixofparents<-matrix(NA,nrow=(2^length(elements)-1),ncol=n)
  
  cutoff<-0
  
  for (r in 1:length(elements)){
    possparents<-combinations(length(elements),r,elements)
    heighty<-nrow(possparents)
    matrixofparents[(1:heighty)+cutoff,1:r]<-possparents
    cutoff<-cutoff+heighty
  }
  return(matrixofparents)
}

partitionlistnumberofparents<-function(needednodetable,numberofparentsvec,n,updatenodes){
  numberofpartitionparentsvec<-list()
  
  for (i in updatenodes){
    
    if(length(numberofparentsvec[[i]])==1) {numberofpartitionparentsvec[[i]]<-0}
    else {
      nrows<-length(numberofparentsvec[[i]])
      npar<-numberofparentsvec[[i]][nrows]
      P_local<-c()
      for (j in 1:(nrows-1))  { # last row has no possibilities
        np<-npar-numberofparentsvec[[i]][j]
        P_local<-c(P_local,rep(c(1:np),choose(np,c(1:np))))
      }
      numberofpartitionparentsvec[[i]]<-P_local
    }
  }
  return(numberofpartitionparentsvec)
}

partitionmapneedednodebannedrow<-function(numberofparents,numberofparentsvec,n,updatenodes) {
  needednodebannedrow<-list()
  for (i in updatenodes) { 
    j<-numberofparents[i]
    needednodebannedrow[[i]]<-rep(1:2^j,2^(j-numberofparentsvec[[i]])-1)
  }
  return(needednodebannedrow)
}

partitionlist<-function(parenttable,numberofparentsvec,n,updatenodes) {
  needednodetable<-list()
  for  (i in updatenodes) {
    cutoff<-0
    nrows<-nrow(parenttable[[i]])
    ncols<-ncol(parenttable[[i]])
    matrixofparents<-matrix(NA,nrow=(3^ncols-2^ncols),ncol=ncols)
    if(nrows>1) {
      for (j in 1:(nrows-1)){ # last row has no possibilities
        parentnodes <- parenttable[[i]][j,1:numberofparentsvec[[i]][j]]
        elements<-setdiff(c(1:ncols),parentnodes)
        newpart<-parentlistnonempty(elements,ncols)
        heighty<-nrow(newpart)
        matrixofparents[(1:heighty)+cutoff,]<-newpart
        cutoff<-cutoff+heighty
      }
    }
    needednodetable[[i]]<-as.matrix(matrixofparents)
  }
  return(needednodetable)
}

neededparentsmapping<-function(parenttable,numberofparentsvec,needednodetable,
                               numberofpartitionparentsvec,needednodebannedrow,n,updatenodes) {
  
  maps<-list()
  for (i in updatenodes){
    nrows<-nrow(needednodetable[[i]])
    ncols<-ncol(needednodetable[[i]])
    maps[[i]]<-list()
    P_local <- numeric(3^ncols) # we leave some elements empty
    P_localinv <- numeric(3^ncols) # here too
    if(nrows>1){
      for (j in 1:nrows)  {# the needed nodes
        needednodes <- needednodetable[[i]][j,c(1:numberofpartitionparentsvec[[i]][j])]
        k <- needednodebannedrow[[i]][j] # the banned nodes row
        if(k>1) { # if there is at least one banned node
          bannednodes <- parenttable[[i]][k,c(1:numberofparentsvec[[i]][k])]
          P_local[j]<-sum(3^bannednodes)/3+2*sum(3^needednodes)/3+1
        } else {P_local[j]<-2*sum(3^needednodes)/3+1}
        P_localinv[P_local[j]]<-j # the inverse mapping
      }
    }
    maps[[i]]$forward<-P_local
    maps[[i]]$backwards<-P_localinv
  }
  return(maps)
}

needed.poset<-function(parenttable,numberofparentsvec,needednodebannedrow,neededrowmaps,n,
                       updatenodes=c(1:n)) { 
  posetparents<-list()
  for (i in updatenodes) {
    nrows<-length(needednodebannedrow[[i]])
    numpar<-numberofparentsvec[[i]][nrow(parenttable[[i]])]
    if(nrows>1) {
      posetneededtable<-matrix(NA,nrow=nrows,ncol=numpar)
      offsets<-rep(1,nrows)
      if(nrows>1) {
        for(j in (nrows:(2^numpar))){
          k <- needednodebannedrow[[i]][j] # the banned nodes row, there should be at least one banned node
          bannednodes <- parenttable[[i]][k,c(1:numberofparentsvec[[i]][k])]
          # we can either delete one of the banned nodes
          children1<-neededrowmaps[[i]]$backwards[neededrowmaps[[i]]$forward[j]-3^bannednodes/3]
          # move one of the banned nodes to a needed noded
          children2<-neededrowmaps[[i]]$backwards[neededrowmaps[[i]]$forward[j]+3^bannednodes/3]
          children<-c(children1,children2)
          posetneededtable[cbind(children,offsets[children])]<-j
          offsets[children]<-offsets[children]+1
        }
      }
    }
    posetparents[[i]]<-list()
    if (numpar==1) {
      posetneededtable<-matrix(NA,nrow=1,ncol=1)
      posetneededtable[1,1]<-1
      posetparents[[i]]$sizes<-c(1)
      posetparents[[i]]$table<-posetneededtable
    }
    else if (numpar>1) {
      posetparents[[i]]$sizes<-offsets-1
      posetparents[[i]]$table<-posetneededtable}
    else {
      posetparents[[i]]$sizes<-vector("integer",length=0)
      posetparents[[i]]$table<-NULL}
  }
  return(posetparents)
}

plus1allowed.partition<-function(posetparenttable,scoretable,numberofparentsvec,rowmaps,n,plus1lists=NULL,
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

