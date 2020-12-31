#samples maximum scoring DAG from a given order, returns incidence matrix of this DAG
#added samplenodes for background nodes, so we sample only for nodes which can have parents
samplescoreplus1.max<-function(matsize,samplenodes,scores,plus1lists,maxmatrices,scoretable,parenttable,
                               numberofparentsvec,aliases) {
  incidence<-matrix(0,ncol=matsize,nrow=matsize) # store the adjacency matrix
  sampledscore<-0
  for (i in samplenodes){
    if(is.null(plus1lists)) {
      krow<-maxmatrices$maxrow[[i]][scores$therow[i]]
      parentset<-aliases[[i]][parenttable[[i]][krow,1:numberofparentsvec[[i]][krow]]] #take right parent set
      incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
      sampledscore<-sampledscore+scoretable[[i]][krow]
    } else {
      klist<-which.max(maxmatrices$maxmatrix[[i]][scores$therow[i],scores$allowedlists[[i]]])
      abslist<-scores$allowedlists[[i]][klist]
      krow<-maxmatrices$maxrow[[i]][scores$therow[i],abslist]
      parentrow<-plus1lists$aliases[[i]][abslist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)] #take right parent set
      parentset<-parentrow[which(parentrow>0)] # removing NAs
      incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
      sampledscore<-sampledscore+scoretable[[i]][[abslist]][krow]
    }
  }
  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)
}

#samples a DAG from a given order, returns incidence matrix of this DAG
samplescoreplus1<-function(matsize,samplenodes,scores,plus1lists,scoretable,bannedscores,parenttable,
                           numberofparentsvec,aliases) {
  incidence<-matrix(numeric(matsize*matsize),nrow=matsize) # store the adjacency matrix
  sampledscore<-0
  for (i in samplenodes){
    if (scores$therow[i]==1) {
      allowedrows<-c(1:nrow(parenttable[[i]]))
      } else {
          bannednodes<-parenttable[[i]][scores$therow[i],1:numberofparentsvec[[i]][scores$therow[i]]]
          tablesize<-dim(parenttable[[i]]) # just to remove some arguments
          if(tablesize[1]==1||length(bannednodes)==tablesize[2]){ # no parents are allowed
          allowedrows<-c(1) # there is only one score
            } else {
              allowedrows<-c(2:tablesize[1])
              for (j in 1:tablesize[2]) {
                # working columnwise allows R to speed up
                bannedrows<-which(parenttable[[i]][allowedrows,j]%in%bannednodes)
                if(length(bannedrows)>0) { allowedrows<-allowedrows[-bannedrows] }
              }
            allowedrows<-c(1,allowedrows)
        }
      }
    if(is.null(plus1lists)) {
      
      allowedscores<-matrix(scoretable[[i]][allowedrows,1],nrow=length(allowedrows))
      k<-sample.int(length(allowedscores),1,prob=exp(allowedscores-scores$totscores[i]))
      krow<-allowedrows[k]
      parentset<-aliases[[i]][parenttable[[i]][krow,which(parenttable[[i]][krow,]>0)]] #take right parent set
      incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
      sampledscore<-sampledscore+scoretable[[i]][krow,1]
      
    } else {
              #first sample plus1 list
              klist<-sample.int(length(scores$allowedlists[[i]]),1,
                                      prob=exp(bannedscores[[i]][scores$therow[i],scores$allowedlists[[i]]]-scores$totscores[i]))
              abslist<-scores$allowedlists[[i]][klist]
              #then sample a row from allowed rows
              k<-sample.int(length(allowedrows),1,prob=exp(scoretable[[i]][[abslist]][allowedrows]-scores$totscores[i]))
              krow<-allowedrows[k]
              parentrow<-plus1lists$aliases[[i]][abslist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
              parentset<-parentrow[which(parentrow>0)] # removing NAs
              incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
              sampledscore<-sampledscore+scoretable[[i]][[abslist]][krow]
            }
    }

  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)
}

#samples a DAG from a given partition, returns incidence matrix of this DAG
samplescore.partition.plus1<-function(matsize,samplenodes,scores,scoretable,scoresallowed,scoresneeded,scoretab,parenttable,
                                      needednodetable,numberofparentsvec,
                                      needednodebannedrow,numberofpartitionparentsvec,
                                      plus1lists) {
  incidence<-matrix(numeric(matsize*matsize),nrow=matsize) # store the adjacency matrix
  sampledscore<-0
  logscorevec<-vector(length=matsize)
  for (i in samplenodes) {
    tablesize<-dim(parenttable[[i]]) # just to remove some arguments
    if (scores$neededrow[i]==0&&scores$allowedrow[i]==0){ #no parents are allowed
      sampledscore<-sampledscore+scoretable[[i]][[1]][1,1]
        } else if (scores$neededrow[i]==-1) { #bgnodes + last partition
           bannednodes<-parenttable[[i]][scores$allowedrow[i],1:numberofparentsvec[[i]][scores$allowedrow[i]]]
           if(tablesize[1]==1||length(bannednodes)==tablesize[2]) {#no allowed paernts in main table
             allowedscores<- scoretab[[i]][1,scores$plus1neededlists[[i]]]
             k<-sample.int(length(allowedscores),1,prob=exp(allowedscores-scores$totscores[i]))
             klist<-scores$plus1neededlists[[i]][k]
             sampledscore<-sampledscore+scoretable[[i]][[klist]][1]
             } else {
                  allowedrows<-c(2:tablesize[1])
                  for (j in 1:tablesize[2]) {
                   # working columnwise allows R to speed up
                  bannedrows<-which(parenttable[[i]][allowedrows,j]%in%bannednodes)
                  if(length(bannedrows)>0) { allowedrows<-allowedrows[-bannedrows] }
                    }
                  allowedrows<-c(1,allowedrows)
                  klist<-sample.int(length(scores$plus1neededlists[[i]]),1,
                                    prob=exp(scoresallowed[[i]][scores$allowedrow[i],scores$plus1neededlists[[i]]]-scores$totscores[i]))
                  abslist<-scores$plus1neededlists[[i]][klist]
                  #then sample a row from allowed rows
                  k<-sample.int(length(allowedrows),1,prob=exp(scoretable[[i]][[abslist]][allowedrows]-scores$totscores[i]))
                  krow<-allowedrows[k]
                  parentrow<-plus1lists$aliases[[i]][abslist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
                  parentset<-parentrow[which(parentrow>0)] # removing NAs
                  incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
                  sampledscore<-sampledscore+scoretable[[i]][[abslist]][krow]
                }
            } else {  
                if (scores$allowedrow[i]==1){
                    allowedrows2<-c(1:nrow(parenttable[[i]]))
                  } else {
                      bannednodes<-parenttable[[i]][scores$allowedrow[i],1:numberofparentsvec[[i]][scores$allowedrow[i]]]
                      if(tablesize[1]==1||length(bannednodes)==tablesize[2]){ # no parents are allowed
                        allowedrows2<-c(1) # there is only one score
                          } else {
                              allowedrows2<-c(2:tablesize[1])
                              for (j in 1:tablesize[2]) { # working columnwise allows R to speed up
                                bannedrows<-which(parenttable[[i]][allowedrows2,j]%in%bannednodes)
                                if(length(bannedrows)>0){
                                  allowedrows2<-allowedrows2[-bannedrows]
                                  }
                                }
                              allowedrows2<-c(1,allowedrows2)
                            }
                      }
    if (scores$neededrow[i]==0) {
      klist2<-sample.int(length(scores$plus1neededlists[[i]]),1,
                         prob=exp(scoresallowed[[i]][scores$allowedrow[i],scores$plus1neededlists[[i]]]-scores$totscores[i]))
      abslist2<-scores$plus1neededlists[[i]][klist2]
      k2<-sample.int(length(allowedrows2),1,prob=exp(scoretable[[i]][[abslist2]][allowedrows2]-scores$totscores[i]))
      krow2<-allowedrows2[k2]
      parentrow<-plus1lists$aliases[[i]][abslist2,c(1,parenttable[[i]][krow2,!is.na(parenttable[[i]][krow2,])]+1)]
      parentset<-parentrow[which(parentrow>0)] # removing NAs
      incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
      sampledscore<-sampledscore+scoretable[[i]][[abslist2]][krow2]
    } else {
      requirednodes<-needednodetable[[i]][scores$neededrow[i],1:(numberofpartitionparentsvec[[i]][scores$neededrow[i]])]
      bannedrow<-needednodebannedrow[[i]][scores$neededrow[i]]
      bannednodes<-parenttable[[i]][bannedrow,1:numberofparentsvec[[i]][bannedrow]]
      allowedrows1<-c(2:tablesize[1])
      if(bannedrow>1) { 
        for (j in 1:tablesize[2]) { # working columnwise allows R to speed up
        bannedrows<-which(parenttable[[i]][allowedrows1,j]%in%bannednodes)
        if(length(bannedrows)>0) {
          allowedrows1<-allowedrows1[-bannedrows]
          }
        }
      }
      notrequiredrows<-allowedrows1
      for (j in 1:tablesize[2]) { # now we remove the allowable rows instead
        requiredrows<-which(parenttable[[i]][notrequiredrows,j]%in%requirednodes)
        if(length(requiredrows)>0) {
          notrequiredrows<-notrequiredrows[-requiredrows]
        }
      }
      allowedscores2<-matrix(scoretab[[i]][allowedrows2,scores$plus1neededlists[[i]]],nrow=length(allowedrows2))
      allowedtablesize2<-dim(allowedscores2)
      allowedrows1<-setdiff(allowedrows1,notrequiredrows) # and keep just the difference!
      allowedscores1<-matrix(scoretab[[i]][allowedrows1,scores$plus1allowedlists[[i]]],nrow=length(allowedrows1))
      allowedtablesize1<-dim(allowedscores1)
      kboarder<-allowedtablesize2[1]*allowedtablesize2[2]
      k<-sample.int(length(allowedscores2)+length(allowedscores1),1,prob=c(exp(allowedscores2-scores$totscores[i]),exp(allowedscores1-scores$totscores[i])))
      if(k>kboarder) {
        colnumber<-ceiling((k-kboarder)/allowedtablesize1[1])
        rownumber<-k-kboarder-(colnumber-1)*allowedtablesize1[1]
        krow<-allowedrows1[rownumber]
        klist<-scores$plus1allowedlists[[i]][colnumber]
        parentrow<-plus1lists$aliases[[i]][klist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
        parentset<-parentrow[which(parentrow>0)] # removing NAs
        incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
        sampledscore<-sampledscore+scoretable[[i]][[klist]][krow]
      }
      else {
        colnumber<-ceiling(k/allowedtablesize2[1])
        rownumber<-k-(colnumber-1)*allowedtablesize2[1]
        krow<-allowedrows2[rownumber]
        klist<-scores$plus1neededlists[[i]][colnumber]
        parentrow<-plus1lists$aliases[[i]][klist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
        parentset<-parentrow[which(parentrow>0)] # removing NAs
        incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
        sampledscore<-sampledscore+scoretable[[i]][[klist]][krow]
      }
    }
    }
  }
  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)

}

#old version of function, kept for just in case
samplescore.partition.plus1.old<-function(n,scores,scoretable,scoresallowed,scoresneeded,scoretab,parenttable,needednodetable,numberofparentsvec,
                                          needednodebannedrow,numberofpartitionparentsvec,plus1lists){
  incidence<-matrix(numeric(n*n),nrow=n) # store the adjacency matrix
  sampledscore<-0
  logscorevec<-vector(length=n)
  for (i in 1:n)
  {
    if (scores$therow1[i]==0&&scores$therow2[i]==0){
      sampledscore<-sampledscore+scoretable[[i]][[1]][1,1]
    } else { tablesize<-dim(parenttable[[i]]) # just to remove some arguments
    if (scores$therow2[i]==1){
      allowedrows2<-c(1:nrow(parenttable[[i]]))
    }
    else{
      bannednodes<-parenttable[[i]][scores$therow2[i],1:numberofparentsvec[[i]][scores$therow2[i]]]
      if(tablesize[1]==1||length(bannednodes)==tablesize[2]){ # no parents are allowed
        allowedrows2<-c(1) # there is only one score
      } else{
        allowedrows2<-c(2:tablesize[1])
        for (j in 1:tablesize[2])
        { # working columnwise allows R to speed up
          bannedrows<-which(parenttable[[i]][allowedrows2,j]%in%bannednodes)
          if(length(bannedrows)>0){
            allowedrows2<-allowedrows2[-bannedrows]
          }
        }
        allowedrows2<-c(1,allowedrows2)
      }
    }
    if (scores$therow1[i]==0) {
      klist2<-sample.int(length(scores$allowedlists2[[i]]),1,
                         prob=exp(scoresneeded[[i]][scores$therow2[i],scores$allowedlists2[[i]]]-scores$totscores[i]))
      abslist2<-scores$allowedlists2[[i]][klist2]
      k2<-sample.int(length(allowedrows2),1,prob=exp(scoretable[[i]][[abslist2]][allowedrows2]-scores$totscores[i]))
      krow2<-allowedrows2[k2]
      parentrow<-plus1lists$aliases[[i]][abslist2,c(1,parenttable[[i]][krow2,!is.na(parenttable[[i]][krow2,])]+1)]
      parentset<-parentrow[which(parentrow>0)] # removing NAs
      incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
      sampledscore<-sampledscore+scoretable[[i]][[abslist2]][krow2]
    } else {
      requirednodes<-needednodetable[[i]][scores$therow1[i],1:(numberofpartitionparentsvec[[i]][scores$therow1[i]])]
      bannedrow<-needednodebannedrow[[i]][scores$therow1[i]]
      bannednodes<-parenttable[[i]][bannedrow,1:numberofparentsvec[[i]][bannedrow]]
      allowedrows1<-c(2:tablesize[1])
      if(bannedrow>1) { 
        for (j in 1:tablesize[2]) { # working columnwise allows R to speed up
          bannedrows<-which(parenttable[[i]][allowedrows1,j]%in%bannednodes)
          if(length(bannedrows)>0) {
            allowedrows1<-allowedrows1[-bannedrows]
          }
        }
      }
      notrequiredrows<-allowedrows1
      for (j in 1:tablesize[2]) { # now we remove the allowable rows instead
        requiredrows<-which(parenttable[[i]][notrequiredrows,j]%in%requirednodes)
        if(length(requiredrows)>0) {
          notrequiredrows<-notrequiredrows[-requiredrows]
        }
      }
      allowedscores2<-matrix(scoretab[[i]][allowedrows2,scores$allowedlists2[[i]]],nrow=length(allowedrows2))
      allowedtablesize2<-dim(allowedscores2)
      allowedrows1<-setdiff(allowedrows1,notrequiredrows) # and keep just the difference!
      allowedscores1<-matrix(scoretab[[i]][allowedrows1,scores$allowedlists1[[i]]],nrow=length(allowedrows1))
      allowedtablesize1<-dim(allowedscores1)
      kboarder<-allowedtablesize2[1]*allowedtablesize2[2]
      k<-sample.int(length(allowedscores2)+length(allowedscores1),1,prob=c(exp(allowedscores2-scores$totscores[i]),exp(allowedscores1-scores$totscores[i])))
      if(k>kboarder) {
        colnumber<-ceiling((k-kboarder)/allowedtablesize1[1])
        rownumber<-k-kboarder-(colnumber-1)*allowedtablesize1[1]
        krow<-allowedrows1[rownumber]
        klist<-scores$allowedlists1[[i]][colnumber]
        parentrow<-plus1lists$aliases[[i]][klist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
        parentset<-parentrow[which(parentrow>0)] # removing NAs
        incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
        sampledscore<-sampledscore+scoretable[[i]][[klist]][krow]
      }
      else {
        colnumber<-ceiling(k/allowedtablesize2[1])
        rownumber<-k-(colnumber-1)*allowedtablesize2[1]
        krow<-allowedrows2[rownumber]
        klist<-scores$allowedlists2[[i]][colnumber]
        parentrow<-plus1lists$aliases[[i]][klist,c(1,parenttable[[i]][krow,!is.na(parenttable[[i]][krow,])]+1)]
        parentset<-parentrow[which(parentrow>0)] # removing NAs
        incidence[parentset,i]<-1 # fill in elements of the adjacency matrix
        sampledscore<-sampledscore+scoretable[[i]][[klist]][krow]
      }
      
      
    }
    }
  }
  DAG<-list()
  DAG$incidence<-incidence
  DAG$logscore<-sampledscore
  return(DAG)
}
