TableDAGscore.alias <- function(parentrows, j, n,alias,param,parentmaps=NULL,numparents=NULL,
                                numberofparentsvec=NULL) {
  
  if (param$type=="bde" & !param$DBN) {
    nrows<-nrow(parentrows)
    parentnodes<- alias[parentrows[nrows,!is.na(parentrows[nrows,])]]
    P_local<-DAGbinarytablescore(j,parentnodes,n,param,parentrows,parentmaps,numparents,numberofparentsvec)
  } else if (param$type=="bdecat" & !param$DBN) {
    nrows<-nrow(parentrows)
    parentnodes<- alias[parentrows[nrows,!is.na(parentrows[nrows,])]]
    P_local<-DAGcattablescore(j,parentnodes,n,param,parentrows,parentmaps,numparents,numberofparentsvec)
  } else {
    nrows<-nrow(parentrows)
    P_local <- numeric(nrows)
    for (i in 1:nrows)  {
      parentnodes <- alias[parentrows[i,!is.na(parentrows[i,])]]
      P_local[i]<-DAGcorescore(j,parentnodes,n,param)
    } 
  }
  
  return(P_local)
}

TableDAGscore.alias.plus1<-function(parentrows, j, n,alias,param,parentmaps=NULL,numparents=NULL,numberofparentsvec=NULL) {

  if (param$type=="bde" & !param$DBN) {
    nrows <- nrow(parentrows)
    parentnodes <- alias[parentrows[nrows,!is.na(parentrows[nrows,])]+1]
    addpar <- alias[1]
    P_local <- DAGbinarytablescoreplus1(j,parentnodes,addpar,n,param,parentrows,parentmaps,numparents,numberofparentsvec)
  } else if (param$type=="bdecat" & !param$DBN) {
    nrows <- nrow(parentrows)
    parentnodes <- alias[parentrows[nrows,!is.na(parentrows[nrows,])]+1]
    addpar <- alias[1]
    P_local <- DAGcattablescoreplus1(j,parentnodes,addpar,n,param,parentrows,parentmaps,numparents,numberofparentsvec)
  } else { 
    nrows<-nrow(parentrows)
    P_local <- numeric(nrows)
    for (i in 1:nrows)  {
      parentnodes <- alias[c(1,parentrows[i,!is.na(parentrows[i,])]+1)]
      P_local[i] <- DAGcorescore(j,parentnodes,n,param)
    } 
  }
  
  return(P_local)
}

listpossibleparents.PC.aliases<-function(skeletonedges,isgraphNEL=FALSE,n,updatenodes=c(1:n)){
  if(isgraphNEL==FALSE){
    l1<-ncol(skeletonedges)
      } else {l1<-length(skeletonedges)}
  listy<-vector("list",l1)
  aliases<-vector("list",l1)
  numparents<-vector("numeric",l1)

  #we keep record of which parent table lengths we already constructed
  table.with.k.parents<-matrix(rep(0,l1*2),nrow=2,ncol=l1)

  for (i in updatenodes){
    if (isgraphNEL==TRUE) {possparents<-skeletonedges[[i]]$edges
    } else{possparents<-which(skeletonedges[,i]==1)}
    aliases[[i]]<-possparents
    l<-length(possparents)
    numparents[i]<-l
    possparents<-c(1:l)
    if (l==0){
      matrixofparents<-matrix(rep(NA,1),1,1)
    } else if (table.with.k.parents[1,l]>0){
      matrixofparents<-listy[[table.with.k.parents[2,l]]]
    } else {
      matrixofparents<-rep(NA,l)
      for (r in 1:l){
        combpossparents<-combinations(l,r,possparents)
        if(r<l){
          for (j in 1:(l-r)){
            combpossparents <- cbind(combpossparents, NA)
          }
        }
        matrixofparents<-rbind(matrixofparents,combpossparents,deparse.level=0)
      }
    }
    listy[[i]] <- matrixofparents
    table.with.k.parents[1,l]<-1
    table.with.k.parents[2,l]<-i
  }
  listz<-list()
  listz$parenttable<-listy
  listz$aliases<-aliases
  listz$numparents<-numparents
  listz$numberofparentsvec<-lapply(numparents,function(x)rep(c(0:x),choose(x,c(0:x))))

  return(listz)
}

scorepossibleparents.alias<-function(parenttable,aliases,n,param,
                                     updatenodes=c(1:n),parentmaps=NULL,
                                     numparents=NULL,numberofparentsvec=NULL){

  listz<-vector("list",n)

  for (i in updatenodes) {
    scoretemp<-TableDAGscore.alias(parenttable[[i]], i, n,aliases[[i]],param,parentmaps[[i]],numparents[i],numberofparentsvec[[i]])
    listz[[i]] <- as.matrix(scoretemp)
  }
  return(listz)
}

PLUS1<-function(n,aliases,updatenodes=c(1:n),blacklistparents=NULL) {
  listz<-list()
  plus1mask<-list()
  plus1parents<-list()
  plus1aliases<-list()
  for (i in updatenodes){
    plus1mask[[i]]<-rep(1,n)
    plus1mask[[i]][c(aliases[[i]],blacklistparents[[i]])]<-0
    plus1parents[[i]]<-which(plus1mask[[i]]==1)
    nrows<-length(plus1parents[[i]])+1
    ncols<-length(aliases[[i]])+1
    plus1aliases[[i]]<-matrix(c(NaN,plus1parents[[i]],rep(aliases[[i]], each = nrows) ),
                              nrow=nrows,ncol=ncols)
  }
  listz$mask<-plus1mask
  listz$parents<-plus1parents
  listz$aliases<-plus1aliases

  return(listz)
}

scorepossibleparents.PLUS1<-function(parenttable,plus1lists,n,param,updatenodes,
                                     parentmaps,numparents,numberofparentsvec){

  listy<-vector("list",n)
  aliases<-plus1lists$aliases

  for (i in updatenodes){ #for every node which needs to be updated
    k<-nrow(aliases[[i]])
    ncols<-ncol(aliases[[i]])
    listz<-vector("list",k)

    for (j in 1:k){ #for every list
      if (j==1) {
        scoretemp<-TableDAGscore.alias(parenttable[[i]], i, n,aliases[[i]][j,which(!is.na(aliases[[i]][j,]))],param,parentmaps[[i]],numparents[i],
                                       numberofparentsvec[[i]])
      } else {
        scoretemp<-TableDAGscore.alias.plus1(parenttable[[i]], i, n,aliases[[i]][j,],param,parentmaps[[i]],numparents[i],numberofparentsvec[[i]])}
        listz[[j]] <- as.matrix(scoretemp)
    }
    listy[[i]]<-listz
  }
  return(listy)
}

bgNodeScore<-function(n,param) {
  totscores<-vector()
  for(i in param$bgnodes) {
    totscores[i]<-DAGcorescore(i,NULL,n,param)
  }
  return(totscores)
}



