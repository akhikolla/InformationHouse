#'Deriving an adjecency matrix of a full DBN
#'
#'This function transforms a compact 2-slice adjacency matrix of DBN into full T-slice adjacency matrix
#'
#' @param DBNmat a square matrix, representing initial and transitional structure of a DBN; the size of matrix is 2*n.dynamic+n.static
#' @param n.dynamic integer, number of dynamic variables in one time slice
#' @param n.slices integer, number of slices in an unrolled DBN
#' @param n.static integer, number of static variables
#' @return an adjacency matrix of an unrolled DBN
#' @examples
#' compact2full(DBNmat, n.dynamic=12, n.slices=5, n.static=3)
#' @export
compact2full<-function(DBNmat,n.dynamic,n.slices,n.static=0) {

  if(n.slices<3) {
    return(DBNmat)
  } else {

    if(all(is.character(colnames(DBNmat)))){
      baseall<-colnames(DBNmat)
      basenames<-colnames(DBNmat)[1:n.dynamic+n.static]
    } else {
      if(n.static!=0) {
        staticnames<-paste("s",1:n.static,sep="")
        basenames<-paste("v",1:n.dynamic,sep="")
        baseall<-c(staticnames,basenames)
      } else {
        basenames<-paste("v",1:n.dynamic,sep="")
        baseall<-basenames
      }
    }
    
    for(i in 3:n.slices) {
    baseall<-c(baseall,paste(basenames,".",i,sep=""))
    }
    
    nbig<-n.slices*n.dynamic+n.static
    DBNbig<-matrix(0,nrow=nbig,ncol=nbig)
    print(baseall)
    colnames(DBNbig)<-baseall
    rownames(DBNbig)<-baseall
    
    
    DBNbig[1:(n.dynamic+n.static),1:n.dynamic+n.static]<-DBNmat[1:(n.dynamic+n.static),1:n.dynamic+n.static] #copying initial structure
   
    intStruct<-DBNmat[1:n.dynamic+n.dynamic+n.static,1:n.dynamic+n.dynamic+n.static] #internal structure
    transStruct<-DBNmat[1:n.dynamic+n.static,1:n.dynamic+n.dynamic+n.static] #transitional structure
    if(n.static>0) {
      bgStrct<-DBNmat[1:n.static,1:n.dynamic+n.dynamic+n.static] #edges from static variables
    }
    for(i in 1:(n.slices-1)) {
      if(n.static>0) {
        DBNbig[1:n.static,1:n.dynamic+i*n.dynamic+n.static]<-bgStrct
      }
      DBNbig[1:n.dynamic+(i-1)*n.dynamic+n.static,1:n.dynamic+i*n.dynamic+n.static]<-transStruct
      DBNbig[1:n.dynamic+i*n.dynamic+n.static,1:n.dynamic+i*n.dynamic+n.static]<-intStruct
    }
    return(DBNbig)
  }
}
#'Deriving a compact adjacency matrix of a DBN 
#'
#'This function transforms an unrolled adjacency matrix of DBN into a compact representation
#'
#' @param DBNmat a square matrix, representing the structure of an unrolled DBN; the size of matrix is n.slices*n.dynamic+n.static; all static variables are assumed to be in the first n.static rows and columns of the matrix
#' @param n.dynamic integer, number of dynamic variables in each time slice
#' @param n.static integer, number of static variables
#' @examples
#' full2compact(DBNunrolled,n.dynamic=12,n.static=3)
#'@export
full2compact<-function(DBNmat,n.dynamic,n.static=0) {
    DBNcompact<-DBNmat[1:(2*n.dynamic+n.static),1:(2*n.dynamic+n.static)]
    return(DBNcompact)
}

#turns internal representation into user-friendly
DBNtransform<-function(DBNmat,param) {
  newDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
  colnames(newDBNmat)<-param$labels.short
  rownames(newDBNmat)<-param$labels.short
  newDBNmat[param$usrinitstr$rows,param$usrinitstr$cols]<-DBNmat[param$intstr$rows,param$intstr$cols]
  newDBNmat[param$usrintstr$rows,param$usrintstr$cols]<-DBNmat[param$intstr$rows,param$intstr$cols]
  newDBNmat[param$usrtrans$rows,param$usrtrans$cols]<-DBNmat[param$trans$rows,param$trans$cols]
  return(newDBNmat)
}
#turns internal representation into user-friendly
DBNtransform.init<-function(DBNmat,param) {
  if(param$bgn>0) {
      newDBNmat<-matrix(0,nrow=param$bgn+param$nsmall,ncol=param$bgn+param$nsmall)
      colnames(newDBNmat)<-param$labels.short[1:param$n]
      rownames(newDBNmat)<-param$labels.short[1:param$n]
      newDBNmat[,1:param$bgn]<-DBNmat[,1:param$bgn+param$nsmall]
      newDBNmat[,1:param$nsmall+param$bgn]<-DBNmat[,1:param$nsmall]
      DBNmat<-newDBNmat
      newDBNmat[1:param$bgn,]<-DBNmat[1:param$bgn+param$nsmall,]
      newDBNmat[1:param$nsmall+param$bgn,]<-DBNmat[1:param$nsmall,]
      return(newDBNmat)
  } else {
    return(DBNmat)
  }
}
#turns user-friendly representation into internal
DBNbacktransform<-function(DBNmat,param) {
  if(!is.null(colnames(DBNmat))) {
    oldnodelabels<-colnames(DBNmat)
    newnodelabels<-oldnodelabels
    newnodelabels[param$intstr$cols]<-oldnodelabels[param$usrtrans$cols]
  }
  newDBNmat<-matrix(0,nrow=param$n+param$nsmall,ncol=param$n+param$nsmall)
  newDBNmat[param$intstr$rows,param$intstr$cols]<-1*(DBNmat[param$usrintstr$rows,param$usrintstr$cols]|DBNmat[param$usrinitstr$rows,param$usrinitstr$cols])
  newDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
  if(!param$split) {
    return(newDBNmat)
  } else {
    res<-list()
    initDBNmat<-DBNmat[1:param$n,1:param$n]
    newinitDBNmat<-DBNmat[1:param$n,1:param$n]
    
    newinitDBNmat[,1:param$bgn+param$nsmall]<-initDBNmat[,1:param$bgn]
    newinitDBNmat[,1:param$nsmall]<-initDBNmat[,1:param$nsmall+param$bgn]
    initDBNmat<-newinitDBNmat
    newinitDBNmat[1:param$bgn+param$nsmall,]<-initDBNmat[1:param$bgn,]
    newinitDBNmat[1:param$nsmall,]<-initDBNmat[1:param$nsmall+param$bgn,]
    res$init<-newinitDBNmat
    
    transDBNmat<-matrix(0,nrow=param$otherslices$n,ncol=param$otherslices$n)
    DBNmat<-DBNcut(DBNmat,n.dynamic=param$nsmall,n.static=param$bgn)
    transDBNmat[param$intstr$rows,param$intstr$cols]<-DBNmat[param$usrintstr$rows,param$usrintstr$cols]
    transDBNmat[param$trans$rows,param$trans$cols]<-DBNmat[param$usrtrans$rows,param$usrtrans$cols]
    res$trans<-transDBNmat
    
    return(res)
  }
}
DBNcut<-function(adj,n.dynamic,n.static){
  adj[,1:(n.dynamic+n.static)]<-0
  return(adj)
}
DBNinit<-function(adj,n.dynamic,n.static){
  adj<-adj[1:(n.static+n.dynamic),1:(n.static+n.dynamic)]
  if(n.static>0) {
    adj[,1:n.static]<-0
  }
  return(adj)
}

#Combining initial and transition DBN structures in one matrix
mergeDBNstr<-function(initStruct,transStruct) {
  if(!is.matrix(initStruct)) {
    initStruct<-graph2m(initStruct)
  }
  if(!is.matrix(transStruct)) {
    transStruct<-graph2m(transStruct)
  }
  n<-ncol(initStruct)
  transStruct[1:n,1:n]<-initStruct
  return(transStruct)
}

#Combining orders for a DBN
mergeDBNord<-function(initorder,transorder) {
 return(c(initorder,transorder))
}

#Combining order scores for a DBN
mergeDBNscore<-function(initscore,transscore) {
  return(initscore+transscore)
}

#this function produces common result for DBN structure learning when samestruct=FALSE
mergeDBNres<-function(result.init,result.trans,scorepar,algo) {
  
  res<-list()
  
  maxtrans<-DBNtransform(result.trans$max$DAG,scorepar)
  maxinit<-DBNtransform.init(result.init$max$DAG,scorepar)
  res$max$DAG<-mergeDBNstr(maxinit,maxtrans)
  res$max$order<-mergeDBNord(result.init$max$order,result.trans$max$order)
  res$max$score<-mergeDBNscore(result.init$max$score,result.trans$max$score)
  
  if(!is.null(result.init$chain)) {
    
    result.init$chain$incidence<-lapply(result.init$chain$incidence,DBNtransform.init,param=scorepar)
    result.trans$chain$incidence<-lapply(result.trans$chain$incidence,DBNtransform,param=scorepar)
    result.trans$chain$incidence<-lapply(result.trans$chain$incidence,DBNcut,n.dynamic=scorepar$nsmall,n.static=scorepar$bgn)
    res$chain$incidence<-mapply(mergeDBNstr,result.init$chain$incidence,result.trans$chain$incidence,SIMPLIFY = FALSE)
    res$chain$DAGscores<-mapply(mergeDBNscore,result.init$chain$DAGscores,result.trans$chain$DAGscores)
    
  if(algo=="order"){
    res$chain$orders<-mapply(mergeDBNord,result.init$chain$orders,result.trans$chain$orders,SIMPLIFY = FALSE)
    res$chain$orderscores<-mapply(mergeDBNscore,result.init$chain$orderscores,result.trans$chain$orderscores)
  } else if (algo=="partition") {
    res$chain$order<-mapply(mergeDBNord,result.init$chain$order,result.trans$chain$order,SIMPLIFY = FALSE)
    res$chain$partitionscores<-mapply(mergeDBNscore,result.init$chain$partitionscores,result.trans$chain$partitionscores)
  }
  }
  
  attr(res,"class")<-"MCMCres"
  
  return(res)
}


#this function produces common result for DBN iterative structure learning when samestruct=FALSE
mergeDBNres.it<-function(result.init,result.trans,scorepar) {
  
  res<-list()
  
  res$init<-result.init
  res$trans<-result.trans
  
  maxtrans<-DBNtransform(result.trans$max$DAG,scorepar)
  maxinit<-DBNtransform.init(result.init$max$DAG,scorepar)
  
  
  for(i in 1:length(res$trans$maxtrace)) {
    res$trans$maxtrace[[i]]$DAG<-DBNtransform(res$trans$maxtrace[[i]]$DAG,scorepar)
    res$trans$maxtrace[[i]]$DAG<-DBNcut(res$trans$maxtrace[[i]]$DAG,n.dynamic=scorepar$nsmall,n.static=scorepar$bgn)
  }
  
  for(i in 1:length(res$init$maxtrace)) {
    res$init$maxtrace[[i]]$DAG<-DBNtransform.init(res$init$maxtrace[[i]]$DAG,scorepar)
    res$init$maxtrace[[i]]$DAG<-DBNinit(res$init$maxtrace[[i]]$DAG,n.dynamic=scorepar$nsmall,n.static=scorepar$bgn)
  }
  
  res$max$DAG<-mergeDBNstr(maxinit,maxtrans)
  res$max$order<-mergeDBNord(result.init$max$order,result.trans$max$order)
  res$max$score<-mergeDBNscore(result.init$max$score,result.trans$max$score)

  endinit<-DBNtransform.init(result.init$endspace,scorepar)
  endtrans<-DBNtransform(result.trans$endspace,scorepar)
  startinit<-DBNtransform.init(result.init$startspace,scorepar)
  starttrans<-DBNtransform(result.trans$startspace,scorepar)
    
  res$endspace<-mergeDBNstr(endinit,endtrans)
  res$startspace<-mergeDBNstr(startinit,starttrans)
  
  return(res)
}






