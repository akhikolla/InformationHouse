#'Deriving an adjacency matrix of a graph
#'
#'This function derives the adjacency matrix corresponding to a graph object
#'
#'@param g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#'@return a square matrix whose dimensions are the number of nodes in the graph g, where element
#' \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j} in the graph \code{g},
#'  and \code{0} otherwise
#' @examples 
#' Asiagraph<-m2graph(Asiamat)
#' Asia.adj<-graph2m(Asiagraph)
#'@export
graph2m<-function(g) {
  l<-length(g@edgeL)
  adj<-matrix(rep(0,l*l),nrow=l,ncol=l)
  for (i in 1:l) {
    adj[g@edgeL[[i]]$edges,i]<-1}
  rownames(adj)<-g@nodes
  colnames(adj)<-g@nodes
  return(t(adj))
}


#'Deriving a graph from an adjacancy matrix
#'
#'This function derives a graph object corresponding to an adjacency matrix
#'
#'@param adj square adjacency matrix with elements in \code{\{0,1\}}, representing a graph
#'@param nodes (optional) labels of the nodes, \code{c(1:n)} are used by default
#'@return object of class \code{\link[graph]{graphNEL}} (package `graph'); if element \code{adj[i,j]} equals \code{1}, then there is a directed edge from node \code{i} to node \code{j} in the graph, and no edge otherwise
#'@examples 
#'m2graph(Asiamat)
#'@export
m2graph<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (is.null(nodes)) {
    if (!all(is.character(colnames(adj)))) {
    V <- c(1:l)
    edL <- vector("list", length=l)
    names(edL) <- sapply(V,toString)
    } else {
      edL <- vector("list", length=l)
      names(edL) <- colnames(adj)
      V<-colnames(adj)
    }
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  gR <- new("graphNEL", nodes=sapply(V,toString), edgeL=edL,edgemode="directed")
  
  return(gR)
  
}

#'Comparing two DAGs
#'
#'This function compares one (estimated) DAG to another DAG (true DAG), returning a vector of 3 values: structural Hamming distance,
#'number of true positive edges and number of false positive edges.
#'
#'@param eDAG an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the DAG which should be compared to a ground truth DAG or an ajacency matrix corresponding to the DAG
#'@param trueDAG an object of class \code{\link[graph]{graphNEL}} (package `graph'), representing the ground truth DAG or an ajacency matrix corresponding to the DAG
#'@return a vector of 5: SHD, number of true positive edges, number of false positive edges, number of false negative edges and true positive rate
#'@examples
#'Asiascore<-scoreparameters(8,"bde",Asia)
#'\dontrun{
#'eDAG<-orderMCMC(Asiascore)
#'compareDAGs(eDAG$max$DAG,Asiamat)
#'}
#'@export
compareDAGs<-function(eDAG,trueDAG) {
  
   if(is.matrix(eDAG)) eDAG<-m2graph(eDAG)
   if(is.matrix(trueDAG)) trueDAG<-m2graph(trueDAG)
   
   skeleton1<-graph2skeleton(eDAG)
   skeleton2<-graph2skeleton(trueDAG)
   
  numedges1<-sum(skeleton1)
  numedges2<-sum(skeleton2)
  
  diff2<-skeleton2-skeleton1
  res<-vector()
  res[1]<-pcalg::shd(eDAG, trueDAG)
  res[2]<-numedges1-sum(diff2<0) #TP
  res[3]<-sum(diff2<0) #FP
  res[4]<-sum(diff2>0)#FN
  res[5]<-res[2]/numedges2 #TPR
  attr(res, "class")<-"compDAGs"
  return(res)
}

#'Comparing two DBNs
#'
#'This function compares one (estimated) DBN structure to another DBN (true DBN). Comparisons for initial and transitional structures are returned separately if \code{equalstruct} equals \code{TRUE}.
#'
#'@param eDBN an object of class \code{\link[graph]{graphNEL}} (or an ajacency matrix corresponding to this DBN), representing the DBN which should be compared to a ground truth DBN 
#'@param trueDBN an object of class \code{\link[graph]{graphNEL}} (or an ajacency matrix corresponding to this DBN), representing the ground truth DBN 
#'@param struct option used to determine if the initial or the transitional structure should be compared; accaptable values are init or trans
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@return a vector of 5: SHD, number of true positive edges, number of false positive edges, number of false negative edges and true positive rate
#'@examples
#'testscore<-scoreparameters(15, "bge", DBNdata, bgnodes=c(1,2,3), DBN=TRUE,
#'                           dbnpar=list(slices=5, stationary=TRUE))
#'\dontrun{
#'DBNfit<-iterativeMCMC(testscore,moveprobs=c(0.11,0.84,0.04,0.01))
#'compareDBNs(DBNfit$max$DAG,DBNmat,struct="trans",n.dynamic=12,n.static=3)
#'}
#'@export
compareDBNs<-function(eDBN,trueDBN,struct=c("init","trans"),n.dynamic,n.static) {
  n<-n.static+n.dynamic
  matsize<-n.static+2*n.dynamic
  
  if(!is.matrix(eDBN)) {
    adj1<-graph2m(eDBN)
  } else {
    adj1<-eDBN
  }
  
  if(!is.matrix(trueDBN)) {
    adj2<-graph2m(trueDBN)
  } else {
    adj2<-trueDBN
  }
  
  adj1<-adj1[1:matsize,1:matsize]
  adj2<-adj2[1:matsize,1:matsize]
  
  skeleton1<-graph2skeleton(eDBN)
  skeleton2<-graph2skeleton(trueDBN)
  
  if(struct=="trans") {
  
    skeleton1[,1:n]<-0
    skeleton2[,1:n]<-0
    
    adj1[,1:n]<-0
    adj2[,1:n]<-0
    
  } else if (struct=="init"){
    
    adj1<-adj1[1:n,1:n]
    adj2<-adj2[1:n,1:n]

    skeleton1<-skeleton1[1:n,1:n]
    skeleton2<-skeleton2[1:n,1:n]
    
    
  } else {
    stop("parameter struct should equal either init or trans!")
  }
  eDAG<-m2graph(adj1)
  trueDAG<-m2graph(adj2)
  
  numedges1<-sum(skeleton1)
  numedges2<-sum(skeleton2)
  
  diff2<-skeleton2-skeleton1
  res<-vector()
  res[1]<-pcalg::shd(eDAG, trueDAG)
  res[2]<-numedges1-sum(diff2<0) #TP
  res[3]<-sum(diff2<0) #FP
  res[4]<-sum(diff2>0)#FN
  res[5]<-res[2]/numedges2 #TPR
  attr(res, "class")<-"compDAGs"
  return(res)
}


#returns a matrix of a CPDAG corresponding to a given DAG
dagadj2cpadj<-function(adj) {
  g<-m2graph(adj)
  cpg<-pcalg::dag2cpdag(g)
  return(graph2m(cpg))
}

#returns a symmetric matrix of a skeleton corresponding to a given CPDAG
#UPPER TRIANGULAR VERSION!!!!
adjacency2skeleton<-function(adj) {
  skel<-1*(adj|t(adj))
  skel<-ifelse(upper.tri(skel)==TRUE,skel,0)
  return(skel)
}

#returns a list of edges corresponding to an adjacency matrix
adjacency2edgel<-function(adj,nodes=NULL) {
  l<-ncol(adj)
  if (is.null(nodes)) {
    if (!all(is.character(colnames(adj)))) {
      V <- c(1:l)
      edL <- vector("list", length=l)
      names(edL) <- sapply(V,toString)
    } else {
      edL <- vector("list", length=l)
      names(edL) <- colnames(adj)
      V<-colnames(adj)
    }
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))}
  else {
    V <- nodes
    edL <- vector("list", length=l)
    names(edL) <- V
    for(i in 1:l)
      edL[[i]] <- list(edges=which(adj[i,]==1))
  }
  return(edL)
  
}


#Deriving an adjacency matrix of the skeleton of a graph
#
#This function derives the skeleton matrix corresponding to a graph object
#
#g graph, object of class \code{\link[graph]{graphNEL}} (package `graph')
#returns a symmetric square matrix whose dimensions are the number of nodes in the graph \code{g},
# where element \code{[i,j]} equals \code{1} if there is a directed edge from node \code{i} to node \code{j},
#  or from node \code{j} to node \code{i}, in the graph \code{g}, and \code{0} otherwise
# examples 
#myDAG<-pcalg::randomDAG(20, prob=0.15, lB = 0.4, uB = 2)
#graph2skeleton.m(myDAG)
graph2skeleton<-function(g,upper=TRUE,outmat=TRUE) {
  if(!is.matrix(g)) {
    adj<-graph2m(g)
    colnames(adj)<-g@nodes
    rownames(adj)<-g@nodes
  } else {
    adj<-g
  }
  skel<-1*(adj|t(adj))
  if(upper) {
  skel<-ifelse(upper.tri(skel)==TRUE,skel,0)
  }
  if(outmat) {
    return(skel)
  } else {
    return(m2graph(skel))
  }
}



