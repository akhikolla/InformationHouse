#' Plotting difference between two DBNs
#' 
#' This function plots an estimated DBN such that the edges which are different to the ground truth DBN are highlighted. 
#' 
#'@param eDBN object of class graphNEL (or its adjacency matrix), representing estimated structure (not necessarily acyclic) to be compared to the ground truth graph
#'@param trueDBN object of class graphNEL (or its adjacency matrix), representing the ground truth structure (not necessarily acyclic)
#'@param struct option used to determine if the initial or the transition structure should be plotted; accaptable values are init or trans
#'@param n.dynamic number of dynamic variables in one time slice of a DBN
#'@param n.static number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix
#'@return plots the graph which includes edges from edag and truedag, however edges which are different in edag compared to truedag are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'dbnscore<-scoreparameters(15,"bge",DBNdata,
#'dbnpar = list(samestruct=TRUE, slices=5, stationary=TRUE),
#'DBN=TRUE,bgnodes=c(1,2,3))
#'\dontrun{
#'orderDBNfit<-iterativeMCMC(dbnscore,chainout = TRUE, mergetype = "skeleton",scoreout=TRUE,alpha=0.4)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="trans",n.dynamic=12,n.static=3)
#'plotdiffs.DBN(orderDBNfit$max$DAG,DBNmat,struct="init",n.dynamic=12,n.static=3)
#'}
#'@export
plotdiffs.DBN<-function(eDBN,trueDBN,struct=c("init","trans"),n.dynamic,n.static=0) {
  
  if(!is.matrix(eDBN)) {
    adj<-graph2m(eDBN)
  } else {
    adj<-eDBN
  }
  if(!is.matrix(trueDBN)) {
    adjt<-graph2m(trueDBN)
  } else {
    adjt<-trueDBN
  }
  
  nodelabs<-colnames(adj)
  statcol<-"lightgrey"
  dyn1col<-"#f7f4f9"
  dyn2col<-"#d4b9da"
  
  if(struct=="init") {
    n<-n.static+n.dynamic
    nodelabs<-nodelabs[1:(n.dynamic+n.static)]
    
    if(n.static>0){
      legadj<-matrix(0,nrow=2,ncol=2)
      colnames(legadj)<-c("static","t=1")
      legadj[1,2]<-1
      legendG<-m2graph(legadj)
      staticnames<-nodelabs[1:n.static]
      legcol<-c(statcol,dyn1col)
      names(legcol)<-c("static","t=1")
    }
    
    dynamicnames<-nodelabs[1:n.dynamic+n.static]
    
    adj<-adj[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    adjt<-adjt[1:(n.dynamic+n.static),1:(n.dynamic+n.static)]
    jointmat<-1*(adj|adjt)
    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL
    
    diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
    for(i in 1:n) {
      for(j in 1:n) {
        bi1<-adj[i,j]+adj[j,i]
        bi2<-adjt[i,j]+adjt[j,i]
        if(bi1==2 | bi2==2) {
          if(bi1!=bi2) {
            if(bi2==2 & bi1==0) { #FN
              diffmat[i,j]<-3
              diffmat[j,i]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
            } else if (bi2==2 & bi1==1) { #ED FN
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
              MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
            } else if (bi2==1 & bi1==2) { #ED FP
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
              
            } else { #FP
              diffmat[i,j]<-2
              diffmat[j,i]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
              FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
            }
          }
        } else {
          if(adj[i,j]!=adjt[i,j]){
            if(adj[j,i]!=adjt[j,i]) {#ED
              if(adj[i,j]==1) {
                diffmat[i,j]<-4 
                jointmat[j,i]<-0
                EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
              } else {
                diffmat[j,i]<-4
                jointmat[i,j]<-0
                EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
              }
            } else if (adj[i,j]==1) { #FP
              diffmat[i,j]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            } else {#FN
              diffmat[i,j]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              
            }
          }
        }
      }
    }
    
    
    jointarcs<-adjacency2edgel(jointmat,nodes=nodelabs)
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                    edgemode = 'directed')
    
    if(n.static!=0) {  
      sg1 = subGraph(dynamicnames, graph.obj)
      sg2 = subGraph(staticnames, graph.obj)
      sgL = list(list(graph=sg1, cluster = TRUE),
                 list(graph=sg2, cluster = TRUE))
      colvector<-c(rep(statcol,n.static),
                   rep(dyn1col,n.dynamic))
      names(colvector)<-nodelabs
    } else {
      colvector<-c(rep(dyn1col,n.dynamic))
      names(colvector)<-nodelabs
    }
    
    graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL)
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle")
    if(!is.null(FPlist)) {
      FP<-apply(FPlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
    }
    if(!is.null(FNlist)) {
      FN<-apply(FNlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
    }
    if(!is.null(EDlist)) {
      ED<-apply(EDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
    }
    if(!is.null(BiFP)) {
      BiFP<-apply(BiFP, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
    }
    if(!is.null(BiFN)) {
      BiFN<-apply(BiFN, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
      
    }
    if(!is.null(MissDlist)) {
      MissD<-apply(MissDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
      
    }
    if(!is.null(ExtrDlist)) {
      ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
      
    }
    
    graph.par(list(graph=list(main="Comparison of DBNs initial structures")))
    
    if(n.static>0) {
      layout(matrix(c(1,1,1,1,3,
                      1,1,1,1,2,
                      1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
      Rgraphviz::renderGraph(graph.plot)
      graph::plot(legendG,attrs=list(graph=list(rankdir="TB"),
                                     edge=list(color="white")),nodeAttrs=list(fillcolor=legcol),
                  main="node codes")
    } else {
      Rgraphviz::renderGraph(graph.plot)
    }
  }
  
  if(struct=="trans") {
    
    n<-n.static+2*n.dynamic
    
    if(n.static>0) {
      legadj<-matrix(0,nrow=3,ncol=3)
      legadj[1,2]<-1
      legadj[2,3]<-1
      colnames(legadj)<-c("static","t=i","t=i+1")
      legendG<-m2graph(legadj)
      legcol<-c(statcol,dyn1col,dyn2col)
      names(legcol)<-c("static","t=i","t=i+1")
      
      colvector = c(rep(statcol,n.static), rep(dyn1col,n.dynamic),rep(dyn2col,n.dynamic))
      names(colvector)<-nodelabs
    } else {
      legadj<-matrix(0,nrow=2,ncol=2)
      legadj[1,2]<-1
      colnames(legadj)<-c("t=i","t=i+1")
      legendG<-m2graph(legadj)
      legcol<-c(dyn1col,dyn2col)
      names(legcol)<-c("t=i","t=i+1")
      colvector = c(rep(dyn1col,n.dynamic),rep(dyn2col,n.dynamic))
      names(colvector)<-nodelabs
    }
    
    adj<-DBNcut(adj[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    adjt<-DBNcut(adjt[1:(n.static+2*n.dynamic),1:(n.static+2*n.dynamic)],n.dynamic,n.static)
    jointmat<-1*(adj|adjt)
    
    arcslist<-adjacency2edgel(jointmat,nodes=nodelabs)
    
    graph.obj = new("graphNEL", nodes = nodelabs, edgeL = arcslist,
                    edgemode = 'directed')
    
    staticnames<-nodelabs[1:n.static]
    dyn1names<-nodelabs[1:n.dynamic+n.static]
    dyn2names<-nodelabs[1:n.dynamic+n.static+n.dynamic]
    
    sgStat = subGraph(staticnames, graph.obj)
    sgDyn1 = subGraph(dyn1names, graph.obj)
    sgDyn2 = subGraph(dyn2names, graph.obj)
    
    sgL = list(list(graph=sgStat, cluster = TRUE, attrs = c(rankdir="TB",rank="sink")),
               list(graph=sgDyn1, cluster = TRUE, attrs = c(rank="same")),
               list(graph=sgDyn2, cluster = TRUE, attrs = c(rank="same")))
    
    
    graph.plot = Rgraphviz::layoutGraph(graph.obj, subGList = sgL)
    graph::nodeRenderInfo(graph.plot)<-list(fill=colvector,shape="circle")
    
    
    #define edges with wrong directions
    FPlist<-NULL
    FNlist<-NULL
    EDlist<-NULL
    MissDlist<-NULL #missing directions
    ExtrDlist<-NULL #extra directions
    
    BiFP<-NULL
    BiFN<-NULL
    
    
    diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
    for(i in 1:n) {
      for(j in 1:n) {
        bi1<-adj[i,j]+adj[j,i]
        bi2<-adjt[i,j]+adjt[j,i]
        if(bi1==2 | bi2==2) {
          if(bi1!=bi2) {
            if(bi2==2 & bi1==0) { #FN
              diffmat[i,j]<-3
              diffmat[j,i]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
            } else if (bi2==2 & bi1==1) { #ED FN
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
              MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
            } else if (bi2==1 & bi1==2) { #ED FP
              diffmat[i,j]<-4
              diffmat[j,i]<-4
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
              ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
              
            } else { #FP
              diffmat[i,j]<-2
              diffmat[j,i]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
              FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
            }
          }
        } else {
          if(adj[i,j]!=adjt[i,j]){
            if(adj[j,i]!=adjt[j,i]) {#ED
              if(adj[i,j]==1) {
                diffmat[i,j]<-4 
                jointmat[j,i]<-0
                EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
              } else {
                diffmat[j,i]<-4
                jointmat[i,j]<-0
                EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
              }
            } else if (adj[i,j]==1) { #FP
              diffmat[i,j]<-2
              FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            } else {#FN
              diffmat[i,j]<-3
              FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
              
            }
          }
        }
      }
    }
    
    
    if(!is.null(FPlist)) {
      FP<-apply(FPlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
    }
    if(!is.null(FNlist)) {
      FN<-apply(FNlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
    }
    if(!is.null(EDlist)) {
      ED<-apply(EDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
    }
    if(!is.null(BiFP)) {
      BiFP<-apply(BiFP, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
    }
    if(!is.null(BiFN)) {
      BiFN<-apply(BiFN, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
      graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
      
    }
    if(!is.null(MissDlist)) {
      MissD<-apply(MissDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
      
    }
    if(!is.null(ExtrDlist)) {
      ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
      graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
      graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
      
    }
    
    graph.par(list(graph=list(main="Comparison of DBNs transition structures")))
    layout(matrix(c(1,1,1,1,3,
                    1,1,1,1,2,
                    1,1,1,1,3), nrow = 3, ncol = 5, byrow = TRUE))
    
    Rgraphviz::renderGraph(graph.plot)
    graph::plot(legendG,attrs=list(graph=list(rankdir="TB"),edge=list(color="white")),
                nodeAttrs=list(fillcolor=legcol),main="node codes")
  }
  
}

#' Plotting difference between two graphs
#' 
#' This function plots an estimated graph such that the edges which are different to the ground truth graph are highlighted. 
#' 
#'@param edag object of class graphNEL (or its adjacency matrix), representing estimated structure (not necessarily acyclic) to be compared to the ground truth graph
#'@param truedag object of class graphNEL (or its adjacency matrix), representing the ground truth structure (not necessarily acyclic)
#'@param clusters (optional) a list of nodes to be represented on the graph as clusters 
#'@return plots the graph which includes edges from edag and truedag, however edges which are different in edag compared to truedag are coloured according to the type of a mistake: false positive with red, false negative with dashed grey, error in direction with magenta
#'@examples
#'Asiascore<-scoreparameters(8,"bde",Asia)
#'Asiamap<-orderMCMC(Asiascore)
#'plotdiffs(Asiamap$max$DAG,Asiamat)
#'Asiacp<-pcalg::dag2cpdag(m2graph(Asiamat))
#'mapcp<-pcalg::dag2cpdag(m2graph(Asiamap$max$DAG))
#'plotdiffs(mapcp,Asiacp)
#'@export
plotdiffs<-function(edag,truedag,clusters=NULL) {
  
  if(!is.matrix(edag)) {
    adj<-graph2m(edag)
  } else {
    adj<-edag
  }
  if(!is.matrix(truedag)) {
    adjt<-graph2m(truedag)
  } else {
    adjt<-truedag
  }
  
  nodelabs<-colnames(adj)
  if(is.null(nodelabs)) {
    nodelabs<-paste("v",1:ncol(adj),sep="")
  }
  jointmat<-1*(adj|adjt)
  n<-nrow(adj)
  
  FPlist<-NULL
  FNlist<-NULL
  EDlist<-NULL
  MissDlist<-NULL #missing directions
  ExtrDlist<-NULL #extra directions
  
  BiFP<-NULL
  BiFN<-NULL
  
  
  diffmat<-matrix(0,nrow=nrow(jointmat),ncol=ncol(jointmat))
  for(i in 1:n) {
    for(j in 1:n) {
      bi1<-adj[i,j]+adj[j,i]
      bi2<-adjt[i,j]+adjt[j,i]
      if(bi1==2 | bi2==2) {
        if(bi1!=bi2) {
          if(bi2==2 & bi1==0) { #FN
            diffmat[i,j]<-3
            diffmat[j,i]<-3
            FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
            FNlist<-rbind(FNlist,c(nodelabs[j],nodelabs[i]))
          } else if (bi2==2 & bi1==1) { #ED FN
            diffmat[i,j]<-4
            diffmat[j,i]<-4
            MissDlist<-rbind(MissDlist,c(nodelabs[i],nodelabs[j]))
            MissDlist<-rbind(MissDlist,c(nodelabs[j],nodelabs[j]))
          } else if (bi2==1 & bi1==2) { #ED FP
            diffmat[i,j]<-4
            diffmat[j,i]<-4
            ExtrDlist<-rbind(ExtrDlist,c(nodelabs[i],nodelabs[j]))
            ExtrDlist<-rbind(ExtrDlist,c(nodelabs[j],nodelabs[j]))
            
          } else { #FP
            diffmat[i,j]<-2
            diffmat[j,i]<-2
            FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
            FPlist<-rbind(FPlist,c(nodelabs[j],nodelabs[i]))
          }
        }
      } else {
        if(adj[i,j]!=adjt[i,j]){
          if(adj[j,i]!=adjt[j,i]) {#ED
            if(adj[i,j]==1) {
              diffmat[i,j]<-4 
              jointmat[j,i]<-0
              EDlist<-rbind(EDlist,c(nodelabs[i],nodelabs[j]))
            } else {
              diffmat[j,i]<-4
              jointmat[i,j]<-0
              EDlist<-rbind(EDlist,c(nodelabs[j],nodelabs[i]))
            }
          } else if (adj[i,j]==1) { #FP
            diffmat[i,j]<-2
            FPlist<-rbind(FPlist,c(nodelabs[i],nodelabs[j]))
          } else {#FN
            diffmat[i,j]<-3
            FNlist<-rbind(FNlist,c(nodelabs[i],nodelabs[j]))
            
          }
        }
      }
    }
  }
  
  
  jointarcs<-adjacency2edgel(jointmat,nodes=nodelabs)
  graph.obj = new("graphNEL", nodes = nodelabs, edgeL = jointarcs,
                  edgemode = 'directed')
  if(is.null(clusters)){
    subGList = NULL } else {
      numsubg<-length(clusters)
      sg<-list()
      subGList<-list()
      for(i in 1:numsubg) {
        sg[[i]] <-subGraph(clusters[[i]],  graph.obj)
        subGList[[i]]<-list(graph = sg[[i]], cluster = TRUE)
      }
    }
  graph.plot = Rgraphviz::layoutGraph(graph.obj,subGList = subGList)
  if(!is.null(FPlist)) {
    FP<-apply(FPlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][FP] = "red"
  }
  if(!is.null(FNlist)) {
    FN<-apply(FNlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][FN] = "grey"
    graph::edgeRenderInfo(graph.plot)[["lty"]][FN] = "dashed"
  }
  if(!is.null(EDlist)) {
    ED<-apply(EDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][ED] = "blue"
  }
  if(!is.null(BiFP)) {
    BiFP<-apply(BiFP, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][BiFP] = "red"
  }
  if(!is.null(BiFN)) {
    BiFN<-apply(BiFN, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][BiFN] = "grey"
    graph::edgeRenderInfo(graph.plot)[["lty"]][BiFN] = "dashed"
    
  }
  if(!is.null(MissDlist)) {
    MissD<-apply(MissDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][MissD] = "blue"
    graph::edgeRenderInfo(graph.plot)[["lty"]][MissD] = "solid"
    
  }
  if(!is.null(ExtrDlist)) {
    ExtrD<-apply(ExtrDlist, 1, paste, collapse = "~")
    graph::edgeRenderInfo(graph.plot)[["col"]][ExtrD] = "blue"
    graph::edgeRenderInfo(graph.plot)[["lty"]][ExtrD] = "solid"
    
  }
  
  graph::edgeRenderInfo(graph.plot)[["lwd"]] = 2
  
  layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))
  Rgraphviz::renderGraph(graph.plot)

  par(mar = c(0,0,0,0))
  plot(1:10,1:10,type="n", axes = FALSE, xlab = "", ylab = "")
  op <- par(cex = 1.7)
  legend("center",legend=c("true positive",
                           "false positive",
                           "false negative",
                           "error in direction"),lty=c(1,1,2,1), lwd=c(2,2,2,2),
         col=c("black","red","grey","blue"),bty="n")
}









