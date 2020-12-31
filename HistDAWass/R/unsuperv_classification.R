# K-means of histogram data based on Wasserstein ----
#' K-means of a dataset of histogram-valued data
#' @description The function implements the k-means for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param verbose A logic value (default is FALSE). If TRUE, details on computations are shown.   
#' @return a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$centers A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references Irpino A., Verde R., Lechevallier Y. (2006). Dynamic clustering of histograms using Wasserstein 
#' metric. In: Rizzi A., Vichi M.. COMPSTAT 2006 - Advances in computational statistics. p. 869-876, 
#' Heidelberg:Physica-Verlag
#' @examples
#' results=WH_kmeans(x = BLOOD,k = 2, rep = 10,simplify = TRUE,
#' qua = 10,standardize = TRUE,verbose=TRUE)
#' @export
WH_kmeans =function (x,k, rep=5, 
                     simplify=FALSE,
                     qua=10,
                     standardize=FALSE,verbose=FALSE){
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  
  init="RPART"
  ind=nrow(x@M)
  vars=ncol(x@M)
  NAMER=get.MatH.rownames(x)
  NAMEV=get.MatH.varnames(x)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ## compute total sum of 
  TOTSSQ=0
  
  tmpd=distributionH()
  for (v in 1:vars){
#    tmp=ComputeFastSSQ(MM[[v]])
    tmp=c_ComputeFASTSSQ(MM[[v]])
    TOTSSQ=TOTSSQ+tmp$SSQ
    
  }
  solutions=list()
  criteria=numeric(0)
  for (repet in 1:rep){
    if(verbose){
    cat(paste("---------> rep  ",repet,"\n"))}
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    cards=rep(0,k)
    SSQ=matrix(0,k,vars)
    GenCrit=Inf;
    if (init=="RPART"){
      rperm=sample(1:ind)
      tmp=rep(c(1:k),ceiling(ind/k))[1:ind]
      IDX=tmp[rperm]
      ##  compute prototypes
      for (clu in 1:k){
        sel=(1:ind)[IDX==clu]
        cards=length(sel)
        for (j in 1:vars){
          #tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp=c_ComputeFASTSSQ(MM[[j]][,c(sel,(ind+1))])
          proto@M[clu,j][[1]]@x=tmp$mx
          proto@M[clu,j][[1]]@p=tmp$mp
          tt=M_STD_H(proto@M[clu,j][[1]])
          proto@M[clu,j][[1]]@m=tt[1]
          proto@M[clu,j][[1]]@s=tt[2]
          # 
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[clu,j][[1]]=tmpH
          SSQ[clu,j]=tmp$SSQ
        }
      }
    }
    GenCrit=sum(SSQ)
    OK=1 
    maxit=100
    treshold=GenCrit*1e-10
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #assign
      tmp=c_DISTA_M(MM,proto)
      
      MD=tmp$dist
      IDX=tmp$wm
      ##########################
      # for (i in 1:ind){
      #     IDX[i]=which.min(MD[i,])
      # }
      # browser()
      ################################
      #recompute
      OldCrit=GenCrit
      
      #check for empty clusters
      moved=numeric(0)
      for (i in 1:k){
        sel=(1:ind)[IDX==i]
        #show(sel)
        if (length(sel)==0){
          cat("empty cluster\n")
          which.max(apply(MD,1,max))
          tomove=which.max(apply(MD,1,max))
          IDX[tomove]=i
          moved=c(moved,tomove)
          MD[tomove,]=rep(0,k)
          #show(MD)
        }
      }
      for (clu in 1:k){
        sel=(1:ind)[IDX==clu]
        
        cards=length(sel)
        for (j in 1:vars){
          #tmp=ComputeFastSSQ(MM[[j]][,c(sel,(ind+1))])
          tmp=c_ComputeFASTSSQ(MM[[j]][,c(sel,(ind+1))])
          proto@M[clu,j][[1]]@x=tmp$mx
          proto@M[clu,j][[1]]@p=tmp$mp
          tt=M_STD_H(proto@M[clu,j][[1]])
          proto@M[clu,j][[1]]@m=tt[1]
          proto@M[clu,j][[1]]@s=tt[2]
          # tmpH=new('distributionH',x=tmp$mx, p=tmp$mp)
          # proto@M[i,j][[1]]=tmpH
          SSQ[clu,j]=tmp$SSQ
        }
      }
      GenCrit=sum(SSQ)
      if (verbose){
      cat(paste(itcount,GenCrit, "\n", sep="---->"))}
      #check criterion
      if (abs(GenCrit-OldCrit)<treshold){
        OK=0
      }
    }
    cardinality=table(IDX)
    names(IDX)=NAMER
    dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
    row.names(proto@M)=paste("Cl",c(1:k),sep=".")
    colnames(proto@M)=NAMEV
    solutions=c(solutions,list(solution=list(IDX=IDX,cardinality=cardinality,
                                             centers=proto,Crit=GenCrit)))
    criteria=c(criteria,GenCrit)
  }
  #plot(solutions[[which.min(criteria)]]$centers,type="DENS")
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            quality=1-min(criteria)/TOTSSQ))
}


# Hierarchical clustering of histogram data based on Wasserstein ----
#' Hierarchical clustering of histogram data
#' @description The function implements a Hierarchical clustering 
#'  for a set of histogram-valued data, based on the L2 Wassertein distance.
#'  Extends the \code{hclust} function of the \pkg{stat} package. 
#' @param x A MatH object (a matrix of distributionH).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param distance A string default "WDIST" the L2 Wasserstein distance (other distances will be implemented)
#' @param method A string, default="complete", is the the agglomeration method to be used.
#'  This should be (an unambiguous abbreviation of) one of "\code{ward.D}", "\code{ward.D2}",
#'   "\code{single}", "\code{complete}", "\code{average}" (= UPGMA), "\code{mcquitty}" 
#'   (= WPGMA), "\code{median}" (= WPGMC) or "\code{centroid}" (= UPGMC).
#' @seealso \code{\link{hclust}} of \pkg{stat} package for further details.
#' @return An object of class hclust which describes the tree produced by
#'  the clustering process. 
#' @references Irpino A., Verde R. (2006). A new Wasserstein based distance for the hierarchical clustering 
#' of histogram symbolic data. In: Batanjeli et al. Data Science and Classification, IFCS 2006. p. 185-192,
#'  BERLIN:Springer, ISBN: 3-540-34415-2
#' @examples
#' results=WH_hclust(x = BLOOD,simplify = TRUE, method="complete")
#' plot(results) # it plots the dendrogram
#' cutree(results,k = 5) # it returns the labels for 5 clusters
#' @importFrom stats hclust quantile as.dist
#' @export
WH_hclust =function (x, 
                     simplify=FALSE,
                     qua=10,
                     standardize=FALSE,
                     distance="WDIST",
                     method="complete"){
  ind=nrow(x@M)
  vars=ncol(x@M)
  
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ##compute distance matrix
  if (distance=="WDIST"){
    d=sqrt(c_Fast_D_Mat(MM))
    }
  
  rownames(d)=rownames(x@M)
  colnames(d)=rownames(x@M)
  hc=hclust(as.dist(d),method=method)
  return(hc)
}


# Adaptive-distance based k-means ----
#' K-means of a dataset of histogram-valued data using adaptive  Wasserstein distances
#' @description The function implements the k-means using adaptive distance for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param schema a number from 1 to 4 \cr
#' 1=A weight for each variable (default) \cr
#' 2=A weight for the average and the dispersion component of each variable\cr
#' 3=Same as 1 but a different set of weights for each cluster\cr
#' 4=Same as 2 but a different set of weights for each cluster 
#' @param init (optional, do not use) initialization for partitioning the data default is 'RPART', other strategies shoul be implemented.
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, 
#' using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @param weight.sys a string. Weights may add to one ('SUM') or their product is equal to 1 ('PROD', default).
#' @param theta a number. A parameter if \code{weight.sys='SUM'}, default is 2.  
#' @param init.weights a string how to initialize weights: 'EQUAL' (default), all weights are the same, 
#' 'RANDOM', weights are initalised at random.
#' @param verbose A logic value (default is FALSE). If TRUE, details on computations are shown.   
#' 
#' @return a list with the results of the k-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares criterion.
#' @slot solution$IDX A vector. The clusters at which the objects are assigned.
#' @slot solution$cardinality A vector. The cardinality of each final cluster.
#' @slot solution$centers A \code{MatH} object with the description of centers.
#' @slot solution$Crit A number. The criterion (Sum od square deviation 
#' from the centers) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @references Irpino A., Rosanna V., De Carvalho F.A.T. (2014). Dynamic clustering of histogram data 
#' based on adaptive squared Wasserstein distances. EXPERT SYSTEMS WITH APPLICATIONS, vol. 41, p. 3351-3366, 
#' ISSN: 0957-4174, doi: http://dx.doi.org/10.1016/j.eswa.2013.12.001
#' @examples
#' results=WH_adaptive.kmeans(x = BLOOD,k = 2, rep = 10,simplify = TRUE,qua = 10,standardize = TRUE)
#' @importFrom stats runif
#' @export

WH_adaptive.kmeans=function (x,k,
                              schema=1, #1=VariableGLOBAL 2=componentGLOBAL 3=Variable x Cluster 4=Components x cluster 
                              init, rep,  
                              simplify=FALSE,
                              qua=10,
                              #empty.action="singleton",
                              standardize=FALSE,
                              weight.sys='PROD',
                              theta=2,
                              init.weights='EQUAL',
                              verbose=FALSE){
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  if (missing(init)){
    init="RPART"
  }
  if (missing(rep)){
    rep=5
  }
  ind=nrow(x@M)
  vars=ncol(x@M)
  if (weight.sys=='PROD'){theta=1}
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  
  ## compute total sum of 
  TOTSSQ=0
  for (v in 1:vars){
    #tmp=ComputeFastSSQ(MM[[v]])
    tmp=c_ComputeFASTSSQ(MM[[v]])
    TOTSSQ=TOTSSQ+tmp$SSQ
  }
  
  solutions=list()
  criteria=numeric(0)
  repet=0
  while(repet<rep){
    repet=repet+1
    if(verbose) cat(paste("---------> rep  ",repet,"\n"))
    #initialize matrix of variables' weights
    if (init.weights=='EQUAL'){
      if (verbose) cat("Weights initialization === EQUAL  \n")
      if (weight.sys=='PROD'){
        lambdas=matrix(1,2*vars,k)}
      else{lambdas=matrix(1/vars,2*vars,k)}
    }
    else{#Random initialization
      if(verbose) cat("Weights initialization === RANDOM  \n")
      m1=matrix(runif((vars*k),0.01,0.99),vars,k)
      m2=matrix(runif((vars*k),0.01,0.99),vars,k)
      m1=m1/matrix(rep(apply(m1,2,sum)),vars,k,byrow = TRUE)
      m2=m2/matrix(rep(apply(m2,2,sum)),vars,k,byrow = TRUE)
      if (weight.sys=='PROD'){
        m1=exp(m1*vars-1)
        m2=exp(m2*vars-1)
        
      }
      
      if (schema==1){
        m1=matrix(m1[,1],nrow = vars,ncol = k)
        m2=m1
        
      }
      if (schema==2){
        m1=matrix(rep(m1[,1],k),vars,k)
        m2=matrix(rep(m2[,1],k),vars,k)
      }
      if (schema==3){m2=m1}
      
      lambdas=matrix(0,2*vars,k)
      colnames(lambdas)=paste('Clust',c(1:k),sep="_")
      
      n1=paste('M_Var.',c(1:vars),sep="_")
      n2=paste('C_Var.',c(1:vars),sep="_")
      nr=list()
      for (nn in 1:vars){
        nr[[nn*2-1]]=n1[[nn]]
        nr[[nn*2]]=n2[[nn]]
      }
      rownames(lambdas)=nr
      lambdas[(c(1:vars)*2-1),]=m1
      lambdas[(c(1:vars)*2),]=m2
    }
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    cards=rep(0,k)
    SSQ=matrix(0,k,vars)
    GenCrit=Inf;
    if (init=="RPART"){
      rperm=sample(1:ind)
      tmp=rep(c(1:k),ceiling(ind/k))[1:ind]
      IDX=tmp[rperm]
      
    }
    memb=matrix(0,ind,k)
    for(indiv in 1:ind){memb[indiv,IDX[indiv]]=1}
    GenCrit=Inf
    OK=1 
    maxit=100
    treshold=1e-5
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #STEP 1 ##  compute prototypes
      for (cluster in 1:k){
        sel=(1:ind)[IDX==cluster]
        
        cards=length(sel)
        for (variables in 1:vars){
          tmp=c_ComputeFASTSSQ(MM[[variables]][,c(sel,(ind+1))])
          proto@M[cluster,variables][[1]]@x=tmp$mx
          proto@M[cluster,variables][[1]]@p=tmp$mp
          tt=M_STD_H(proto@M[cluster,variables][[1]])
          proto@M[cluster,variables][[1]]@m=tt[1]
          proto@M[cluster,variables][[1]]@s=tt[2]
          SSQ[cluster,variables]=tmp$SSQ
        }
      }
      #       TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      #       cat('\n after prototypes SSQ:', TMP_SSQ,'\n')
      #STEP 2 ##  compute weights (Lambda) Fixed Prototypes and partition
      #######################################################################################
      
      TMPRESU=c_STEP_2_ADA_KMEANS(MM,proto, schema,memb)
      distances=array(0,dim=c(vars,k,2))
      distances[,,1]=TMPRESU$distances1
      distances[,,2]=TMPRESU$distances2
      diINDtoPROT=array(0,dim=c(ind,vars,k,2))
      for (cluster in 1:k){
        diINDtoPROT[,,cluster,1]=TMPRESU$dIpro_m[[cluster]]
        diINDtoPROT[,,cluster,2]=TMPRESU$dIpro_v[[cluster]]
      }
      
      
      
      ##  compute weights
      # S2.2) Weights computation
      if (weight.sys=='PROD'){PROSUM=1}else{PROSUM=2}
      lambdas=c_STEP_2_2_WEIGHTS_ADA_KMEANS(TMPRESU$distances1,
                                            TMPRESU$distances2,
                                            PROSUM, schema, theta)
      
      ############################################################################
      #       TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto)
      #       cat('\n after weights SSQ:', TMP_SSQ,'\n')      
      #recompute
      OldCrit=GenCrit
      # S3) affectation (prototype, weights are fixed)
      resu= c_STEP_3_AFFECT_ADA_KMEANS(lambdas, 
                                       TMPRESU$dIpro_m,TMPRESU$dIpro_v, ind, k,  vars)
      IDX=resu$IDX;
      
      #check for empty clusters and restart
      empty=FALSE;
      for (i in 1:k){
        if (sum(IDX==i)==0){
          
          cat("empty cluster\n")
          empty=TRUE
          OK=0
        }
      }
      if (empty){repet=repet-1}
      memb=matrix(0,ind,k)
      for(indiv in 1:ind){memb[indiv,IDX[indiv]]=1}
      #browser()
      TMP_SSQ=resu$SSQ
      #TMP_SSQ=WH.ADPT.FCMEANS.SSQ_FAST(MM,x,memb,1,lambdas,proto,theta=theta)
      #browser()
      #TMP_SSQ=WH.ADPT.FCMEANS.SSQ(x,memb,1,lambdas,proto,theta)
      GenCrit=TMP_SSQ
      if(is.na(GenCrit)){
        cat('isNAN')
      }
      
      if(verbose)   cat(paste(itcount,GenCrit, "\n", sep="---->"))
      #check criterion
      if (abs(GenCrit-OldCrit)<treshold){OK=0}
    }
    if(!empty){  
      cardinality=table(IDX)
      dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
      solutions=c(solutions,list(solution=list(IDX=IDX,cardinality=cardinality,proto=proto, weights=lambdas,Crit=GenCrit)))
      criteria=c(criteria,GenCrit)
    }
  }
  best.solution=c(solutions[[which.min(criteria)]])
  #plot(best.solution$proto,type="DENS")
  memb=matrix(0,ind,k)
  for (i in 1:ind){
    memb[i,best.solution$IDX[i]]=1;
  }
  
  resTOTSQ=c_WH_ADPT_KMEANS_TOTALSSQ(x,memb,m=1,
                                     lambdas=best.solution$weights,proto=best.solution$proto)
  # resTOTSQ=c_WH_ADPT_KMEANS_TOTALSSQ(x,memb,m=1,
  #                                    lambdas=best.solution$weights^theta,proto=best.solution$proto)
  GEN_proto=resTOTSQ$gen
  DET_TSQ=list(TSQ_M=resTOTSQ$TSQ_m,TSQ_V=resTOTSQ$TSQ_v)
  TOTFSSQ=sum(resTOTSQ$TSQ_m+resTOTSQ$TSQ_v)
  DET_TSQ2=list(TSQ2_M=resTOTSQ$TSQ2_m,TSQ2_V=resTOTSQ$TSQ2_v)
  TOTFSSQ2=sum(resTOTSQ$TSQ2_m+resTOTSQ$TSQ2_v)
  BSQ=sum(resTOTSQ$BSQ_m+resTOTSQ$BSQ_v)
  
  quality1=1-min(criteria)/TOTFSSQ
  quality2=1-min(criteria)/TOTFSSQ2
  best.solution=c(best.solution,TOTSSQ=TOTFSSQ2,BSQ=BSQ,WSQ=best.solution$Crit,quality=BSQ/TOTFSSQ2)
  return(best.solution)
  #  return(0)
}

