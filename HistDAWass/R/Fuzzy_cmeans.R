
#' Fuzzy c-means of a dataset of histogram-valued data
#' @description The function implements the fuzzy c-means for a set of histogram-valued data. 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param m A number grater than 0, a fuzziness coefficient (default \code{m}=1.6).
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.   
#' @return a list with the results of the fuzzy c-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' @slot solution A list.Returns the best solution among the \code{rep}etitions, i.e. 
#' the one having the minimum sum of squares deviation.
#' @slot solution$membership A matrix. The membership degree of each unit to each cluster.
#' @slot solution$IDX A vector. The crisp assignement to a cluster.
#' @slot solution$cardinality A vector. The cardinality of each final cluster (after the crisp assignement).
#' @slot solution$Crit A number. The criterion (Sum of square deviation 
#' from the prototypes) value at the end of the run.
#' @slot quality A number. The percentage of Sum of square deviation explained by the model. 
#' (The higher the better)
#' @examples
#' results=WH_fcmeans(x = BLOOD,k = 2,m = 1.5,rep = 10,simplify = TRUE,qua = 10,standardize = TRUE)
#' @export
WH_fcmeans =function (x,k,m=1.6, rep, 
                       simplify=FALSE,
                       qua=10,
                       standardize=FALSE){
  tol=1e-10
  # Check initial conditions and passed parameters
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  if (missing(rep)){
    rep=5
  }
  ind=nrow(x@M)
  vars=ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  
  # end of the preprocessing step
  
  ##------------------------------------------------FCMEANS-----------------------------------------------------
  solutions=list()
  criteria=numeric(0)
  for (repet in (1:rep)){
    #initialize matrix of memberships
    memb=matrix(runif(ind*k,0.01,1.01),ind,k)
    #MIGLIORA QUI
    for (i in (1:ind)){
      memb[i,]=memb[i,]/sum(memb[i,])
    }
    #
    cat(paste("FCM---------> rep  ",repet,"\n"))
    ## initialize clusters and prototypes
    proto=new("MatH",nrows=k,ncols=vars)
    
    ##  compute prototypes
    
    nr=get.MatH.rownames(proto)
    nc=get.MatH.varnames(proto)
    
    tmp3=c_ComputeFastSSQ_Fuzzy(MM,memb, m, nr, nc)
    
    SSQ=tmp3$SSQ;
    proto=tmp3$proto;
    
    GenCrit=sum(SSQ)
    OK=1 
    maxit=100
    #treshold=GenCrit*1e-10
    treshold=1e-3
    itcount=0
    
    while (OK==1){
      itcount=itcount+1
      #computing distances
      
      MD=c_ComputeFast_L2_SQ_WASS_DMAT_simple(MM,proto)
      #compute fuzzy membership
      #MIGLIORA QUI
      MD[which(MD<tol)]=tol
      for (i in 1:ind){
        for (cl in 1:k){
          memb[i,cl]=1/sum((MD[i,cl]/MD[i,])^(1/(m-1)))
        }
      }
      #recompute criterion
      OldCrit=GenCrit
      #MIGLIORA QUI
      for (clust in 1:k){
        for (variables in 1:vars){
          tmp=MM[[variables]][,1:ind]
          tmpP=rowSums(tmp*matrix(memb[,clust]^(m),nrow(tmp),ncol(tmp),byrow = TRUE))/sum(memb[,clust]^(m))
          proto@M[clust,variables][[1]]@x=as.vector(tmpP)
          proto@M[clust,variables][[1]]@p=as.vector(MM[[variables]][,ind+1])
          tmpvals=M_STD_H(proto@M[clust,variables][[1]])
          proto@M[clust,variables][[1]]@m=tmpvals[1]
          proto@M[clust,variables][[1]]@s=tmpvals[2]  
        }
      }
      
      tmpD=c_ComputeFast_L2_SQ_WASS_DMAT_simple(MM,proto)      
      SSQ=sum(tmpD*(memb^m))
      GenCrit=SSQ
      #cat(paste(itcount,GenCrit, "\n", sep="---->"))
      #check criterion
      #      if (abs(GenCrit-OldCrit)<treshold){OK=0}
      if (1-GenCrit/OldCrit<treshold){OK=0}
    }
    #crisp assignment
    IDX=apply(memb,1,which.max)
    cardinality=table(IDX)
    dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
    #    dimnames(cardinality)$IDX=paste("Cl",1:k,sep=".")
    colnames(proto@M)=colnames(x@M)
    rownames(proto@M)=paste("Clust", 1:k, sep = ". ")
    row.names(memb)=get.MatH.rownames(x)
    solutions=c(solutions,list(solution=list(membership=memb, IDX=IDX,cardinality=cardinality,proto=proto,Crit=GenCrit)))
    criteria=c(criteria,GenCrit)
    #plot(proto,type="DENS")
  }
  
  tmp=c_ComputeFast_Fuzzy_TOT_SSQ(MM,
                                  solutions[[which.min(criteria)]]$proto,
                                  solutions[[which.min(criteria)]]$membership,
                                  m)
  
  TOTSSQ=tmp$SSQ
  Crisp_clu=list()
  for (clu in 1:k){
    Crisp_clu[[clu]]=rownames(x@M)[solutions[[which.min(criteria)]]$IDX==clu]
  }
  return(best.solution=list(solution=solutions[[which.min(criteria)]],
                            rep=which.min(criteria),
                            SSQ_det=tmp$SSQ_det,
                            TSQ=tmp$SSQ,
                            WSQ=min(criteria),
                            BSQ=tmp$SSQ-min(criteria),
                            ProtoGEN=tmp$ProtoGEN,
                            quality=1-min(criteria)/TOTSSQ, Crisp_clu=Crisp_clu)
  )
}


#' @title Fuzzy c-means with adaptive distances for histogram-valued data
#' @description Fuzzy c-means of a dataset of histogram-valued data using different adaptive distances based on the L2 Wasserstein metric.
#' 
#' @param x A MatH object (a matrix of distributionH).
#' @param k An integer, the number of groups.
#' @param schema An integer. 1=one weight per variable, 2=two weights per variables (one for each component: the mean and the variability component),  
#' 3=one weight per variable and per cluster, 4= two weights per variable and per cluster.
#' @param m A number grater than 0, a fuzziness coefficient (default \code{m}=1.6).
#' @param rep An integer, maximum number of repetitions of the algorithm (default \code{rep}=5).
#' @param simplify A logic value (default is FALSE), if TRUE histograms are recomputed in order to speed-up the algorithm.
#' @param qua An integer, if \code{simplify}=TRUE is the number of quantiles used for recodify the histograms.
#' @param standardize A logic value (default is FALSE). If TRUE, histogram-valued data are standardized,  variable by variable, using the Wassertein based standard deviation. Use if one wants to have variables with std equal to one.
#' @param init.weights A string. (default='EQUAL'). EQUAL, all variables or components have the same weight; 'RANDOM', a random assignment is done.
#' @param weight.sys  A string. (default='PROD') PROD, Weights product is equal to one. SUM, the weights sum up to one.  
#' @param theta A number. (default=2) A parameter for the system of weights summing up to one. 
#' @param verbose A logic value (default is FALSE). If TRUE  some details are provided.   
#' @return The results of the fuzzy c-means of the set of Histogram-valued data \code{x} into  \code{k} cluster.
#' \item{solution}{A list.Returns the best solution among the \code{rep}etitions, i.e. the ona having the minimum sum of squares deviation.}
#' \item{solution$membership}{A matrix. The membership degree of each unit to each cluster.}
#' \item{solution$IDX}{A vector. The crisp assignement to a cluster.}
#' \item{solution$cardinality}{A vector. The cardinality of each final cluster (after the crisp assignement).}
#' \item{solution$Crit}{A number. The criterion (Sum od square deviation from the prototypes) value at the end of the run.}
#' \item{quality}{A number. The percentage of Sum of square deviation explained by the model. (The higher the better)}
#' @examples
#' results=WH_adaptive_fcmeans(x = BLOOD,k = 2,schema=4,m = 1.5,rep = 3,simplify = TRUE,
#'                            qua = 10,standardize = TRUE,init.weights='EQUAL', weight.sys='PROD')
#' @export
WH_adaptive_fcmeans=function (x,k=5,
                               schema, #1=one weight per variable  
                               #2=two weights per variables 
                               #3=one weight per variable and per cluster 
                               #4= two weights per variable and per cluster
                               #5= two weights per variable multiplied all
                               #6= two weights per variable and per cluster multiplied all
                               m=1.6,rep, 
                               simplify=FALSE,
                               qua=10,
                               standardize=FALSE, init.weights='EQUAL',weight.sys='PROD',theta=2,
                              verbose=FALSE){
  tol=1e-10
  if(weight.sys=='PROD') weight_sys=1
  if(weight.sys=='SUM') weight_sys=2
  # Check initial conditions and passed parameters
  if ((k<2) || (k>=nrow(x@M))){
    str=paste("The number of clusters is not appropriate:\n 2<= k<",
              nrow(x@M), "(no. of individuals in x)\n")
    stop(str)
  }
  if (missing(rep)){
    rep=5
  }
  ind=nrow(x@M)
  vars=ncol(x@M)
  ## we homogeneize data for speeding up the code if required
  tmp=Prepare(x,simplify,qua,standardize)
  MM=tmp$MM
  x=tmp$x
  # end of the preprocessing step
  
  TOTSSQ=0
  for (v in 1:vars){
    TOTSSQ=TOTSSQ+(WH.SSQ(x[,v]))
  }
  if (weight.sys=="PROD") {theta=1}
  ##------------------------------------------------ADAPTIVE FCMEANS2-----------------------------------------------------
  solutions=list()
  criteria=numeric(0)
  # proto_empty=new("MatH",nrows=k,ncols=vars)
  # rownames(proto_empty@M)=paste('Clust',c(1:k),sep="_")
  # colnames(proto_empty@M)=colnames(x@M)
  
  
  for (repet in (1:rep)){
    # select at random k prototypes
    #proto=proto_empty
    RSEL=sample(c(1:ind), size=k, replace = FALSE)
    proto=x[RSEL,]
    
    #initialize matrix of variables' weights
    if (init.weights=='EQUAL'){
      cat("Weights initialization === EQUAL  \n")
      if (weight.sys=='PROD'){
        lambdas=matrix(1,2*vars,k)
      }
      else{
        if (schema==5 | schema==6){
          lambdas=matrix(1/(2*vars),2*vars,k)
        }else{
          lambdas=matrix(1/vars,2*vars,k)
        }
      }
    }
    else{#Random initialization
      cat("Weights initialization === RANDOM  \n")
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
    cat(paste("AFCM--",weight.sys,"--Cl. ",k,"---------> rep  ",repet,"\n"))
    ## initialize memberships
    memb=matrix(1,nrow=ind,ncol=k)
    tmpR2=c_DISTA_ADA(x,MM,proto,vars,ind,k,memb,m,schema)
    diINDtoPROT_M=tmpR2$dista_to_pro_M
    diINDtoPROT_V=tmpR2$dista_to_pro_V
    # memberships computation (prototype, weights are fixed)
    
    if (weight.sys=='PROD'){
      memb=c_MEMB_comp(ind, vars, k, 
                       lambdas,
                       diINDtoPROT_M,
                       diINDtoPROT_V,
                       m,1)
      #TMP_SSQ=WH.ADPT.FCMEANS.SSQ_FAST(MM,x,memb,m,lambdas,proto,theta=1)
      #TMP_SSQ2=c_WH_ADPT_FCMEANS_SSQ(x,memb,m, lambdas, proto, 1)
      TMP_SSQ=c_WH_ADPT_FCMEANS_SSQ_FAST_NEW(diINDtoPROT_M,
                                             diINDtoPROT_V,
                                             memb,m,
                                             lambdas, 1)
    }else{
      memb=c_MEMB_comp(ind, vars, k, 
                       lambdas,
                       diINDtoPROT_M,
                       diINDtoPROT_V,
                       m,theta)
      #TMP_SSQ=WH.ADPT.FCMEANS.SSQ_FAST(MM,x,memb,m,lambdas,proto,theta=theta)
      #TMP_SSQ2=c_WH_ADPT_FCMEANS_SSQ(x,memb,m, lambdas, proto, theta)
      TMP_SSQ=c_WH_ADPT_FCMEANS_SSQ_FAST_NEW(diINDtoPROT_M,
                                             diINDtoPROT_V,
                                             memb,m,
                                             lambdas, theta)
    }
    
    ## initialize clusters and prototypes
    
    GenCrit=Inf
    OK=1 
    maxit=100
    
    treshold=1e-5
    itcount=0
    
    while (OK==1){
      OldCrit=GenCrit
      itcount=itcount+1
      ##S1)  representation step Stage 1) compute initial prototypes (fix membership and weights)
      for (cluster in 1:k){
        for (variables in (1:vars)){
          tmp=MM[[variables]][,1:ind]
          tmpP=rowSums(tmp*matrix(memb[,cluster]^(m),nrow(tmp),ncol(tmp),byrow = TRUE))/sum(memb[,cluster]^(m))
          
          proto@M[cluster,variables][[1]]@x=as.vector(tmpP)
          proto@M[cluster,variables][[1]]@p=as.vector(MM[[variables]][,ind+1])
          tmpvals=M_STD_H(proto@M[cluster,variables][[1]])
          proto@M[cluster,variables][[1]]@m=tmpvals[1]
          proto@M[cluster,variables][[1]]@s=tmpvals[2]
          
          
        }
      }
      ##S2)  representation step Stage 2) compute variable weigths (fix membership and prototypes)
      tmpR2=c_DISTA_ADA(x,MM,proto,vars,ind,k,memb,m,schema)
      diINDtoPROT_M=tmpR2$dista_to_pro_M
      diINDtoPROT_V=tmpR2$dista_to_pro_V
      distances_M=tmpR2$dista_M;
      distances_V=tmpR2$dista_V;
      
      # S2.2) Weights computation
      #NEW
      lambdas=c_ADA_F_WHEIGHT(list(distances_M,distances_V),k, vars, ind, schema, weight_sys,
                              theta, m)
      #  
      # S3) memberships computation (prototype, weights are fixed)
      if(weight.sys=='PROD'){
        #NEW
        memb=c_MEMB_comp(ind, vars, k, 
                         lambdas,
                         diINDtoPROT_M,
                         diINDtoPROT_V,
                         m,1)
        
        
        # TMP_SSQ=WH.ADPT.FCMEANS.SSQ_FAST(MM,x,memb,m,lambdas,proto,theta=1)
        #   TMP_SSQ2=c_WH_ADPT_FCMEANS_SSQ(x,memb,m, lambdas, proto, 1)
        TMP_SSQ=c_WH_ADPT_FCMEANS_SSQ_FAST_NEW(diINDtoPROT_M,
                                               diINDtoPROT_V,
                                               memb,m,
                                               lambdas, 1)
      }
      else{     
        
        #NEW
        memb=c_MEMB_comp(ind, vars, k, 
                         lambdas,
                         diINDtoPROT_M,
                         diINDtoPROT_V,
                         m,theta)
        #TMP_SSQ=WH.ADPT.FCMEANS.SSQ_FAST(MM,x,memb,m,lambdas,proto,theta=theta)
        #TMP_SSQ2=c_WH_ADPT_FCMEANS_SSQ_FAST(MM,x,memb,m, lambdas, proto, theta)
        TMP_SSQ=c_WH_ADPT_FCMEANS_SSQ_FAST_NEW(diINDtoPROT_M,
                                               diINDtoPROT_V,
                                               memb,m,
                                               lambdas, theta)
      }
      # browser()
      ## compute criterion (fix membership and weights)
      
      GenCrit=TMP_SSQ
      
     if (verbose) {cat(paste(itcount,GenCrit, "\n", sep="---->"))}
      #check criterion
      if (1-GenCrit/OldCrit<treshold){OK=0}
    }
    if (verbose){cat(paste(itcount,GenCrit, "\n", sep="---->"))}
    #crisp assignment
    IDX=apply(memb,1,which.max)
    cardinality=table(IDX)
    dimnames(cardinality)$IDX=paste("Cl",dimnames(cardinality)$IDX,sep=".")
    #    dimnames(cardinality)$IDX=paste("Clust",1:k,sep=".")
    colnames(proto@M)=colnames(x@M)
    rownames(proto@M)=paste("Clust", 1:k, sep = ". ")
    
    
    solutions=c(solutions,list(solution=list(var.weights=lambdas, membership=memb,
                                             IDX=IDX,cardinality=cardinality,proto=proto,Crit=GenCrit)))
    criteria=c(criteria,GenCrit)
    #plot(proto,type="DENS")
  }
  BestSOL=solutions[[which.min(criteria)]]
  if (weight.sys=='PROD'){
    resTOTSQ=ComputeFast_Fuzzy_Adaptive_TOT_SSQ(MM,proto=BestSOL$proto,
                                                memb=BestSOL$membership,m=m,
                                                lambdas=BestSOL$var.weights,theta=1)
  }
  else{
    resTOTSQ=ComputeFast_Fuzzy_Adaptive_TOT_SSQ(MM,proto=BestSOL$proto,
                                                memb=BestSOL$membership,m=m,
                                                lambdas=BestSOL$var.weights,theta=theta)
    
  }
  GEN_proto=resTOTSQ$ProtoGEN
  DET_TSQ=resTOTSQ$SSQ
  TOTFSSQ=sum(resTOTSQ$SSQ)
  quality=1-min(criteria)/TOTFSSQ
  
  return(best.solution=list(solution=BestSOL,
                            rep=which.min(criteria),
                            TSQ=TOTFSSQ,
                            WSQ=min(criteria),
                            BSQ=TOTFSSQ-min(criteria),
                            quality=1-min(criteria)/TOTFSSQ,
                            schema=schema, weight.sys=weight.sys))
}

WH.ADPT.FCMEANS.SSQ=function(x,memb,m,lambdas,proto,theta){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  SSQ=0
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        SSQ=SSQ+((memb[indiv,cluster])^m)*(lambdas[(variables*2-1),cluster]^theta*tmpD_mean
                                           +lambdas[(variables*2),cluster]^theta*tmpD_centered)
        
      }
    }
  }
  return(SSQ)
}



# S2.2) Weights computation
ADA_F_DIST=function(distances,k,vars,ind,schema,weight.sys,theta=2,m=1.6){
  lambdas=matrix(0,2*vars,k)
  for (variables in (1:vars)){
    for (cluster in 1:k){
      #product
      if (weight.sys=='PROD'){
        if (schema==1){#one weigth for one variable
          if (cluster<=1){
            num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
            denom=max(sum(distances[variables,,1]),1e-10)
            lambdas[(variables*2-1),cluster]=num/denom
            lambdas[(variables*2),cluster]=num/denom
          }else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==2){#two weigths for the two components for each variable
          if (cluster<=1){
            num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
            denom=max(sum(distances[variables,,1]),1e-10)
            lambdas[(variables*2-1),cluster]=num/denom
            num=(prod(apply(distances[,,2],MARGIN=c(1),sum)))^(1/vars)
            denom=max(sum(distances[variables,,2]),1e-10)
            lambdas[(variables*2),cluster]=num/denom}
          else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==3){#a weigth for one variable and for each cluster
          num=(prod(distances[,cluster,1]))^(1/vars)
          denom=max(distances[variables,cluster,1],1e-10)
          lambdas[(variables*2-1),cluster]=num/denom
          lambdas[(variables*2),cluster]=num/denom
        }
        if (schema==4){#two weigths for the two components for each variable and each cluster
          num=(prod(distances[,cluster,1]))^(1/vars)
          denom=max(distances[variables,cluster,1],1e-10)
          lambdas[(variables*2-1),cluster]=num/denom
          num=(prod(distances[,cluster,2]))^(1/vars)
          denom=max(distances[variables,cluster,2],1e-10)
          lambdas[(variables*2),cluster]=num/denom
        }
        if (schema==5){#two weigths for the two components for each variable 
          #multiplied all equal one
          if (cluster<=1){
            num1=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/(2*vars))
            denom1=max(sum(distances[variables,,1]),1e-10)
            
            num2=(prod(apply(distances[,,2],MARGIN=c(1),sum)))^(1/(2*vars))
            denom2=max(sum(distances[variables,,2]),1e-10)
            lambdas[(variables*2-1),cluster]=(num1*num2)/denom1
            lambdas[(variables*2),cluster]=(num1*num2)/denom2}
          else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==6){#two weigths for the two components for each variable and each cluster
          #multiplied all equal one
          num=(prod(distances[,cluster,1])*prod(distances[,cluster,2]))^(1/(2*vars))
          denom=max(distances[variables,cluster,1],1e-10)
          lambdas[(variables*2-1),cluster]=num/denom
          #num=(prod(distances[,cluster,2]))^(1/vars)
          denom=max(distances[variables,cluster,2],1e-10)
          lambdas[(variables*2),cluster]=num/denom
          
        }
      }else{
        #sum ugual to 1
        if (schema==1){#one weigth for one variable 1W GS Global-sum
          if (cluster<=1){
            num=rep(sum(distances[variables,,1]),vars)
            den=apply(distances[,,1],1,sum)
            lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
            lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
          }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
          lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
        }
        if (schema==2){#two weigths for the two components for each variable 2W GS Global-sum
          if (cluster<=1){
            num=rep(sum(distances[variables,,1]),vars)
            den=apply(distances[,,1],1,sum)
            lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
            num=rep(sum(distances[variables,,2]),vars)
            den=apply(distances[,,2],1,sum)
            lambdas[(variables*2),cluster]=1/sum((num/den)^(1/(theta-1)))
          }
          else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==3){#a weigth for one variable and for each cluster 1W LS Local-sum
          num=rep(distances[variables,cluster,1],vars)
          den=distances[,cluster,1]
          
          lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
        }
        if (schema==4){#two weigths for the two components for each variable and each cluster 2W LS Local-sum
          num=rep(distances[variables,cluster,1],vars)
          den=distances[,cluster,1]
          
          lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          
          num=rep(distances[variables,cluster,2],vars)
          den=distances[,cluster,2]
          lambdas[(variables*2),cluster]=1/sum((num/den)^(1/(theta-1)))
        }
        if (schema==5){#two weigths for the two components for each variable 2W GS Global-sum
          if (cluster<=1){
            unonum=sum(distances[variables,,1])
            unoden=apply(distances[,,1],MARGIN=c(1),sum)
            uno=sum((unonum/unoden)^(1/(theta-1)))
            dueden=apply(distances[,,2],MARGIN=c(1),sum)
            due=sum((unonum/dueden)^(1/(theta-1)))
            lambdas[(variables*2-1),cluster]=1/(uno+due)
            
            unonum=sum(distances[variables,,2])
            tre=sum((unonum/unoden)^(1/(theta-1)))
            quattro=sum((unonum/dueden)^(1/(theta-1)))
            lambdas[(variables*2),cluster]=1/(tre+quattro)
          }
          else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==6){#two weigths for the two components for each variable and each cluster 2W LS Local-sum
          unonum=distances[variables,cluster,1]
          unoden=distances[,cluster,1]
          uno=sum((unonum/unoden)^(1/(theta-1)))
          dueden=distances[,cluster,2]
          due=sum((unonum/dueden)^(1/(theta-1)))
          lambdas[(variables*2-1),cluster]=1/(uno+due)
          
          unonum=distances[variables,cluster,2]
          tre=sum((unonum/unoden)^(1/(theta-1)))
          quattro=sum((unonum/dueden)^(1/(theta-1)))
          lambdas[(variables*2),cluster]=1/(tre+quattro)
        }
      }
    }
  }
  return(lambdas=lambdas)
}

DISTA_ADA=function(x,MM,proto,vars,ind,k,memb,m,schema){ 
  distances=array(0,dim=c(vars,k,2))
  diINDtoPROT=array(0,dim=c(ind,vars,k,2))
  #####################################
  
  DM=array(0,c(ind,vars,k))
  DV=DM
  for (cl in 1:k){
    tmp=ComputeFast_L2_SQ_WASS_DMAT(MM,proto[cl,],DETA=T)
    DM[,,cl]=tmp$DM
    DV[,,cl]=tmp$DV
  }
  #  ComputeFast_L2_SQ_WASS_DMAT(MM,prot,DETA=T)
  
  for (cluster in 1:k){
    for (variables in (1:vars)){
      for (indiv in 1:ind){
        #         tmpD=ComputeFast_L2_SQ_WASS_D(cbind(MM[[variables]][,indiv],proto@M[cluster,variables][[1]]@x,
        #                                             MM[[variables]][,(ind+1)]))
        tmpD=(DM[indiv,variables,cluster]+DV[indiv,variables,cluster])
        
        if (schema==1){#one weigth for one variable
          distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*memb[indiv,cluster]^m
          diINDtoPROT[indiv,variables,cluster,1]=tmpD
        }
        if (schema==2|schema==5){#two weigths for the two components for each variable
          #tmpD_mean=(x@M[indiv,variables][[1]]@m-proto@M[cluster,variables][[1]]@m)^2
          tmpD_mean=DM[indiv,variables,cluster]
          #tmpD_centered=tmpD-tmpD_mean
          tmpD_centered=DV[indiv,variables,cluster]
          distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*memb[indiv,cluster]^m
          distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*memb[indiv,cluster]^m
          diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
          diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
        }
        if (schema==3){#a weigth for one variable and for each cluster
          distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD*memb[indiv,cluster]^m
          diINDtoPROT[indiv,variables,cluster,1]=tmpD
        }
        if (schema==4|schema==6){#two weigths for the two components for each variable and each cluster
          #tmpD_mean=(x@M[indiv,variables][[1]]@m-proto@M[cluster,variables][[1]]@m)^2
          tmpD_mean=DM[indiv,variables,cluster]
          #tmpD_centered=tmpD-tmpD_mean
          tmpD_centered=DV[indiv,variables,cluster]
          
          distances[variables,cluster,1]=distances[variables,cluster,1]+tmpD_mean*memb[indiv,cluster]^m
          distances[variables,cluster,2]=distances[variables,cluster,2]+tmpD_centered*memb[indiv,cluster]^m
          diINDtoPROT[indiv,variables,cluster,1]=tmpD_mean
          diINDtoPROT[indiv,variables,cluster,2]=tmpD_centered
        }
      }
    }
  }
  
  return(res=list(distances=distances,diINDtoPROT=diINDtoPROT))
}
MEMB_comp=function(ind,vars,k,lambdas,diINDtoPROT,m,theta){
  memb=matrix(0,ind,k)
  for (indiv in 1:ind){
    d1=rep(0,k)
    for (cluster in 1:k){
      for (variables in 1:vars){
        d1[cluster]=d1[cluster]+lambdas[(variables*2-1),cluster]^theta*diINDtoPROT[indiv,variables,cluster,1]+lambdas[(variables*2),cluster]^theta*diINDtoPROT[indiv,variables,cluster,2]
      }
    }
    for (cluster in 1:k){
      if (d1[cluster]>0){
        memb[indiv,cluster]=1/sum((d1[cluster]/d1)^(1/(m-1)))}
      else{memb[indiv,cluster]=1}
    }
  }
  return(memb)
}



