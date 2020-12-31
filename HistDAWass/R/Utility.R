# res=Prepare(x,simplify,qua,standardize)
Prepare =function (x, simplify, qua, standardize){
 resu=c_Prepare(x,simplify,qua,standardize)
  return(resu)
}

# compute fast SSQ
ComputeFastSSQ=function(subMM){
  resu=c_ComputeFASTSSQ(subMM)
  return(resu)
}

# compute fast DIST
ComputeFast_L2_SQ_WASS_D=function(subMM){
  Dist=c_MM_L2_SQ_WASS_D(subMM)
  return(Dist)
}


WH.ADPT.KMEANS.TOTALSSQ=function(x,memb,m,lambdas,proto){
  vars=ncol(x@M)
  ind=nrow(x@M)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  #Compute general weights for the global PROTOTYPE
  mu=apply(memb^m,2,sum)
  protoGEN=x[1,]
  for (variables in (1:vars)){
    protoGEN@M[1,variables][[1]]=WH.vec.mean(proto[,variables],mu)
  }
  #Compute general weights for the global PROTOTYPE2
  Mat_of_means=matrix(0,ind,vars)
  CentredD=x
  for (i in 1:ind){
    for (j in 1:vars){
      Mat_of_means[i,j]=x@M[i,j][[1]]@m
      CentredD@M[i,j][[1]]=x@M[i,j][[1]]-x@M[i,j][[1]]@m
    }
  }
  W_centers=matrix(0,ind,vars)
  W_dist_cent=matrix(0,ind,vars)
  for (i in 1:ind){
    for (j in 1:vars){
      W_centers[i,j]=lambdas[(j*2-1),which.max(memb[i,])]
      W_dist_cent[i,j]=lambdas[(j*2),which.max(memb[i,])]
    }
  }
  protoGEN2=x[1,]
  mprot=matrix(0,1,vars)
  for (variables in (1:vars)){
    mprot[1,variables]=sum(Mat_of_means[,variables]*W_centers[,variables])/sum(W_centers[,variables])
    protoGEN2@M[1,variables][[1]]=WH.vec.mean(CentredD[,variables],W_dist_cent[,variables])+mprot[1,variables]
  }
  #Compute BETWEEN
  BSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (variables in (1:vars)){
    GENPROTm=protoGEN2@M[1,variables][[1]]@m
    GENPROTcent=protoGEN2@M[1,variables][[1]]-GENPROTm
    for (clu in 1:k){
      LOCPROTm=proto@M[clu,variables][[1]]@m
      LOCPROTcent=proto@M[clu,variables][[1]]-LOCPROTm
      BSQ_clu[1,variables,clu]=sum(memb[,clu])*lambdas[(variables*2-1),clu]*(GENPROTm-LOCPROTm)^2
      BSQ_clu[2,variables,clu]=sum(memb[,clu])*lambdas[(variables*2),clu]*WassSqDistH(GENPROTcent,LOCPROTcent)
    }
  }
  #Compute the total fuzzy sum of SQUARES
  TSQ_clu=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu[1,variables,cluster]=TSQ_clu[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu[2,variables,cluster]=TSQ_clu[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
      }
    }
  }
  TSQ_clu2=array(0,dim = c(2,vars,k), dimnames=list(c('Mean','Variability'),colnames(x@M),rownames(proto@M)))
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        tmpD=WassSqDistH(x@M[indiv,variables][[1]],protoGEN2@M[1,variables][[1]],details=T)
        tmpD_mean=as.numeric(tmpD[2])
        tmpD_centered=as.numeric(tmpD[1]-tmpD[2])
        TSQ_clu2[1,variables,cluster]=TSQ_clu2[1,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2-1),cluster]*tmpD_mean
        TSQ_clu2[2,variables,cluster]=TSQ_clu2[2,variables,cluster]+((memb[indiv,cluster])^m)*lambdas[(variables*2),cluster]*tmpD_centered
      }
    }
  }
  return(RES=list(protoGEN=proto,TSQdetailed=TSQ_clu, TSQ=sum(TSQ_clu),TSQdetailed2=TSQ_clu2, TSQ2=sum(TSQ_clu2),
                  BSQdetailed=BSQ_clu, BSQ=sum(BSQ_clu)))
}

ComputeFastSSQ_Fuzzy=function(subMM,memb,m){
  ind=ncol(subMM)-1
  rr=nrow(subMM)
  SSQ=0
  W=memb^m
  W2=W/sum(W)
  m1=subMM[,1:ind]%*%as.matrix(W2)
  
  m1=as.vector(m1)
  p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
  cm=(m1[2:rr]+m1[1:(rr-1)])/2
  rm=(m1[2:rr]-m1[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
    ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    SSQ=SSQ+(sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2)))*W[indiv]
  }
  return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
}

ComputeFast_Fuzzy_TOT_SSQ=function(MM,proto,memb,m){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
#computeGenProt
mus=apply(memb^m,2,sum)
mus=mus/sum(mus)
GP=new('MatH',1,var)
for (variable in 1:var){
  dG=MM[[variable]][,1]*0
  for (clu in 1:k){
    dG=dG+proto@M[clu,variable][[1]]@x*mus[clu]
  }
  tmp=new('distributionH',x=dG,p=MM[[variable]][,(ind+1)])
  GP@M[1,variable][[1]]=tmp
}
#compute TOTALSSQ
SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
for (variable in 1:var){
  rr=nrow(MM[[variable]])
  p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
  cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
  rm=(GP@M[1,variable][[1]]@x[2:rr]-GP@M[1,variable][[1]]@x[1:(rr-1)])/2
  for (indiv in 1:ind){
    ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
    ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
    dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
    dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
    for (clu in 1:k){
      SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+memb[indiv,clu]^m*dc
      SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+memb[indiv,clu]^m*dv
    }
  }
}
return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}

ComputeFast_Fuzzy_Adaptive_TOT_SSQ=function(MM,proto,memb,m,lambdas,theta){
  ind=ncol(MM[[1]])-1
  var=length(MM)
  k=ncol(memb)
  #computeGenProt
  
  GP=new('MatH',1,var)
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    ci=(MM[[variable]][2:rr,1:ind]+MM[[variable]][1:(rr-1),1:ind])/2
    ri=(MM[[variable]][2:rr,1:ind]-MM[[variable]][1:(rr-1),1:ind])/2
    wi=(MM[[variable]][2:rr,(ind+1)]-MM[[variable]][1:(rr-1),(ind+1)])
    LMc=memb^m*matrix(lambdas[(variable*2-1),]^theta,ind,k, byrow=TRUE)
    LMc=LMc/sum(LMc)
    LMv=memb^m*matrix(lambdas[(variable*2),]^theta,ind,k, byrow=TRUE)
    LMv=LMv/sum(LMv)
    Cen=MM[[variable]][,1]*0
    MG=0
    McG=ci[,1]*0
    RcG=ri[,1]*0
    for (clu in 1:k){
      for (i in 1:ind){
        mui=sum(ci[,i]*wi)
        MG=MG+LMc[i,clu]*mui
        Cen=Cen+(MM[[variable]][,i]-mui)*LMv[i,clu]
      }
    }
    Cen=Cen+MG
    tmp=new('distributionH',x=Cen,p=MM[[variable]][,(ind+1)])
    GP@M[1,variable][[1]]=tmp
  }
  #compute G by protos
  
  
  #compute TOTALSSQ
  SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
  for (variable in 1:var){
    rr=nrow(MM[[variable]])
    p1=MM[[variable]][2:rr,ind+1]-MM[[variable]][1:(rr-1),ind+1]
    cm=(GP@M[1,variable][[1]]@x[2:rr]+GP@M[1,variable][[1]]@x[1:(rr-1)])/2
    rm=(GP@M[1,variable][[1]]@p[2:rr]-GP@M[1,variable][[1]]@p[1:(rr-1)])/2
    for (indiv in 1:ind){
      ci=(MM[[variable]][2:rr,indiv]+MM[[variable]][1:(rr-1),indiv])/2
      ri=(MM[[variable]][2:rr,indiv]-MM[[variable]][1:(rr-1),indiv])/2
      dist=sum(p1*((ci-cm)^2+1/3*(ri-rm)^2))
      dc=(sum(p1*ci)-GP@M[1,variable][[1]]@m)^2
      dv=dist-dc
      for (clu in 1:k){
        SSQ_det[1,variable,clu]=SSQ_det[1,variable,clu]+
          lambdas[(variable*2-1),clu]^theta*(memb[indiv,clu]^m)*dc
        SSQ_det[2,variable,clu]=SSQ_det[2,variable,clu]+
          lambdas[(variable*2),clu]^theta*(memb[indiv,clu]^m)*dv
      }
    }
  }
  return(list(SSQ=sum(SSQ_det),SSQ_det=SSQ_det,ProtoGEN=GP))
}

#'From real data to distributionH.
#' 
#' @param data a set of numeric values.
#' @param algo (optional) a string. Default is "histogram", i.e. the function "histogram"
#' defined in the \code{\link[histogram]{histogram}}  package. \cr If "base" 
#' the \code{\link[graphics]{hist}} function is used. \cr
#' "FixedQuantiles" computes the histogram using as breaks a fixed number of quantiles.\cr
#' "ManualBreaks" computes a histogram where braks are provided as a vector of values.\cr
#' "PolyLine" computes a histogram using a piecewise linear approximation of the empirical
#' cumulative distribution function using the "Ramer-Douglas-Peucker algorithm", 
#'  \url{http://en.wikipedia.org/wiki/Ramer-Douglas-Peucker_algorithm}. 
#'  An \code{epsilon} parameter is required.
#'  The data are scaled in order to have a standard deviation equal to one.
#' @param type (optional) a string. Default is "combined" and generates 
#' a histogram having regularly spaced breaks (i.e., equi-width bins) and 
#' irregularly spaced ones. The choice is done accordingly with the penalization method described in 
#' \code{\link[histogram]{histogram}}. "regular" returns equi-width binned histograms, "irregular" returns
#' a histogram without equi-width histograms. 
#' @param qua a positive integer to provide if \code{algo="FixedQuantiles"} is chosen. Default=10.
#' @param breaks a vector of values to provide if  \code{algo="ManualBreaks"} is chosen.
#' @param epsilon a number between 0 and 1 to provide if \code{algo="PolyLine"} is chosen. Default=0.01.
#' @return A \code{distributionH} object, i.e. a distribution.
#' @importFrom histogram histogram
#' @importFrom stats quantile sd
#' @export
#' @examples
#' data=rnorm(n = 1000,mean = 2,sd = 3)
#' mydist=data2hist(data)
#' plot(mydist)
#' @seealso \code{\link[histogram]{histogram}} function
data2hist<-function(data, 
                    algo="histogram",
                    type="combined",
                    qua=10,
                    breaks=numeric(0),
                    epsilon=0.01){
  a=switch(algo,"base"=1, "histogram"=2, "FixedQuantiles"=3,"ManualBreaks"=4 ,
           "PolyLine"=5,2)
  t=switch(type,"regular"=1, "irregular"=2, "combined"=3, 3)
  if (a==1){
    h <- hist(data)
    x=h$breaks
    counts=h$counts
    counts[which(counts==0)]=0.001
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==2){
   
    if (t==1){
      h<-histogram(data,type="regular", verbose=FALSE,plot=FALSE)
    }
    if (t==2){
      h<-histogram(data,type="irregular",verbose=FALSE,plot=FALSE)
    }
    if (t==3){
      h<-histogram(data,type="combined",verbose=FALSE,plot=FALSE)
    }
    x=h$breaks
    counts=h$counts
    counts[which(counts==0)]=0.001
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==3){
    p=c(0:qua)/qua
    x=rep(0,qua+1)
    for (i in 0:qua){      
      x[i+1]=as.numeric(quantile(data,probs = p[i+1]))
      if (i>0){
        if(x[i+1]<=x[i]) x[i+1]=x[i]+(1e-10)
      }
    }    
  }
  if (a==4){
    #checking limits
    
    ini=min(data)
    end=max(data)
    rr=end-ini
    tol=min(1e-10,(rr*(1e-10)))
    end=end+tol
    breaks=sort(unique(c(breaks,ini,end)))
    breaks=breaks[which(breaks>=ini)]
    breaks=breaks[which(breaks<=end)]
    #end check
    data.cut = cut(data, breaks, right=FALSE)
    counts=as.vector(table(data.cut))
    x=breaks
    p=c(0,cumsum(counts/sum(counts)))
  }
  if (a==5){
    data=sort(data)
    stdev=sd(data)
    
    cums=c(0:(length(data)-1))/(length(data)-1)
    points=cbind(data/stdev,cums)
    resu=DouglasPeucker(as.matrix(points),epsilon=epsilon)
    x=as.vector(resu[,1])*stdev
    p=as.vector(resu[,2])
  }
  p[1]=0
  p[length(p)]=1
  mydist<- distributionH(x,p)
  return(mydist)
}
#'Ramer-Douglas-Peucker algorithm for curve fitting with a PolyLine
#' 
#' @param points a 2D matrix with the coordinates of 2D points
#' @param epsilon an number between 0 and 1. Recomended 0.01.
#' @return A matrix with the points of segments of a Poly Line.
#' @export
#' @seealso \code{\link[HistDAWass]{data2hist}} function
#' @export
DouglasPeucker=function(points,epsilon){
  dmax=0
  index=0
  end=nrow(points)
  ResultList=numeric(0)
  if (end<3) return (ResultList=rbind(ResultList,points))
  for (i in 2:(end-1)){
    d=ShortestDistance(points[i,], line=rbind(points[1,],points[end,]))
    if (d>dmax){
      index=i
      dmax=d
    }
  }
  #if dmax is greater than epsilon recursively apply
  if (dmax>epsilon){
   # print(dmax)
    recResults1=DouglasPeucker(points[1:index,],epsilon)
    recResults2=DouglasPeucker(points[index:end,],epsilon)
    ResultList=rbind(ResultList,recResults1,recResults2)
   
  }
  else
  {
    ResultList=rbind(ResultList,points[1,],points[end,])
  }
  ResultList=as.matrix(ResultList[!duplicated(ResultList),])
  colnames(ResultList)=c("x","p")
  return(ResultList)
}
#' Shortes distance from a point o a 2d segment
#' 
#' @param p coordinates of a point
#' @param line a 2x2 matrix with the coordinates of two points defining a line
#' @return A numeric value, the Euclidean distance of point \code{p} to the \code{line}.
#' @export
#' @seealso \code{\link[HistDAWass]{data2hist}} function and \code{\link[HistDAWass]{DouglasPeucker}} function
#' @export
ShortestDistance=function(p, line){
  x1=line[1,1]
  y1=line[1,2]
  x2=line[2,1]
  y2=line[2,2]
  x0=p[1]
  y0=p[2]
  d=abs((y2-y1)*x0-(x2-x1)*y0+x2*y1-y2*x1)/sqrt((y2-y1)^2+(x2-x1)^2)
  return(as.numeric(d))
}

##Adaptive distance computation
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

###############################################################
ComputeFast_L2_SQ_WASS_DMAT=function(MM,prot,DETA=FALSE){
  ind=ncol(MM[[1]])-1
  rr=nrow(MM[[1]])
  Dist=matrix(0,length(MM),ind)
  DM=Dist
  DV=Dist
  for (v in 1:length(MM)){
    rr=nrow(MM[[v]])
    MPro=t(t(prot@M[1,v][[1]]@x))%*%matrix(1,1,ind)
    p1=MM[[v]][2:rr,ind+1]-MM[[v]][1:(rr-1),ind+1]
    
    c1=(MM[[v]][2:rr,1:ind]+MM[[v]][1:(rr-1),1:ind])/2
    r1=(MM[[v]][2:rr,1:ind]-MM[[v]][1:(rr-1),1:ind])/2
    c2=(MPro[2:rr,]+MPro[1:(rr-1),])/2
    r2=(MPro[2:rr,]-MPro[1:(rr-1),])/2
    
    Dist[v,]=p1%*%(((c1-c2)^2)+ 1/3*((r1-r2)^2))
    if (DETA==TRUE){
      DM[v,]=(p1%*%(c1-c2))^2
      DV[v,]=Dist[v,]-DM[v,]
    }
  }
  if (DETA==FALSE){
    return(t(Dist))}
  else{
    return(list(Dist=t(Dist), DM=t(DM),DV=t(DV)))
  }
}

########################
# compute fast DIST
ComputeFast_L2_SQ_WASS_D=function(subMM){
  Dist=c_MM_L2_SQ_WASS_D(subMM)
  return(Dist)
}

##  Weights computation
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
            if ((num==0)&(denom==0)){browser()}
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
            if ((num==0)&(denom==0)){browser()}
            lambdas[(variables*2-1),cluster]=num/denom
            num=(prod(apply(distances[,,2],MARGIN=c(1),sum)))^(1/vars)
            denom=max(sum(distances[variables,,2]),1e-10)
            if ((num==0)&(denom==0)){browser()}
            lambdas[(variables*2),cluster]=num/denom}
          else{
            lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
            lambdas[(variables*2),cluster]=lambdas[(variables*2),1]
          }
        }
        if (schema==3){#a weigth for one variable and for each cluster
          num=(prod(distances[,cluster,1]))^(1/vars)
          denom=max(distances[variables,cluster,1],1e-10)
          if ((num==0)&(denom==0)){browser()}
          lambdas[(variables*2-1),cluster]=num/denom
          lambdas[(variables*2),cluster]=num/denom
        }
        if (schema==4){#two weigths for the two components for each variable and each cluster
          num=(prod(distances[,cluster,1]))^(1/vars)
          denom=max(distances[variables,cluster,1],1e-10)
          if ((num==0)&(denom==0)){browser()}
          lambdas[(variables*2-1),cluster]=num/denom
          num=(prod(distances[,cluster,2]))^(1/vars)
          denom=max(distances[variables,cluster,2],1e-10)
          if ((num==0)&(denom==0)){browser()}
          lambdas[(variables*2),cluster]=num/denom
        }
        if (schema==5){#two weigths for the two components for each variable 
          #multiplied all equal one
          if (cluster<=1){
            num1=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/(2*vars))
            denom1=max(sum(distances[variables,,1]),1e-10)
            if ((num1==0)&(denom1==0)){browser()}
            num2=(prod(apply(distances[,,2],MARGIN=c(1),sum)))^(1/(2*vars))
            denom2=max(sum(distances[variables,,2]),1e-10)
            if ((num2==0)&(denom2==0)){browser()}
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
          if ((num==0)&(denom==0)){browser()}
          lambdas[(variables*2-1),cluster]=num/denom
          #num=(prod(distances[,cluster,2]))^(1/vars)
          denom=max(distances[variables,cluster,2],1e-10)
          if ((num==0)&(denom==0)){browser()}
          lambdas[(variables*2),cluster]=num/denom
          
        }
      }else{
        #sum ugual to 1
        if (schema==1){#one weigth for one variable 1W GS Global-sum
          if (cluster<=1){
            num=rep(sum(distances[variables,,1]),vars)
            den=apply(distances[,,1],1,sum)
            if ((num==0)&(den==0)){browser()}
            lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
            lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
          }else{lambdas[(variables*2-1),cluster]=lambdas[(variables*2-1),1]
          lambdas[(variables*2),cluster]=lambdas[(variables*2),1]}
        }
        if (schema==2){#two weigths for the two components for each variable 2W GS Global-sum
          if (cluster<=1){
            num=rep(sum(distances[variables,,1]),vars)
            den=apply(distances[,,1],1,sum)
            if ((num==0)&(den==0)){browser()}
            lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
            num=rep(sum(distances[variables,,2]),vars)
            den=apply(distances[,,2],1,sum)
            if ((num==0)&(den==0)){browser()}
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
          if ((num==0)&(den==0)){browser()}
          lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
        }
        if (schema==4){#two weigths for the two components for each variable and each cluster 2W LS Local-sum
          num=rep(distances[variables,cluster,1],vars)
          den=distances[,cluster,1]
          if ((num==0)&(den==0)){browser()}
          lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          
          num=rep(distances[variables,cluster,2],vars)
          den=distances[,cluster,2]
          if ((num==0)&(den==0)){browser()}
          lambdas[(variables*2),cluster]=1/sum((num/den)^(1/(theta-1)))
        }
        if (schema==5){#two weigths for the two components for each variable 2W GS Global-sum
          if (cluster<=1){
            unonum=sum(distances[variables,,1])
            unoden=apply(distances[,,1],MARGIN=c(1),sum)
           # if ((unonum==0)&(unoden==0)){browser()}
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
WH.ADPT.FCMEANS.SSQ_FAST=function(MM,x,memb,m,lambdas,proto,theta){
  ind=ncol(MM[[1]])-1
  vars=length(MM)
  k=ncol(memb)
  lambdas[is.na(lambdas)]=0
  
  DM=array(0,c(ind,vars,k))
  DV=DM
  for (cl in 1:k){
    tmp=ComputeFast_L2_SQ_WASS_DMAT(MM,proto[cl,],DETA=T)
    DM[,,cl]=tmp$DM
    DV[,,cl]=tmp$DV
  }
  
  SSQ=0
  for (indiv in 1:ind){
    for (cluster in 1:k){
      for (variables in (1:vars)){
        
        #   tmpD=WassSqDistH(x@M[indiv,variables][[1]],proto@M[cluster,variables][[1]],details=T)
        tmpD_mean=DM[indiv,variables,cluster]
        tmpD_centered=DV[indiv,variables,cluster]
        SSQ=SSQ+((memb[indiv,cluster])^m)*(lambdas[(variables*2-1),cluster]^theta*tmpD_mean
                                           +lambdas[(variables*2),cluster]^theta*tmpD_centered)
        
      }
    }
  }
  return(SSQ)
}

#compute a vector of quantiles faster----
compQ_vect=function(object,vp){
  if(is.unsorted(vp)) vp=sort(vp)
  x=object@x
  pro=object@p
  q=COMP_Q_VECT(object@x, object@p, vp) 
  return(q)
}


MEAN_VA=function(MAT,wei){
  
res=MEDIA_V(MAT, wei)
return(res)
}

# an alternative to table
Table2=function(x,N, full=TRUE){
  howm=rep(0,N)
  names(howm)=paste0(c(1:N))
  for (i in 1:N){
    howm[i]=length(which(x==i))
  }
  if (full){
    return(howm)}else{
      return(howm[which(howm>0)])
    }
}