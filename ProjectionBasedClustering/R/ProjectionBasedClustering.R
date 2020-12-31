ProjectionBasedClustering=function(k,DataOrDistances,BestMatches,LC,StructureType=TRUE,PlotIt=FALSE,method='euclidean'){
#Cls=DBSclustering(k,DataOrDistances,BestMatches,LC,StructureType=TRUE,PlotIt=F,method='euclidean')
# automated Clustering approach of the DataBionicSwarm with abstact U distances
# INPUT
# k                   number of classes, how many to you see in the 3d landscape?
# DataOrDistances                Matrix of DataOrDistances that will be used. One DataPoint per row
# BestMatches         Array with positions of Bestmatches=ProjectedPoints
# LC
# OPTIONAL  
# StructureType           compact structure of clusters assumed, =FALSE: connected structure of clusters assumed
# PlotIt              Plots Dendrogramm
# method              do not change 
# OUTPUT 
# Cls                 vector with selected classes of the bestmatches
#  
# author: MT 07/17
#  
#  NOTE: das ist eine eigene Idee: Nehme die Distantz in LoetschUltsch2014 definiert und stecke sie in
#		Distanz=ShortestGraphPaths -> Cluster sind immer in sich geschlossen in 2D
#							-> Politische Karte irrelevant -> Nahteil: Falls Projektion Fehler hat (Punke des einen Clusters innerhalb des anderen Clusters)
#							-> Hat Clusterung fehler und CLusterung ist nichtmehr dichte basiert


  if(missing(DataOrDistances))
    stop('Input is missing')
	
  if(!is.matrix(DataOrDistances)){
    warning('Input of DataOrDistances or Distances is not a matrix. Trying to transform it to a matrix')
    DataOrDistances=as.matrix(DataOrDistances)
  }
  
  if (!isSymmetric(unname(DataOrDistances)))
    string='Data'
  else
    string='Distances'
  

  if(sum(!is.finite(DataOrDistances))!=0){
    warning(paste0('Some entries in ',string,' are not finite or missing values. Computations may not work'))
  }
  if(!is.numeric(DataOrDistances)){
    warning(paste0(string,' is not numeric. Trying to transform it to a numeric type.'))
    mode(DataOrDistances)="numeric"
  }
 
  if(isSymmetric(unname(DataOrDistances))){
    InputD=DataOrDistances
  }else{
    if (requireNamespace('parallelDist',quietly = TRUE)) {
      InputD=as.matrix(parallelDist::parallelDist(DataOrDistances,method=method))
    }else{
      InputD=as.matrix(dist(DataOrDistances,method = method))
    }
  }
  
    GOutput=Delaunay4Points(BestMatches, Grid = LC, IsToroid=T,PlotIt=F)
    Dist=ShortestGraphPathsC(GOutput,InputD)
  if(StructureType){
    pDist=as.dist(Dist)
    hc <- hclust(pDist,method="ward.D")
    m="Compact Projection-based Clustering"
  }else{
    ind=which(GOutput==0,arr.ind=T)
    Dist2=Dist*GOutput
    Dist2[ind]=max(Dist)*2
    pDist=as.dist(Dist2)
    hc <- hclust(pDist,method="single")
    m="Connected Projection-based Clustering"
  }
    
  Cls=cutree(hc,k) 

  
  counter = 0
  bool = T
  NumberOfClassesSet = k
  while (counter < 0.05 * length(Cls)) {
    counter = counter + 1
    uniqueClasses <- sort(na.last = T, unique(Cls))
    numberOfClasses <- length(uniqueClasses)
    countPerClass <- rep(0, numberOfClasses)
    for (i in 1:numberOfClasses) {
      inClassI <-
        sum(Cls == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in cls
      countPerClass[i] = inClassI
    }
    
    if (max(countPerClass) > 0.95 * length(Cls)) {
      k = k + 1
      Cls = cutree(hc, k)
      uniqueClasses <- sort(na.last = T, unique(Cls))
      numberOfClasses <- length(uniqueClasses)
      countPerClass <- rep(0, numberOfClasses)
      for (i in 1:numberOfClasses) {
        inClassI <-
          sum(Cls == uniqueClasses[i]) # counts all occurances of uniqueClass[i] in cls
        countPerClass[i] = inClassI
      }
      ind = order(countPerClass, decreasing = F, na.last = T)
      n = numberOfClasses - NumberOfClassesSet
      for (l in 1:n) {
        indc = which(Cls == uniqueClasses[ind[l]])
        Cls[indc] = 99
      }
    } else{
      bool = FALSE
    }
    if (!bool)
      break
    # print(counter)
  }
  
  if(PlotIt){
    x=as.dendrogram(hc)
    if(requireNamespace('dendextend')){
      #what is the ordering of the cluster in dendrogram
      # from left to right
      Cls_tmp=Cls[order.dendrogram(x)]
      #count frequency in that ordering
      uniqueClasses <- unique(Cls_tmp)
      numberOfClasses <- length(uniqueClasses)
      #countPerClass <- rep(0, numberOfClasses)
      countPerClass=list()
      for (i in uniqueClasses) {
        inClassI <- sum(Cls_tmp == uniqueClasses[i])
        countPerClass[[i]] = inClassI
      }
      names(countPerClass)=uniqueClasses
      countPerClass=unlist(countPerClass)
      #get the right number of colors
      cols=ProjectionBasedClustering::DefaultColorSequence[1:numberOfClasses]
     #what would be the ordering of datra based on frequency
      data_order=order(countPerClass,decreasing = TRUE) #from highest frequency
      #what would be the orders of the branches
      unique_reordered=uniqueClasses[data_order]
      # fit that order to the colors
      cols_order = match(table = unique_reordered,uniqueClasses)
      cols=cols[cols_order]
      #branch colors with specific set of colors based on cluster frequency
      x=dendextend::set(x,"branches_k_color", k = k,cols)
    }
    plot(x, main=m,xlab="No. of Data Points N", ylab="Ultrametric Portion of Distance",sub=" ",leaflab ="none")
    axis(1,col="black",las=1)
  }

  return(Cls)
}