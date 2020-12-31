MDS = function(DataOrDistances,method='euclidean',OutputDimension=2,PlotIt=FALSE,Cls){
  # Classical multidimensional scaling of a data matrix. Also known as principal coordinates analysis
  # projection=MDS(Data) 
  # INPUT 
  # DataOrDistances[1:n,1:d]      array of data: n cases in rows, d variables in columns, matrix is not symmetric
  #                           or distance matrix, in this case matrix has to be symmetric
  # OPTIONAL
  # method                    method specified by distance string: 
  #                          'euclidean','cityblock=manhatten','cosine','chebychev','jaccard','minkowski','manhattan','binary' 
  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  #                           of the dimension chosen by cmdscale and OutputDimension
  #
  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  
  # OUTPUT is a list with following elements:
  # ProjectedPoints[1:n,OutputDimension]                    n by OutputDimension matrix containing coordinates of the Projection
  #
  # Eigenvalues[1:p]                                the eigenvalues of MDSvalues*MDSvalues'
  # Stress                                          Shephard-Kruskal Stress  
  # author: MT 06/2015
  if(missing(DataOrDistances))
		stop('No DataOrDistances given')
	DataOrDistances;
	if(!is.matrix(DataOrDistances))
		stop('DataOrDistances has to be a matrix, maybe use as.matrix()')
		
  if(isSymmetric(unname(DataOrDistances))){
    DataDists=DataOrDistances
    AnzVar=ncol(DataOrDistances)
    AnzData=nrow(DataOrDistances)
  }else{ #!isSymmetric
    AnzVar=ncol(DataOrDistances)
    AnzData=nrow(DataOrDistances)
    DataDists=as.matrix(dist(x=DataOrDistances,method = method))
  }# end if(isSymmetric(DataOrDistances))
  
  res=cmdscale(d=DataDists, k = OutputDimension, eig = TRUE, add = FALSE, x.ret = FALSE)
  ProjectedPoints=as.matrix(res$points)

  # calculation of Stress value
    outDistances <- as.matrix(dist(x = ProjectedPoints,method=method) )
  # Distances = Distances';
  # Shephard-Kruskal Stress 
  Stress = KruskalStress(DataDists, outDistances)
  if(PlotIt){
     if(missing(Cls)){
		AnzData=nrow(ProjectedPoints)
		Cls=rep(1,AnzData)
	}  
    
    string=paste0('MDS with Shephard-Kruskal Stress ',round(Stress,2))
    #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string)
      PlotProjectedPoints(ProjectedPoints,Cls,main=string)
  }         
  return(list(ProjectedPoints=ProjectedPoints,Eigenvalues=res$eig,Stress=Stress))

}