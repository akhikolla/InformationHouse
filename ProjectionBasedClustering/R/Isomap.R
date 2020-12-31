Isomap  = function(Distances,k,OutputDimension=2,PlotIt=FALSE,Cls){
  # Isomap procetion as introduced in 2000 by Tenenbaum, de Silva and Langford
  # projection=Isomap(Data,k=10)

  # INPUT 
  # Distances[1:d,1:d]      DistanceMatrix
  # OPTIONAL

  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  #                           of the dimension chosen by cmdscale and OutputDimension
  # k                       number of k nearest neighbors              

  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  
  # OUTPUT is a list with following elements:
  # ProjectedPoints[1:n,OutputDimension]   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
  #Note if Data fragmented choose an higher k
  
  # author: MT 06/2015
	if(missing(Distances))
		stop('No Distances given')
  Distances;
	if(!is.matrix(Distances))
		stop('Distances has to be a matrix, maybe use as.matrix()')
	
  if (!requireNamespace('vegan')) {
    message(
      'Subordinate projection package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Cls = rep(1, nrow(Distances)),
        Object = "Subordinate projection package is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
  	
  if(missing(k)) stop('k nearest neighbor value missing')
 
    res=vegan::isomap(Distances,ndim = 2,k=40,fragmentedOK=T, path = "shortest")
  # requireRpackage("RDRToolbox",biocite=T)
  # 
  # 
  # ProjectedPoints=matrix(unlist(Isomap(data=Data,dim=OutputDimension,k=k,verbose=T)),ncol=OutputDimension)
    ProjectedPoints=res$points
if(PlotIt){
  if(missing(Cls)){
		AnzData=nrow(Distances)
		Cls=rep(1,AnzData)
	}  
  
  string=paste0('Isomap projection with k ',k)
  #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string)
    PlotProjectedPoints(ProjectedPoints,Cls,main=string)
}                    
return(list(ProjectedPoints=ProjectedPoints))

}