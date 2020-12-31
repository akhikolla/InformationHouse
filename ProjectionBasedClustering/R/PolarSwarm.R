PolarSwarm=function(DataOrDistances,method='euclidean',PlotIt=FALSE,Cls){
  
  if(missing(DataOrDistances))
    stop('No DataOrDistances given')
  DataOrDistances;
  if(!is.matrix(DataOrDistances))
    stop('DataOrDistances has to be a matrix, maybe use as.matrix()')
  
  
  if (!requireNamespace('DatabionicSwarm')) {
    message(
      'Subordinate projection package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Cls = rep(1, nrow(DataOrDistances)),
        Object = "Subordinate projection package is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }
  
  ModelObject=DatabionicSwarm::Pswarm(DataOrDistances,method = method)
  
  ProjectedPoints=ModelObject$ProjectedPoints
  if(PlotIt){
    if(missing(Cls)){
      AnzData=nrow(DataOrDistances)
      Cls=rep(1,AnzData)
    }  
    PlotProjectedPoints(ProjectedPoints,Cls,main="Polar Swarm Projection of the DBS algorithm")
  } 
  return(list(ProjectedPoints=ProjectedPoints,ModelObject=ModelObject))
}