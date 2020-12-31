SammonsMapping = function(DataOrDistances,method='euclidean',OutputDimension=2,PlotIt=FALSE,Cls){
  # verbesserte MDS durch Normieren durch die Dimension des Eingaberaums
  # projection=SammonsMapping(Data)
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
  # ProjectedPoints[1:n,OutputDimension]                   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
  #
  # Stress 	                                       The final stress achieved.
  # author: MT 06/2015
  
  if (!requireNamespace('MASS')) {
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
  
  if (missing(DataOrDistances))
    stop('No DataOrDistances given')
  DataOrDistances
  
  if (!is.matrix(DataOrDistances))
    stop('DataOrDistances has to be a matrix, maybe use as.matrix()')
  
  
  if (isSymmetric(DataOrDistances)) {
    DataDists = DataOrDistances
    AnzVar = ncol(DataOrDistances)
    AnzData = nrow(DataOrDistances)
  } else{
    #!isSymmetric
    AnzVar = ncol(DataOrDistances)
    AnzData = nrow(DataOrDistances)
    DataDists = as.matrix(dist(x = DataOrDistances, method = method))
  }# end if(isSymmetric(DataOrDistances))
  
  res = MASS::sammon(
    d = DataDists,
    y = cmdscale(d = DataDists, k = OutputDimension),
    k = OutputDimension
  ) #MT: noch zu entscheiden, wohin mit ICA, erstmal shortcut
  ProjectedPoints = as.matrix(res$points)
  Stress = res$stress
  
  if (PlotIt) {
    if (missing(Cls)) {
      AnzData = nrow(DataDists)
      Cls = rep(1, AnzData)
    }
    
    string = paste0('Sammons Mapping with Stress ', round(Stress, 2))
    #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string)
    PlotProjectedPoints(ProjectedPoints, Cls, main = string)
  }
  return(list(ProjectedPoints = ProjectedPoints, Stress = Stress))
  
}