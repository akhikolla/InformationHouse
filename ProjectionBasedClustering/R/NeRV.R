NeRV <- function(Data, lambda = 0.1, neighbors = 20, iterations = 10, cg_steps = 2, cg_steps_final = 40,randominit=T, OutputDimension=2,PlotIt=FALSE,Cls){
  # NeRV(Data,lambda)
  # Use the NeRV projection with matrix Data and lambda. Lambda controls the trustworthiness-continuity tradeoff
  #
  # INPUT
  # Data		Matrix of the Data to be projected
  #
  # OPTIONAL
  # lambda			    Controls the trustworthiness-continuity tradeoff. Default = 0.1
  # neighbors		    Set the number of nearest neighbours that each point should have. Should be positive. Default = 20
  # iterations 	    The number of iterations to perform. Default = 10
  # cg_steps		    The number of conjugate gradient steps to perform per iteration in NeRV's optimization scheme. Default = 2
  #	cg_steps_final	The number of conjugate gradient steps to perform on the final iteration in NeRV's optimization scheme. Default = 40
  # randominit =TRUE: Random init, ==FALSE PCA init
  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  
  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  # 
  # OUTPUT
  # ProjectedPoints[1:n,OutputDimension]                   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
  #
  #
  # AUTOR: MT 07/16,
  # Note: wrapper for c_nerv() of Felix Pape, which wrapes C++ SourceCode of Kaski 2010


  stepsPerRound=cg_steps
  stepsOnLastRound=cg_steps_final

  ProjectedPoints=c_NeRV(Data, lambda, neighbors, iterations, stepsPerRound, stepsOnLastRound, randominit, OutputDimension, PCA)
  if(PlotIt){
    AnzData=nrow(Data)
    if(missing(Cls))  Cls=rep(1,AnzData)
    string=paste0('NeRV with lambda ',round(lambda,4),' and neighbors ',neighbors)
    PlotProjectedPoints(ProjectedPoints,Cls,main=string)
  }   
  return(ProjectedPoints=ProjectedPoints)
  }
