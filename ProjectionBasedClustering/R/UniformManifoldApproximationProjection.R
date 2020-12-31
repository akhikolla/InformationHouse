UniformManifoldApproximationProjection=function(DataOrDistances,k,Epochs,OutputDimension=2,Algorithm='umap_pkg',PlotIt=FALSE,Cls,...){
  

  
  if(missing(DataOrDistances))
    stop('No DataOrDistances given')
  DataOrDistances;
  if(!is.matrix(DataOrDistances))
    stop('DataOrDistances has to be a matrix, maybe use as.matrix()')
  
  switch(Algorithm,
         umap_pkg={
           
           if (!requireNamespace('umap')) {
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
           
  custom.settings = umap::umap.defaults
  custom.settings$n_components =OutputDimension
  
  if(isSymmetric(unname(DataOrDistances))){
    custom.settings$input="dist"
  }else{ #!isSymmetric
    custom.settings$input="data"
  }# end if(isSymmetric(DataOrDistances))
  
  if(!missing(k))
    custom.settings$n_neighbors =k
  
  
  if(!missing(Epochs))
    custom.settings$n_epochs = Epochs
  
  dots=list(...)
  Names=names(dots)
  
  if(length(Names)>0){
    for(i in 1:length(Names)){
      tryCatch({
        custom.settings[[Names[i]]] <- dots[[i]]
      },error=function(e){
        warning(e)
      })
    }
    print(custom.settings)
  }
 
  proj=umap::umap(d=DataOrDistances,config = custom.settings,method = "naive")
  ProjectedPoints=proj$layout
         }, 
  uwot_pkg={
    if (!requireNamespace('uwot')) {
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
    
    if(isSymmetric(unname(DataOrDistances))){
      DataOrDistances=as.dist(DataOrDistances)
    }else{ #!isSymmetric
      #nothing
    }# end if(isSymmetric(DataOrDistances))
    custom.settings=list()
    #defaults of uwot
    if(missing(k))
      custom.settings$n_neighbor =15
    
    
    if(missing(Epochs))
      custom.settings$n_epochs = 500
    
    proj=uwot::umap(X=DataOrDistances,n_neighbors = custom.settings$n_neighbor,n_epochs = custom.settings$n_epochs,n_components = OutputDimension,...)
    
    ProjectedPoints=proj

  })
 
  
  
  if(PlotIt){
    if(missing(Cls)){
      AnzData=nrow(DataOrDistances)
      Cls=rep(1,AnzData)
    }  
    string=paste0('Uniform Manifold Approximation Projection with knn = ',custom.settings$n_neighbor,", and Epochs = ",custom.settings$n_epochs)

    PlotProjectedPoints(ProjectedPoints,Cls,main=string)
  } 
  
  return(list(ProjectedPoints=ProjectedPoints,ModelObject=proj,Setting=custom.settings))
}