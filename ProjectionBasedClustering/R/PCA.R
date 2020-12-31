PCA  = function(Data,OutputDimension=2,Scale=FALSE,Center=FALSE,PlotIt=FALSE,Cls){
  # Performs a principal components analysis on the given data matrix
  # projection=PCA(Data) 
  # INPUT 
  # Data[1:n,1:d]      array of data: n cases in rows, d variables in columns
  # OPTIONAL

  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  #                           of the dimension chosen by cmdscale and OutputDimension
  # Scale                     a logical value indicating whether the variables should be scaled to
  #                           have unit variance before the analysis takes place. 
  #                           The default is FALSE for consistency with S, but in general scaling
  #                           is advisable. Alternatively, a vector of length equal the 
  #                           number of columns of x can be supplied. The value is passed to scale.#
  # Center                    a logical value indicating whether the variables should be shifted to
  #                           be zero centered. Alternately, a vector of length equal the number 
  #                           of columns of x can be supplied. The value is passed to scale
  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  
  # OUTPUT is a list with following elements:
  # ProjectedPoints[1:n,OutputDimension]   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
  # Rotation                      the matrix of variable loadings (i.e., a matrix whose columns contain the eigenvectors)
  # sDev                          the standard deviations of the principal components (i.e., the square roots of 
  #                               the eigenvalues of the covariance/correlation matrix, 
  #                               though the calculation is actually done with the singular values of the data matrix)
  # TransformedData[1:n,1:d]      matrix  with PCA transformed Data
  # Center                        the centering used, or FALSE
  # Scale                         the scaling used, or FALSE
  
  # author: MT 06/2015
  
	if(missing(Data))
		stop('No Data given')
	Data;
	if(!is.matrix(Data))
		stop('Data has to be a matrix, maybe use as.matrix()')
		
    AnzVar=ncol(Data)
    AnzData=nrow(Data)
  
res <- prcomp(x=Data,retx=T,scale. =Scale,tol = 0,center=Center)
TransData=as.matrix(res$x)
ProjectedPoints=TransData[,1:OutputDimension]
Rotation <- res$rotation
sDev <- res$sdev
CenterOut <- res$center
ScaleOut <- res$scale

if(PlotIt){
    if(missing(Cls)){
		AnzData=nrow(Data)
		Cls=rep(1,AnzData)
	}  
  
  string=paste0('PCA projection with sDev-X= ',round(sDev[1],2),' ;-Y ',round(sDev[2],2))
  #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string)
    
 PlotProjectedPoints(ProjectedPoints,Cls,main=string)
}                    
return(list(ProjectedPoints=ProjectedPoints,Rotation = Rotation, sDev = sDev, TransformedData=TransData,Center = CenterOut, Scale = ScaleOut))

}