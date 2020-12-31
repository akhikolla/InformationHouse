ICA = function(Data,OutputDimension=2,Contrastfunction="logcosh",Alpha=1,Iterations=200,PlotIt=FALSE,Cls){
  # Independent Component Analysis
  # projection=ICA(Data) 
  # INPUT 
  # Data[1:n,1:d]      array of data: n cases in rows, d variables in columns, matrix is not symmetric
  #                           or distance matrix, in this case matrix has to be symmetric
  # OPTIONAL
  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  #                           of the dimension chosen by cmdscale and OutputDimension
  # Contrastfunction          Maximierung der Negentropie ?ber geeignete geeignete Kontrastfunktion
  #                           Default: 'logcosh' G(u)=1/a*log cosh(a*u)
  #                           'exp': G(u)=-exp(u^2/2)
  #
  # Alpha                     constant with 1<=alpha<=2 used in approximation to neg-entropy when fun == "logcosh"
  # Iterations                maximum number of iterations to perform.
  #
  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  
# OUTPUT is a list with following elements:
# ProjectedPoints[1:n,OutputDimension]               n by OutputDimension matrix containing coordinates of the Projection:
#                                           with ICA transformed Data called Source, columns of Soruce contain the independent components
#
# Mixing[1:OutputDimension,1:d]             Mischungsmatrix s.d gilt Data=MixingMatrix*ProjectedPoints      
# Unmixing                                  Entmischungsmatrix mit Data*Unmixing=ProjectedPoints
# PCMatrix                                  pre-whitening matrix that projects data onto the first n.comp principal components.
#
# Note: Uses the R and C code implementation of the FastICA algorithm of Aapo Hyvarinen et al. 
# (http://www.cs.helsinki.fi/u/ahyvarin/) 
# Negentropie: Entropiedifferenz zu einer entsprechenden normalverteilten Zufallsvariable 
#               J(y)=|E(G(y)-E(G(v)))|^2
# author: MT 06/2015
  
  if (!requireNamespace('fastICA')) {
    message(
      'Subordinate projection package is missing. No computations are performed.
            Please install the package which is defined in "Suggests".'
    )
    return(
      list(
        Cls = rep(1, nrow(Data)),
        Object = "Subordinate projection package is missing.
                Please install the package which is defined in 'Suggests'."
      )
    )
  }

	if(missing(Data))
		stop('No Data given')
	Data;
	if(!is.matrix(Data))
		stop('Data has to be a matrix, maybe use as.matrix()')
		
    AnzVar=ncol(Data)
    AnzData=nrow(Data)
res=fastICA::fastICA(X=Data,n.comp=OutputDimension,fun = Contrastfunction,alg.typ = "parallel",alpha = Alpha,
            method = "C",row.norm = FALSE, maxit = Iterations,tol = 0.0001, verbose = TRUE) 

ProjectedPoints=res$S
if(PlotIt){
   if(missing(Cls)){
		AnzData=nrow(Data)
		Cls=rep(1,AnzData)
	}  
  
  string=paste0('ICA projection with ',Contrastfunction, ' approximation')
  #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string,Xlabel='independent component 1',Ylabel='independent component 2')
  PlotProjectedPoints(ProjectedPoints,Cls,main=string)
 } 

return(list(ProjectedPoints=ProjectedPoints,Mixing=res$A,Unmixing=res$W,PCMatrix=res$K))

}