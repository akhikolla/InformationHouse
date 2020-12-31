ProjectionPursuit = function(Data,OutputDimension=2,Indexfunction="logcosh",Alpha=1,Iterations=200,PlotIt=FALSE,Cls){
  # Projection Pursuit
  # projection=ProjectionPursuit(Data)
  # In the absence of a generative model for the data the algorithm can be used to find the projection
  # pursuit directions. Projection pursuit is a technique for finding 'interesting' directions in multidimensional
  #  datasets
  # INPUT
  # Data[1:n,1:d]      array of data: n cases in rows, d variables in columns, matrix is not symmetric
  #                           or distance matrix, in this case matrix has to be symmetric
  # OPTIONAL
  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  #                           of the dimension chosen by cmdscale and OutputDimension
  # Indexfunction             Kriterium f?r Minimierung
  #                           Default: 'logcosh' G(u)=1/a*log cosh(a*u) (ICA)
  #                           'exp': G(u)=-exp(u^2/2)
  #                           'kernel'  1/(1* pi )*exp(r/2)
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
  requireNamespace('fastICA')
  if (missing(Data))
    stop('No Data given')
  Data
  
  if (!is.matrix(Data))
    stop('Data has to be a matrix, maybe use as.matrix()')
  
  AnzVar = ncol(Data)
  AnzData = nrow(Data)
  warning('FactICA dokumentation on indexfunction unclear')
  res = fastICA::fastICA(
    X = Data,
    n.comp = OutputDimension,
    fun = Indexfunction,
    alg.typ = "deflation",
    alpha = Alpha,
    method = "R",
    row.norm = FALSE,
    maxit = Iterations,
    tol = 0.0001,
    verbose = TRUE
  )
  
  ProjectedPoints = res$S
  if (PlotIt) {
    if (missing(Cls))
      Cls = rep(1, AnzData)
    
    string = paste0('Projection pursuit with Indexfunction ', Indexfunction)
    
    PlotProjectedPoints(
      ProjectedPoints,
      Cls = Cls,
      xlab = 'independent component 1',
      ylab = 'independent component 2',
      main = string
    )
    # requireNamespace("dbt.ClassAnalysis")
    # dbt.ClassAnalysis::ClassPlot(
    #   ProjectedPoints[, 1],
    #   ProjectedPoints[, 2],
    #   Cls = Cls,
    #   Title = string,
    #   Xlabel = 'independent component 1',
    #   Ylabel = 'independent component 2'
    # )
    
  }
  
  return(
    list(
      ProjectedPoints = ProjectedPoints,
      Mixing = res$A,
      Unmixing = res$W,
      PCMatrix = res$K
    )
  )
  
}