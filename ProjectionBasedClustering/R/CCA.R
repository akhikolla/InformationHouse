CCA = function(DataOrDistances,Epochs,OutputDimension=2,method='euclidean', alpha0 = 0.5, lambda0,PlotIt=FALSE,Cls){
  #  CCA Projects data vectors using Curvilinear Component Analysis.
  #  
  #  res = CCA(V$Data, Epochs=50,OutputDimension=2,PlotIt = T)
  #   
  # INPUT
  # DataOrDistances[1:n,1:d]      array of data: n cases in rows, d variables in columns, matrix is not symmetric
  #                           or distance matrix, in this case matrix has to be symmetric
  #  epochs                  (scalar) training length
  #
  # OPTIONAL
  # OutputDimension           data is projected onto a R^p where P is the maximum ( default ==2)
  # method                    method specified by distance string: 
  #                          'euclidean','cityblock=manhatten','cosine','chebychev','jaccard','minkowski','manhattan','binary' 
  #
  #    alpha0                (scalar) initial step size, 0.5 by default
  #    lambda0              (scalar) initial radius of influence, 3*max(std(D)) by default
  # PlotIt                    bool, defaut=FALSE, if =TRUE: ClassPlot of every current Position of Databots will be made.
  #                           OutputDimension>2 only the first two dimensions will be shown
  # cls                       vector, Classifikation of Data if available, ClassPlots will be colorized
  # 
  # OUTPUT is a list with following elements:
  # ProjectedPoints[1:n,OutputDimension]                   n by OutputDimension matrix containing coordinates of the Projection: A matrix of the fitted configuration.
  #
  # Error                                         CCA error
  
  # Note: 
  #  Unknown values (NaN's) in the data: projections of vectors with
  #  unknown components tend to drift towards the center of the
  #  projection distribution. Projections of totally unknown vectors are
  #  set to unknown (NaN).
  # 
  #  Implementation in Matlab:
  #   Reference: Demartines, P., Herault, J., "Curvilinear Component
  #     Analysis: a Self-Organizing Neural Network for Nonlinear
  #     Mapping of Data Sets", IEEE Transactions on Neural Networks, 
  #     vol 8, no 1, 1997, pp. 148-154.
  
  #   Contributed to SOM Toolbox 2.0, February 2nd, 2000 by Juha Vesanto
  #   Copyright (c) by Juha Vesanto
  #   http://www.cis.hut.fi/projects/somtoolbox/
  #
  #   juuso 171297 040100
  #
  #	Ue?bersetzung in R: Florian Lerch
  # 1.Editor: 06/2015 MT
  #
  #  %%%%%%%%%%%%%%%%%%%%%%%%%
  #  %% subfunctions
  #
  if (missing(DataOrDistances))
    stop('No DataOrDistances given')
  DataOrDistances
  
  if (!is.matrix(DataOrDistances))
    stop('DataOrDistances has to be a matrix, maybe use as.matrix()')
  
  
  if (missing(Epochs)){
    warning('scalar value for number of eppochs is missing. Setting epochs=20 which may not be prerable in order to continue the algorithm. There is no default setting for this parameter!')
    Epochs=20
  }else{
    epochs = Epochs
  }
  if (missing(lambda0))
    lambda0 = NULL
  
  if (isSymmetric(unname(DataOrDistances))) {
    Mdist = DataOrDistances
    AnzVar = ncol(DataOrDistances)
    AnzData = nrow(DataOrDistances)
  } else{
    #!isSymmetric
    AnzVar = ncol(DataOrDistances)
    AnzData = nrow(DataOrDistances)
    #DataDists=as.matrix(DistanceMatrix(X=DataOrDistances,method = method))
    Mdist = NULL
  }# end if(isSymmetric(DataOrDistances))
  
  squareform=function(X)
# y=squareform(x)
# analog zu matlab übernommen aus Rpaket
# Format or generate a distance matrix.
#
# INPUT
# X         numeric vector or matrix
#
# OUTPUT
# y      Returns a matrix if x is a vector, and a vextor if x is a matrix.

# Autor: [Package pracma version 1.6.4 Index], kopiert von MT
#
# Example:
# x <- 1:6
# y <- squareform(x)
# #  0  1  2  3
# #  1  0  4  5
# #  2  4  0  6
# #  3  5  6  0
# all(squareform(y) == x)
# # TRUE
{
requireNamespace('pracma')
  if(is.matrix(X)){
    if(isSymmetric(X)){
      return(X[upper.tri(X)])
    }else{
      return(pracma::squareform(X))
    }
  }else{
    return(pracma::squareform(X))
  }
}
  
  D = DataOrDistances
  P = OutputDimension
  
  
  potency_curve = function(v0, vn, l)
    return(v0 * (vn / v0) ^ ((0:(l - 1)) / (l - 1)))
  
  cca_error = function(P, Mdist, lambda) {
    noc = nrow(P)
    odim = ncol(P)
    noc_x_1 = rep(1, noc)
    odim_x_1 = rep(1, odim)
    
    error = 0
    for (i in 1:noc) {
      known = which(!is.nan(Mdist[, i]))
      if (length(known) > 0) {
        y = t(as.matrix(P[i, ]))
        Dy = P[known, ] - y[noc_x_1[known], ]
        dy = sqrt(rowSums(Dy ^ 2))
        fy = exp(-dy / lambda)
        error = error + sum(((Mdist[known, i] - dy) ^ 2) * fy)
      }
    }
    error = error / 2
    return(error)
  }
  
  # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  # Hier beginnt die eigentliche Funktion
  
  noc = nrow(D)   # Anzahl der Zeilen also Datensätze
  dim = ncol(D)   # Anzahl der Spalten also Attribute
  
  noc_x_1 = rep(1, noc)
  
  me = matrix(rep(0, dim), nrow = 1)
  st = matrix(rep(0, dim), nrow = 1)
  
  # durchschnitt und standardabweichung
  for (i in 1:dim) {
    me[i] = mean(D[which(is.finite(D[, i]) == TRUE), i])
    st[i] = sd(D[which(is.finite(D[, i]) == TRUE), i])
  }
  
  # mit zufälligen Werten initialisieren
  P = ((2 * matrix(runif(noc * P), noc) - 1) * st[noc_x_1, 1:P]) + me[noc_x_1, 1:P]
  
  dummy = nrow(P)
  odim = ncol(P)
  odim_x_1 = matrix(1, odim, 1)
  
  train_len = epochs * noc
  
  sample_inds = ceiling(runif(train_len, 0, noc))
  
  # falls keine Distanzmatrix angegeben: selber bilden
  if (is.null(Mdist))
    Mdist = as.matrix(dist(D, diag = TRUE, upper = TRUE))
  else {
    # falls Mdist in Squareform, bringe sie wieder in quadratische Form
    if (nrow(Mdist) == 1)
      Mdist = squareform(Mdist)
    if (nrow(Mdist) != noc)
      stop('Mutual distance matrix size and data set size do not match')
  }
  
  alpha = potency_curve(alpha0, alpha0 / 100, train_len)
  if (is.null(lambda0))
    lambda0 = max(st) * 3
  lambda = potency_curve(lambda0, 0.01, train_len)
  #
  #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  Action
  #
  k = 0
  print(sprintf('iterating: %d / %d epochs', k, epochs))
  
  for (i in 1:train_len) {
    ind = sample_inds[i]  # sample index
    dx = Mdist[, ind] # mutual distances in input space
    known = which(!is.nan(dx)) # known distances
    if (length(known) > 0) {
      # sample vector's projection
      y = t(as.matrix(P[ind, ]))
      # distances in output space
      Dy = P[known, ] - y[noc_x_1[known], ]
      dy = sqrt(rowSums(Dy ^ 2))
      # relative effect
      dy[which(dy == 0)] = 1 # to get rid of div-by-zero's
      
      fy = as.matrix(exp(-dy / lambda[i]) * (dx[known] / dy - 1))
      
      #      % Note that the function F here is e^(-dy/lambda))
      #      % instead of the bubble function 1(lambda-dy) used in the
      #      % paper.
      #
      #      % Note that here a simplification has been made: the derivatives of the
      #      % F function have been ignored in calculating the gradient of error
      #      % function w.r.t. to changes in dy.
      #
      #      % update
      P[known, ] = P[known, ] + alpha[i] * fy[, odim_x_1] * Dy
    }
    
    #     track
    if (i %% noc == 0) {
      k = k + 1
      print(sprintf('iterating: %d / %d epochs', k, epochs))
    }
    
  }
  #
  #  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  #  clear up
  #
  #  calculate error
  error = cca_error(P, Mdist, lambda[train_len])
  print(sprintf('%d iterations, error %f', epochs, error))
  
  #  set projections of totally unknown vectors as unknown
  unknown = which(sum(t(is.nan(D))) == dim)
  
  P[unknown, ] = NaN
  ProjectedPoints = P
  if (PlotIt) {
    if (missing(Cls)) {
      AnzData = nrow(ProjectedPoints)
      Cls = rep(1, AnzData)
    }
    
    string = paste0('CCA with error ', round(error, 4), ' and epochs ', Epochs)
    #ClassPlot(ProjectedPoints[,1],ProjectedPoints[,2],Cls=Cls,Title=string)
    PlotProjectedPoints(ProjectedPoints, Cls, main = string)
    
  }
  
  return(list(ProjectedPoints = ProjectedPoints, Error = error))
}
