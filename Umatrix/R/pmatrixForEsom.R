pmatrixForEsom=function(Data,Weights,Lines,Columns, Radius=NULL,PlotIt=FALSE, Toroid = T) {
#  PMatrix =  pmatrixForEsom(Data,EsomNeurons,Radius,Umatrix,PlotIt)
#INPUT
# Data[1:n,1:d]                      A \code{[n,k]} matrix containing the data
# EsomNeurons[1:Lines,1:Columns,1:weights]    Information stored as a List of weights in a 2D matrix,
#                                   oder das Rlistenobjekt aus ReadEsomNeurons
# Radius                            The radius for measuring the density within the hypersphere
# Umatrix[1:Lines,1:Columns]        Matrix with U-Matrix Height values
#Optional
# PlotHist
# Weights[1:Lines*1:Columns,1:weights]    If EsomNeurons is missing, Information stored as a List of weights
# Lines
# Columns

#OUTPUT
# UstarMatrix[1:Lines,1:Columns]
#
# MT 03/2015 aus matlab uebernommen
# MT01/17: Rollback von vorheriger Version, Neue Arguemente nach hinten verschoben
# FL01/17: Parameter EsomNeurons entfernt. Alte Aufrufe muessen gegebenenfalls zu Weights konvertiert werden
# FL07/17: Vorschlag Radius stammt jetzt aus EM und nicht aus ABC Analyse

#x=as.matrix(dist(Data))

if(missing(Lines)|missing(Columns))
  stop("Lines and Columns are necessary")
if(missing(Weights))
  stop("Weights are necessary")

EsomNeurons=ListAsEsomNeurons(Weights,Lines,Columns)

#########
# FL: Radius mittels EM schaetzen
##########
if(is.null(Radius)){
  ##############
  # Bestmatches aus Weights und Data generieren
  ##############
  #dataForRadius = 1:nrow(Data)
  #if(length(dataForRadius)>1000) dataForRadius = sample(dataForRadius, 1000)
  
  V = esomProjection(Weights, Data, Columns, Lines)
  #V = esomProjection(Weights, Data[dataForRadius,], Columns, Lines)
  BestMatches = cbind(1:nrow(V), V)

  ##########
  # Bestimme alle Distanzen und schlage radius durch EM vor
  ###########
  distancesAsMatrix = as.matrix(dist(Data))
  distances=distancesAsMatrix[lower.tri(distancesAsMatrix, diag = FALSE)]
  PradiusMinimum = signif(min(distances), digits=2)
  PradiusMaximum = signif(max(distances), digits=2)

  V = Delaunay4BestMatches(BestMatches, MatrixOrSize = c(Lines,Columns), IsToroid = Toroid)
  Delaunay = V$Delaunay
  ToroidDelaunay = V$ToroidDelaunay

  DelaunayDistances = Delaunay*distancesAsMatrix # nur distanzen zwischen nachbarn
  neighbourDistances = DelaunayDistances = DelaunayDistances[DelaunayDistances!=0]

  Radius = 1
  tryCatch({
    V = EMGauss(neighbourDistances, 2)
    V <- BayesDecisionBoundaries(V$Means, V$SDs, V$Weights)
    if(length(V)>1){
      print("Multiple Decision boundaries found. max(Distances)/2 chosen a Radius")
      Radius <- signif(max(distances)/2,2)
    }
    else Radius <- signif(V * 0.8, 2)
  }, error = function(e){
      print("No Decision boundaries found. max(Distances)/2 chosen a Radius")
      Radius <- signif(max(distances)/2,2)
  })

}

d=dim(EsomNeurons)[3]
UmatrixLines=dim(EsomNeurons)[1]
UmatrixCols=dim(EsomNeurons)[2]

if(is.null(d)){#EsomNeurons als liste
  stop('esom EsomNeurons has to be an array[1:Lines,1:Columns,1:Weights], use ListAsEsomNeurons')
}

AnzData=nrow(Data)
AnzVariablen=ncol(Data)

norm=2
distances=array(0,c(UmatrixLines,UmatrixCols,AnzData))

PMatrix <- matrix(0,UmatrixLines,UmatrixCols)
for ( i in 1:UmatrixLines ){
  for (j in 1:UmatrixCols){
    aux <- t(Data) - EsomNeurons[i,j,] # columns
    distances[i,j,]= sqrt(colSums(aux^norm)) #quadrierte Differenzen zu  EsomNeurons(i)
    PMatrix[i,j] <- sum(distances[i,j,] <= Radius)
  }
}

# plotten
if(PlotIt){
  p <- plotMatrix(PMatrix, ColorStyle="Pmatrix", Toroid=TRUE, Nrlevels = 10, TransparentContours = T) +
    ggtitle('P-Matrix')
  print(p)
#	optNrOfBins = OptimalNoBins(EsomNeuronsAnzInKugel)
#	minData = min(EsomNeuronsAnzInKugel,na.rm = TRUE)
#	maxData = max(EsomNeuronsAnzInKugel,na.rm = TRUE)
#	i = maxData-minData
#	optBreaks = seq(minData, maxData, i/optNrOfBins) # bins in fixed intervals
#	hist(EsomNeuronsAnzInKugel, breaks=optBreaks,xlab='# in spheres',main='distribution of # in spheres')
} # if
#plotPmxTopView(PMatrix,tiling=T)
  return(PMatrix)
}
