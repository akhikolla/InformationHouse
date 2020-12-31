GeneratePmatrix=function(Data,EsomNeurons,...) {
#  PMatrix =  GeneratePmatrix(Data,EsomNeurons,Radius,Umatrix,PlotIt)
#INPUT
# Data[1:n,1:d]                      A \code{[n,k]} matrix containing the data
# EsomNeurons[1:Lines,1:Columns,1:weights]    Information stored as a List of weights in a 2D matrix,
#                                   oder das Rlistenobjekt aus ReadEsomNeurons
# Radius                            The radius for measuring the density within the hypersphere
# Umatrix[1:Lines,1:Columns]        Matrix with U-Matrix Height values


#OUTPUT
# UstarMatrix[1:Lines,1:Columns]
#
# MT 03/2015 aus matlab uebernommen
# MT01/17: Rollback von vorheriger Version, Neue Arguemente nach hinten verschoben
  if(missing(Data))
    stop("Data is necessary.")
  
x=as.matrix(dist(Data))

# if(missing(LC))
#   stop("Lines and Columns are necessary")
if(missing(EsomNeurons))
  stop("EsomNeurons are necessary.")

Lines=dim(EsomNeurons)[1]
Columns=dim(EsomNeurons)[2]


dots=list(...)

if(is.null(dots[["PlotIt"]]))
  PlotIt=FALSE
else
  PlotIt=TRUE


if(is.null(dots[["Radius"]])){

x=x[lower.tri(x, diag = FALSE)]
para=quantile(x,c(0.2)) #geschaetzter paretorRadius
	if(requireNamespace('ABCanalysis')){
		xx=ABCanalysis::ABCRemoveSmallYields(x,0.5)
		x=xx$SubstantialData
		res=suppressWarnings(ABCanalysis::ABCanalysis(x))
		Radius=1/(min(x[res$Aind])/max(x[res$Cind]))*para  #Verhaeltnis vermutliche inner/Inter Clusterdistanz
		#print(min(x[res$Aind])/max(x[res$Cind]))
		#print(par)
		#print(Radius)
	}else{
		warning('Package ABCanalysis is not available.Please install it.')
		Radius=quantile(x,c(0.4))
	}
}else{
  Radius=dots$Radius
}
#dimensions=dim(EsomNeurons)
#AnzEsomNeurons=dimensions[1]
#AnzVariablen=dimensions[2]

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
#EsomNeuronsAnzInKugel=vector(mode='numeric',AnzEsomNeurons)
PMatrix <- matrix(0,UmatrixLines,UmatrixCols)
for ( i in 1:UmatrixLines ){
  for (j in 1:UmatrixCols){
    aux <- t(Data) - EsomNeurons[i,j,] # columns
    distances[i,j,]= sqrt(colSums(aux^norm)) #quadrierte Differenzen zu  EsomNeurons(i)
		#EsomNeuronsAnzInKugel[i]=sum(distances[i,j,]<=Radius)
    PMatrix[i,j] <- sum(distances[i,j,] <= Radius)
  }
}

if(PlotIt){
  if(is.null(dots[["Bestmatches"]]))
    Bestmatches=matrix(1,2,2)#dummy
  else
    Bestmatches=dots$Bestmatches
  
  plotTopographicMap(PMatrix,Bestmatches, Colormap=DataVisualizations::PmatrixColormap, Tiled=TRUE,...)
  # p <- plotTopographicMap(PMatrix,Bestmatches, Colormap=DataVisualizations::PmatrixColormap, Tiled=TRUE) +
  # #   ggtitle('P-Matrix')
  # print(p)
}
  return(PMatrix)
}
