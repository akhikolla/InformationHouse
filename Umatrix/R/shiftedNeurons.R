shiftedNeurons <- function(Weights, Lines, Columns, rows, cols){
#  requireRpackage("abind")
requireNamespace('abind')
  
  EsomNeurons = ListAsEsomNeurons(Weights, Lines, Columns)
  nrows = dim(EsomNeurons)[1]
  ncols = dim(EsomNeurons)[2]

  rows = -rows
  cols = -cols

  rows <- rows%%nrows
  cols <- cols%%ncols

  if(rows != 0){
    upperPart <- EsomNeurons[1:(nrows-rows),,]
    lowerPart <- EsomNeurons[(nrows-rows+1):nrows,,]
    EsomNeurons <- abind::abind(lowerPart, upperPart, along = 1)
  }

  if(cols != 0){
    leftPart <- EsomNeurons[,1:(ncols-cols),]
    rightPart <- EsomNeurons[,(ncols-cols+1):ncols,]
    EsomNeurons <- abind::abind(rightPart, leftPart, along = 2)
  }

  Weights = EsomNeuronsAsList(EsomNeurons)
  return(Weights)
}
