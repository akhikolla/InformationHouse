Projection2Bestmatches=function(ProjectedPoints){
  n=nrow(ProjectedPoints)
  if(n>4096/8)
    minNeurons=n*8
  else
    minNeurons=4096
  
  if (requireNamespace('GeneralizedUmatrix',quietly = TRUE)) {
  coordsres=GeneralizedUmatrix::XYcoords2LinesColumns(X=ProjectedPoints[,1],Y=ProjectedPoints[,2],PlotIt=F,minNeurons = minNeurons)
  }else{
    message("Please Install GeneralizedUmatrix package.")
    return(list(Bestmatches="GeneralizedUmatrix has to be installed for this output",LC=NULL))
  }
  return(list(Bestmatches=coordsres$GridConvertedPoints,LC=coordsres$LC))
}