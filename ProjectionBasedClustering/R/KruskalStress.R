KruskalStress <- function(InputDistances, OutputDistances){
  # KruskalStress(inDist, outDist)
  # Calculates the stress as defined by Kruskal for 2 distance matrices
  # INPUT
  # InputDistances    Distancematrix of the original Data
  # OutputDistances   Distancematrix of the projected Data
  #
  # OUTPUT
  # outStress         the Kruskal stress.
  stress <- sum((InputDistances-OutputDistances)^2)
  ssdOut <- sum(OutputDistances^2)
  outStress <- sqrt(stress/ssdOut)
  
  return(outStress)
  
}