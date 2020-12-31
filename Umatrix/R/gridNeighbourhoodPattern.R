gridNeighbourhoodPattern <- function(Radius){
# gridNeighbourhoodPattern(Radius)
# returns index for neighbours relative to (0,0) and their distances
#
# INPUT
# Radius		      Radius in which neighbours should be searched
# 
# OUTPUT
# Neighbourhood(1:m,1:3)      positions of all weights in the Neighbourhood of c(0,0) together with their distance
# author: Florian Lerch
# gridNeighbourhoodPattern(3)

  # the pattern is starting from point 0,0
  startPoint = c(0,0) 
  
  # get all neighbours
  Neighbourhood <- c()

  # make sure that Radius is an integer
  Radius <- floor(Radius)
  
  # calculate only a quarter of the sphere and get the rest through symmetry
  for (y in (startPoint[1] - Radius):(startPoint[1])){
    for (x in (startPoint[2] - Radius):(startPoint[2])){
      if ((y - startPoint[1])^2 + (x - startPoint[2])^2 <= Radius^2) {
        ySym <- startPoint[1] - (y - startPoint[1])
        xSym <- startPoint[2] - (x - startPoint[2])
        
        # add the new found neighbours and the points symmetric to them into the Neighbourhood
	nn1 <- c(y, x, sqrt(y^2+x^2))
	nn2 <- c(y, xSym, sqrt(y^2+xSym^2))
	nn3 <- c(ySym , x, sqrt(ySym^2+x^2))
	nn4 <- c(ySym, xSym, sqrt(ySym^2+xSym^2))

        newNeighbours <-  matrix(c(nn1, nn2, nn3, nn4),ncol=3,byrow=TRUE)
        Neighbourhood <- rbind(Neighbourhood,newNeighbours)
      }
    }
  }
  
  # remove duplicated entries
  Neighbourhood <- unique(Neighbourhood)
  
  # sort by distance, this makes it easier to remove the farther neighbour on a Toroid
  Neighbourhood <- Neighbourhood[order(Neighbourhood[,3]),] 

  Neighbourhood
}

