neighbourhoodThroughPattern <- function(Index, Pattern, Columns, Lines, Toroid=TRUE, RadiusBiggerThanTorus=TRUE){
# neighbourhoodThroughPattern(Index, Pattern, Columns, Lines, Toroid)
# gets a Pattern for neighbourhood and returns all indices within that pattern starting from
# Index
#
# INPUT
# Index			position of a weightvector for which the neighbours should be returned. Index in form (lines,columns)
# Pattern			The neighbourhoodpattern. (see gridNeighbourhoodPattern)
# Columns			Width of the grid
# Lines			Height of the grid
# OPTIONAL
# Toroid			Is the Pattern build on a Toroid?
# RadiusBiggerThanTorus		If the currently used radius is bigger than min(Columns/2,Lines/2)
#				        there needs to be an expensive check for neighbours on two
#				        sides over the Toroid
# OUTPUT
# Neighbourhood(1:m,1:2)       indices of all weights in the Neighbourhood of Index together with their distance
# author: Florian Lerch
# neighbourhoodThroughPattern(c(3,5), Pattern, 80, 50)

  # move points according to difference
  Neighbourhood <- addRowWiseC(Pattern,c(Index,0))

  if(Toroid==FALSE){ # just keep points within the borders 
    Neighbourhood <- subset(Neighbourhood, Neighbourhood[,1] >= 1 & Neighbourhood[,2] >= 1 & 
                              Neighbourhood[,1] <= Lines & Neighbourhood[,2] <= Columns)
    rows = Neighbourhood[,1]
    cols = Neighbourhood[,2]
  }
  else if(Toroid==TRUE){ 
    # use modulo operator to get all entries into the borders
    rows = (Neighbourhood[,1]) %% Lines 
    cols = (Neighbourhood[,2]) %% Columns 
    
    rows[(rows == 0)] = Lines
    cols[(cols == 0)] = Columns
  }

  indices <- (rows-1)*Columns + cols

  Neighbourhood <- cbind(indices,Neighbourhood[,3])

  #RadiusBiggerThanTorus = F
  # on a Toroid its possible that a value exists two times in the Neighbourhood
  if(RadiusBiggerThanTorus){
    Neighbourhood = Neighbourhood[!duplicated(Neighbourhood[,1]),]
  }

  Neighbourhood
}


