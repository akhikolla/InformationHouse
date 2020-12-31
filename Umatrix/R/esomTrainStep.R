esomTrainStep <- function(WeightVectors, DataPoint, Radius, LearningRate, Pattern, NeighbourhoodFunction="cone", Columns=80, Lines=50, Toroid=TRUE, RadiusBiggerThanTorus=TRUE){
# esomTrainStep(WeightVectors, DataPoint, Radius, LearningRate, Pattern, NeighbourhoodFunction, Columns, Lines, Toroid)
# Does the training for one epoch of the esom. (Internal Function)
#
# INPUT
# WeightVectors(1:m,1:n)	      WeightVectors that will be trained
#				      n weights with m components each
# DataPoint(1:m)              Vector with which the WeightVectors will be trained
# NeighbourhoodFunction                         Method of training / kind of neighbourhood
#				      Options are: "cone", "mexicanhat" and "gauss"
# Radius			      the current Radius that should be used to find neighbours to the bm
# Columns                           Width of the grid
# Lines                          Height of the grid
# LearningRate			      A scaling factor for learning
# Pattern			      a "Pattern" of neighbouring indices relative to a position (0,0). It is used
#				      used for performance reasons and you get it by using gridNeighbourhoodPattern
# Toroid			      should the grid be considered as a Toroid
# RadiusBiggerThanTorus		      is the currently used Radius bigger than min(Columns/2,Lines/2)
#				      if this is the case, there needs to be an expensive check for neighbours on two
#				      sides over the Toroid
#
# OUTPUT
# WeightVectors			      the adjusted Weight Vectors
# author: Florian Lerch
# esomTrainStep(WeightVectors, DataPoint, Radius, LearningRate, Pattern)

  # euclidean distance for bestmatch search
  bm = bestmatchC(WeightVectors,DataPoint)

  # get position on grid out of index in WeightVectors array
  posByWeightInd <- function(ind,Columns){
    row = ((ind-1) %/% Columns) + 1
    col = ((ind-1) %% Columns) + 1
    c(row,col)
  }
  bmPos <- posByWeightInd(bm,Columns)

  # get neighbourhood of bm
  neighbourhood <- neighbourhoodThroughPattern(bmPos, Pattern, Columns, Lines, Toroid, RadiusBiggerThanTorus=RadiusBiggerThanTorus)
  
  # train vectors of neighbourhood
  if(nrow(neighbourhood) > 0){
    if(NeighbourhoodFunction == "cone")
      esomTrainedWeightVectorsConeC(WeightVectors, DataPoint, neighbourhood[,1], neighbourhood[,2], Radius, LearningRate)
    else if(NeighbourhoodFunction == "gauss")
      esomTrainedWeightVectorsGaussC(WeightVectors, DataPoint, neighbourhood[,1], neighbourhood[,2], Radius, LearningRate)
    else if(NeighbourhoodFunction == "mexicanhat")
      esomTrainedWeightVectorsMexicanHatC(WeightVectors, DataPoint, neighbourhood[,1], neighbourhood[,2], Radius, LearningRate)
    #ncols = ncol(WeightVectors)

    #oldWeightVectors = WeightVectors[neighbourhood[,1],]

    #newWeightVectors <- matrix(nrow=nrow(neighbourhood), ncol=ncols)

    # train the vectors
    #for(i in 1:nrow(neighbourhood)){
      # calculate the trained vector
    #  if(NeighbourhoodFunction == "cone"){
	#newWeightVectors[i,] <- esomTrainedWeightVectorConeC(oldWeightVectors[i,], DataPoint, neighbourhood[i,2], Radius, LearningRate)
   #   }
    #  else if(NeighbourhoodFunction == "gauss"){
	#newWeightVectors[i,] <- esomTrainedWeightVectorGauss(oldWeightVectors[i,], DataPoint, neighbourhood[i,2], Radius, LearningRate)
   #   }
  #    else if(NeighbourhoodFunction == "mexicanhat"){
	#newWeightVectors[i,] <- esomTrainedWeightVectorMexicanHat(oldWeightVectors[i,], DataPoint, neighbourhood[i,2], Radius, LearningRate)
   #   }
    #  else if(NeighbourhoodFunction == "nothing"){
	#newWeightVectors[i,] <- oldWeightVectors[i,]
  #    }

   # }

    # update the weightvector
    #WeightVectors[neighbourhood[,1],] <- newWeightVectors[1:nrow(neighbourhood),]

  }

  WeightVectors
}
