esomTrainOnline<-function(WeightVectors,Data, StartRadius=10, EndRadius=1, Columns=80, Lines=50, StartLearningRate=1, EndLearningRate=1,
                          NeighbourhoodFunction="cone", Toroid=TRUE, Epochs=10, NeighbourhoodCooling="linear",
                          LearningRateCooling="linear", ShinyProgress=NULL){
# esomTrainOnline(WeightVectors, Data)
# Trains the WeightVectors based on Data
#
# INPUT
# WeightVectors(1:m,1:n)	      WeightVectors that will be trained
#				      n weights with m components each
# Data(1:m,1:n)		      vectors to be projected with WeightVectors
#				      n datapoints with m components each
# OPTIONAL
# StartRadius			      Start Value for the Radius in which will be searched for neighbours
# EndRadius			      End Value for the Radius in which will be searched for neighbours
# Lines                          Height of the grid
# Columns                           Width of the grid
# StartLearningRate                   startvalue for LearningRate
# EndLearningRate		      endvalue for LearningRate
# NeighbourhoodFunction                         Method of training / kind of neighbourhood
# Toroid			      should the grid be considered as a Toroid
# Epochs			      number of Epochs in which every DataPoint will be used for training once
# NeighbourhoodCooling		      cooling method for radius. "linear" is the only available option at the moment.
# LearningRateCooling		      cooling method for learningRate. "linear" is the only available option at the moment.
# ShinyProgress          generate progress output for shiny if Progress Object is given
# OUTPUT
# result: WeightVectors(1:m,1:n)      the adjusted Weight Vectors
# author: Florian Lerch
# esomTrainOnline(WeightVectors, Data)


  distance = "euclidC"

  # load the Rcpp package
  #if(!require(Rcpp)){
  #  install.packages('Rcpp')
  #  library(Rcpp)
  #}
  #else{
  #  library(Rcpp)
  #}

  # load the cpp Functions that will be used
  #path = paste0(.pfad,'/Umatrix/')
  #source(paste0(path,'R/RcppExports.R'))
  #sourceCpp(paste0(path,'src/addRowWiseC.cpp'))

  #sourceCpp(paste0(path,'src/bestmatchC.cpp'))
  #sourceCpp(paste0(path, 'src/esomTrainedWeightVectorsConeC.cpp'))
  #sourceCpp(paste0(path, 'src/esomTrainedWeightVectorsGaussC.cpp'))
  #sourceCpp(paste0(path, 'src/esomTrainedWeightVectorsMexicanHatC.cpp'))


  # initialize Radius and LearningRate
  Radius = StartRadius
  LearningRate = StartLearningRate

  # initialize Progress Object for shiny
  if(!is.null(ShinyProgress))
    ShinyProgress$set(message = "Train Esom", value = 0)

  for(i in 1:Epochs){
		# permutation of Data
		ind <- sample(1:nrow(Data),nrow(Data))
		Data <- Data[ind,]

    # cool down Radius and LearningRate at the beginning of the epoch
    if(NeighbourhoodCooling=="linear") Radius <- coolDownLinear(StartRadius,EndRadius,Epochs,i)
		else if(NeighbourhoodCooling == "Lead In Lead Out") Radius <- coolDownLeadInLeadOut(StartRadius,EndRadius,Epochs,i)
    else stop("The neighbourhoodcooling is not recognized")

    # cool down learning rate
    if(LearningRateCooling=="linear") LearningRate <- coolDownLinear(StartLearningRate,EndLearningRate,Epochs,i)
		else if(LearningRateCooling=="Lead In Lead Out") LearningRate <- coolDownLeadInLeadOut(StartLearningRate,EndLearningRate,Epochs,i)
    else stop("The learningratecooling ist not recognized")

    # give feedback to user about the current status of execution
    print(paste0("Epoch: ",i," started"))

    # calculate neighbourhood pattern
    Pattern <- gridNeighbourhoodPattern(Radius = Radius)
    # find out if the same neighbour over both sides of the torus are possible
    RadiusBiggerThanTorus <- ((Radius*2) >= Columns) | ((Radius*2) >= Lines)

    for(DataPoint in 1:nrow(Data)){
      # execute the actual epoch
      WeightVectors <- esomTrainStep(WeightVectors, Data[DataPoint,], Radius, LearningRate, Pattern, NeighbourhoodFunction, Columns, Lines, Toroid, RadiusBiggerThanTorus)
    }

    # perecent of necessary epochs done
    if(!is.null(ShinyProgress))
      ShinyProgress$inc(1/Epochs, detail = paste("Epoch",i,"done"))
  }

  # give user feedback that training ist finished
  print("---- Esom Training Finished ----")

  WeightVectors
}
