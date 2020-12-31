esomTrain <- function(Data, Lines=50, Columns=82, Epochs=24, Toroid=TRUE, NeighbourhoodFunction="gauss", StartLearningRate=0.5, EndLearningRate=0.1,
		 StartRadius=24, EndRadius=1, NeighbourhoodCooling="linear", LearningRateCooling="linear", ShinyProgress = NULL,
		 ShiftToHighestDensity = TRUE, InitMethod="uni_min_max", Key = NULL, UmatrixForEsom=T){
# V <- esom(Data, Lines=50, Columns=82, Epochs=24, Toroid=TRUE, NeighbourhoodFunction="gauss", StartLearningRate=0.5, EndLearningRate=0.1, StartRadius=24, EndRadius=1, NeighbourhoodCooling="linear", LearningRateCooling="linear")
# BestMatches <- V$BestMatches
# Weights <- V$Weights
# initialize a grid, train it and project the data with the esom method
#
# INPUT
# Data(1:n,1:m)		Data that will be used for training, initialization and projection
#				  n datapoints with m attributes each
# OPTIONAL
# Lines                          Height of the grid
# Columns                           Width of the grid
# Epochs			      number of Epochs in which every trainingDataPoint will be used for training once
# Toroid			 should the grid be considered as a Toroid
# NeighbourhoodFunction               Method of training / kind of Neighbourhood.
#                                    Possible values are: "cone", "mexicanhat" and "gauss"
# StartLearningRate                   startvalue for LearningRate
# EndLearningRate		      endvalue for LearningRate
# StartRadius			      Start Value for the Radius in which will be searched for Neighbours
# EndRadius			      End Value for the Radius in which will be searched for Neighbours
# NeighbourhoodCooling		      cooling method for radius. "linear" is the only available option at the moment.
# LearningRateCooling		      cooling method for learningRate. "linear" is the only available option at the moment.
# ShinyProgress          generate progress output for shiny if Progress Object is given
# ShiftToHighestDensity   shifts the umatrix so that the point of highest density will be at the center
# InitMethod		name of the method that will be used to choose initializations
#				      Valid Inputs: uni_min_max  <- uniform distribution with minimum and maximum
#								    from sampleData
#						    uni_mean_2std <- uniform distribution based on mean and standard
#								     deviation of sampleData
#						    norm_mean_2std <- normal distribuation based on mean and standard
#								     deviation of sampleData
# Key(1:n)    Key Matching the data
# OUTPUT
# list with
# BestMatches(1:m,2)		      Projected BestMatches
# NeuronWeights(1:m,1:n)              List of trained weights
# EsomNeurons(1:Lines, 1:Columns, 1:m)   Trained Weights in Wts format (DEPRECATED!)
# author: Florian Lerch
# details: A Toroid can be considered an area, where opposing "borders" are connected

  distance = "euclidC"

  grid <- esomInit(Data,Lines,Columns, InitMethod)

    x <- esomTrainOnline(grid, Data, StartRadius, EndRadius,
		       Columns,Lines,
		       StartLearningRate, EndLearningRate,
                       NeighbourhoodFunction, Toroid,
		       Epochs, NeighbourhoodCooling,
		       LearningRateCooling, ShinyProgress = ShinyProgress)
  if(!is.null(ShinyProgress))
    ShinyProgress$set(message = "Generating Umatrix", value = 0)


  #EsomNeurons = ListAsEsomNeurons(x, Lines,Columns)
  Weights = x
  if(ShiftToHighestDensity & Toroid){
    # find point of highest density
    print("shift to point of highest density")

    Radius = mean(dist(Data))
    pmatrix = pmatrixForEsom(Data, Weights=Weights, Lines=Lines, Columns=Columns, Radius = Radius)

    pos <- which(pmatrix == max(pmatrix), arr.ind=T)[1,]

    Weights = shiftedNeurons(Weights, Lines, Columns, -pos[1], -pos[2])
    x <- Weights
  }

  # load cpp function for projection
  #path = paste0(SubversionDirectory(),'PUB/dbt/Umatrix/')
  #source(paste0(path,'R/RcppExports.R'))
  #sourceCpp(paste0(path,'src/bestmatchesC.cpp'))

  #projection <- esomProjection(x, Data, Columns, Lines)
  projection <- bestmatchesC(x, Data, Columns)
  if(!is.null(Key)) projection = cbind(Key, projection)
  else{
    warning("No Keys given. 1:nrow(Data) will be used as keys for Bestmatches.")
    projection = cbind(1:nrow(projection), projection)
  }

  Umatrix = NULL
  if(UmatrixForEsom)
    Umatrix = umatrixForEsom(x, Lines, Columns, Toroid)

  list(BestMatches = projection, Weights = x, Lines = Lines, Columns = Columns,
       Toroid = Toroid, Umatrix = Umatrix)
}
