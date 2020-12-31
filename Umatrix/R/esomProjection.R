esomProjection<-function(WeightVectors,Data,Columns,Lines){
# esomProjection(WeightVectors, Data, Columns, Lines)
# calculates the actual projection for Data with the WeightVectors
#
# INPUT
# WeightVectors(1:m,1:n)	      WeightVectors that will be trained
#				      n weights with m components each
# Data(1:m,1:n)		              vectors to be projected with WeightVectors
#				      n datapoints with m components each
# Lines                               Height of the grid
# Columns                             Width of the grid
#
# OUTPUT
# BestMatches(1:n,1:2)            matrix that contains n points with their positions on the grid
# author: Florian Lerch
# details: reimplemented from Databionic ESOM Tools (http://databionic-esom.sourceforge.net)
# esomProjection(WeightVectors, Data, 80, 50)
  # get the position on a grid of the closest vector out of WeightVectors to x (based on euclidean distance)

  closestVector <- function(x){
    distances <- apply(WeightVectors,1,function(y)sum((y-x)^2))

    # index of the closest Vector in WeightVectors
    bm = which.min(distances)

    # get the gridPositions out of the index in WeightVectors
    row = ((bm-1) %/% Columns) + 1
    col = ((bm-1) %% Columns) + 1

    c(row,col)
  }

  # get the closestVector for each vector out of Data
  projection <- t(apply(Data,1,closestVector))
  projection
}

