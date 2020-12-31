
CrossValidation <- function(X,Y,F) {

  DrugNumber=nrow(Y)
  Index=1:DrugNumber

  FoldedIndex=NULL
  FF=F-1
  for (i in 1:FF) {
    FoldedIndex[[i]]=sample(Index,floor(DrugNumber/F))
    Index=setdiff(Index,FoldedIndex[[i]])
    if (i==FF) {
      FoldedIndex[[i+1]]=Index
    }
  }
  TestingIndex=NULL
  TrainingIndex=NULL
  TrainingData=NULL
  TestingData=NULL
  OutputTrain=NULL
  OutputTest=NULL
  for (i in 1:F){
    TestingIndex[[i]]=FoldedIndex[[i]]
    TrainingIndex[[i]]=setdiff(1:DrugNumber,TestingIndex[[i]])
    ## Input
    TrainingData[[i]]=X[TrainingIndex[[i]],,drop=FALSE]  #rows of training data of matrix X
    TestingData[[i]]=X[TestingIndex[[i]],,drop=FALSE]      #rows of testing data of matrix X
    ## Output
    OutputTrain[[i]]=Y[TrainingIndex[[i]],,drop=FALSE]        #rows of training data of matrix Y
    OutputTest[[i]]=Y[TestingIndex[[i]],,drop=FALSE]           #rows of testing data of matrix Y
  }
  result=list(TrainingData,TestingData,OutputTrain,OutputTest,FoldedIndex)
  return(result)
}
