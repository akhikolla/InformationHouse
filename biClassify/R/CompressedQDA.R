formGroupCov <- function(TrainData, TrainCat, Group){
  if(Group == 1){
    Train <- TrainData[TrainCat == 1, ]
  }
  else if(Group == 2){
    Train <- TrainData[TrainCat == 2, ]
  }
  return(cov(Train))
}

QDAclassify <- function(TrainData, TrainCat, TestData, Cov1 = NULL, Cov2 = NULL, gamma = 1E-5){
  #----- Form Group Covariances --------
  if(is.null(Cov1)){
    Cov1 <- formGroupCov(TrainData = TrainData,
                        TrainCat = TrainCat,
                        Group = 1)
    Cov1 <- Cov1 + gamma * diag(ncol(Cov1))
  }

  if(is.null(Cov2)){
    Cov2 <- formGroupCov(TrainData = TrainData,
                        TrainCat = TrainCat,
                        Group = 2)
    Cov2 <- Cov2 + gamma * diag(ncol(Cov2))
  }

  #------ find Group sizes --------------
  n1 <- nrow(TrainData[TrainCat == 1, ])
  n2 <- nrow(TrainData[TrainCat == 2, ])
  n <- n1 + n2

  #------ Compute Group Means ------------
  Means1 <- colMeans(TrainData[TrainCat == 1, ])
  Means2 <- colMeans(TrainData[TrainCat == 2, ])

  #------ Compute Mahalnobis Distances -------
  logDet1 <- determinant(Cov1, logarithm = T)$modulus
  logDet2 <- determinant(Cov2, logarithm = T)$modulus

  Mahala1 <- t(solve(as.matrix(Cov1)) %*% (t(TestData) - as.numeric(Means1)))
  Mahala1 <- rowSums(t(t(TestData) - Means1) * Mahala1)
  Mahala1 <- Mahala1 + logDet1 - 2 * log(n1 / n)

  Mahala2 <- t(solve(as.matrix(Cov2)) %*% (t(TestData) - as.numeric(Means2)))
  Mahala2 <- rowSums(t(t(TestData) - Means2) * Mahala2)
  Mahala2 <- Mahala2 + logDet2 - 2*log(n2 / n)

  #------ Compare Mahalanobis Distances ----------
  Labels <- as.numeric(Mahala1 < Mahala2)
  Labels[Labels == 0] <- 2
  return(Labels)
}

compressQDApredict <- function(TrainData, TrainCat, TestData, m1, m2, s = 0.01, gamma = 1E-5){
  #------- Compute Group Means -------------
  Means1 <- colMeans(TrainData[TrainCat == 1, ])
  Means2 <- colMeans(TrainData[TrainCat == 2, ])

  #------- Compress Data --------
  compTrain <- compressData(X1 = TrainData[TrainCat == 1, ],
                            X2 = TrainData[TrainCat == 2, ],
                            m1 = m1,
                            m2 = m2,
                            s = s)
  compTrain <- as.matrix(rbind(t(t(as.matrix(compTrain$Cat1))+Means1),
                               t(t(as.matrix(compTrain$Cat2))+Means2)))
  compCat <- c(rep(1, m1), rep(2, m2))


  #------- Predict Labels -------------
  Labels <- QDAclassify(TrainData = compTrain,
                        TrainCat = compCat,
                        TestData = TestData,
                        gamma = gamma)
  return(Labels)
}

subsampleQDA <- function(TrainData, TrainCat, TestData, m1, m2, gamma = 1E-5){
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)

  Cat1Index <- sample(1:n1, m1)
  Cat2Index <- sample(1:n2, m2)

  #------ subset Data -------
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  SubX1 <- TrainX1[Cat1Index, ]
  SubX2 <- TrainX2[Cat2Index, ]

  subTrain <- as.matrix(rbind(SubX1, SubX2))
  subCat <- c(rep(1, m1), rep(2, m2))

  #------- Perform QDA on subset Data ------
  Labels <- QDAclassify(TrainData = subTrain,
                        TrainCat = subCat,
                        TestData = TestData,
                        gamma = gamma)
  return(Labels)
}





#' @title Quadratic Discriminant Analysis (QDA)
#' @description A wrapper function for the various QDA implementations available in this package.
#' @param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
#' @param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
#' @param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
#' @param Method A string of characters which determinds which version of QDA to use. Must be either "Full", "Compressed", or "Subsampled". 
#' @param Mode A string of characters which determines how the reduced sample paramters will be inputted for each method. Must be either "Research", "Interactive", or "Automatic". Default is "Automatic".
#' @param m1 The number of class 1 compressed samples to be generated. Must be a positive integer.
#' @param m2 The number of class 2 compressed samples to be generated. Must be a positive integer.
#' @param m The number of total compressed samples to be generated. Must be a positive integer.
#' @param s The sparsity level used in compression. Must satify 0 < s < 1.
#' @param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
#' @description Generates class predictions for \code{TestData}.
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Compressing large sample data for discriminant analysis'' arXiv preprint arXiv:2005.03858 (2020).
#' @details Function which handles all implementations of LDA. 
#' @examples 
#' TrainData <- QDA_Data$TrainData
#' TrainCat <- QDA_Data$TrainCat
#' TestData <- QDA_Data$TestData
#' plot(TrainData[,2]~TrainData[,1], col = c("blue","orange")[as.factor(TrainCat)])
#' 
#' #----- Full QDA -------
#' QDA(TrainData = TrainData,
#'     TrainCat = TrainCat,
#'     TestData = TestData,
#'     Method = "Full",
#'     gamma = 1E-5)
#'   
#' #----- Compressed QDA -------  
#'  m1 <- 700
#'  m2 <- 300
#'  s <- 0.01
#'  QDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Compressed",
#'      Mode = "Research",
#'      m1 = m1,
#'      m2 = m2,
#'      s = s,
#'      gamma = 1E-5)
#'      
#'  QDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Compressed",
#'      Mode = "Automatic",
#'      gamma = 1E-5)
#'  
#'  #----- Sub-sampled QDA ------
#'  m1 <- 700
#'  m2 <- 300
#'  QDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Subsampled",
#'      Mode = "Research",
#'      m1 = m1,
#'      m2 = m2,
#'      gamma = 1E-5)
#'      
#'  QDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Subsampled",
#'      Mode = "Automatic",
#'      gamma = 1E-5)
#'      
#' @return \item{Predictions}{(m x 1) Vector of predicted class labels for the data points in \code{TestData}.}  
#' @export
QDA <- function(TrainData, TrainCat, TestData, Method = "Full", Mode = "Automatic", m1 = NULL, m2 = NULL, m = NULL, s = NULL, gamma = 1E-5){
  if(Mode == "Automatic"){
    m1 <- sum(TrainCat == 1)/10
    m2 <- sum(TrainCat == 2)/10
    m <- length(TrainCat)/10
    s <- 1/sqrt(length(TrainCat))
  }
  
  if(Method == "Full"){
    return (QDAclassify(TrainData = TrainData,
                        TrainCat = TrainCat,
                        TestData = TestData,
                        gamma = gamma))
  }
  
  else if(Method == "Compressed"){
    if(Mode == "Interactive"){
      m1 <- readline(prompt="Please enter the number of compressed group 1 samples: ")
      m2 <- readline(prompt="Please enter the number of compressed group 2 samples: ")
      s <- readline(prompt="Please enter sparsity level s used in compression: ")
      
      m1 <- as.integer(m1)
      m2 <- as.integer(m2)
      s <- as.numeric(s)
      testParameters(m1 = m1, m2 = m2, s = s)
    }
    
    #----- Produce output ----------
    output <- compressQDApredict(TrainData = TrainData,
                                 TrainCat = TrainCat,
                                 TestData = TestData,
                                 m1 = m1, 
                                 m2 = m2,
                                 s = s,
                                 gamma = gamma)
    return(output)
  }
  
  else if(Method == "Subsampled"){
    if(Mode == "Interactive"){
      m1 <- readline(prompt="Please enter the number of group 1 sub-samples: ")
      m2 <- readline(prompt="Please enter the number of group 2 sub-samples: ")
     
      
      m1 <- as.integer(m1)
      m2 <- as.integer(m2)
    }
    
    output <- subsampleQDA(TrainData = TrainData,
                            TrainCat = TrainCat,
                            TestData = TestData,
                            m1 = m1,
                            m2 = m2,
                            gamma = gamma)
    return(output)
  }
  else{
    return ("Error: No version of QDA was run. Please check that Method is correct.")
  }
}

