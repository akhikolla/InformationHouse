# ------------ Compressed and Subset LDA Code -----------------

# m is the number of rows in the sketching matrix
# n is the number of rows in the data matrix
# p is the probability parameter in the sparse bernoulli
#' @importFrom stats rbinom
rademacherSketchMatrix <- function(n, m, s){
  if(n == 0) stop("n is set to 0")
  if(s == 0) stop("sparsity level s is set to 0")
  
  len = rbinom(n = 1, size = n*m, prob = s) #initialize number of non-zero entries
  
  # --- Give non-zero entries depending on type
  x <- sample(c(-1, 1), size = len, replace = T) #Sample len number of +/- 1s
  Indices <- sample(x = 1:(n*m) , size = len) #Select where non-zero entries go
  
  Q <- Matrix::sparseVector(x = x, i = Indices, length = n*m)
  Q <- Matrix::Matrix(data = Q,
                      nrow = m,
                      ncol = n,
                      sparse = T )
  
  return(Q / sqrt(n * s))
}

#' @importFrom stats rbinom
#' @importFrom stats rnorm
gaussianSketchMatrix <- function(n, m, s){
  if(n == 0) stop("n is set to 0")
  if(s == 0) stop("sparsity level s is set to 0")
  
  len = rbinom(n = 1, size = n*m, prob = s) #initialize number of non-zero entries
  
  # --- Give non-zero entries depending on type
  x <- rnorm(len, mean = 0, sd = 1) #Sample len number of +/- 1s
  Indices <- sample(x = 1:(n*m) , size = len) #Select where non-zero entries go
  
  Q <- Matrix::sparseVector(x = x, i = Indices, length = n*m)
  Q <- Matrix::Matrix(data = Q,
                      nrow = m,
                      ncol = n,
                      sparse = T )
  
  return(Q / sqrt(n * s))
}

countSketchMatrix <- function(n, m){
  # --- Random Sign Flips ---
  randSigns <- sample(c(1, -1), size = n, replace = TRUE) # a n-vector of +/- 1
  
  # --- row allocations ---
  RowAllocation <- sample(1:m, size = n, replace = TRUE)
  
  # ---
  Q <- Matrix::Matrix(data = 0,
                      nrow = m,
                      ncol = n,
                      sparse = T )
  for(i in 1:m){
    Q[i, RowAllocation == i] <- randSigns[RowAllocation == i]
  }
  return(sqrt(m/n)*Q)
}


createSketchMatrix <- function(n, m, s, type = "Rademacher"){
  # --- Create compresison matrix depending on type ---
  if(type == "Rademacher"){ # Rademacher sketch
    Q <- rademacherSketchMatrix(n = n, m = m, s = s)
  }
  
  else if(type == "Gaussian"){ # Gaussian sketch
    Q <- gaussianSketchMatrix(n = n, m = m, s = s)
  }
  
  else if(type == "Count"){ # Count Sketch
    Q <- countSketchMatrix(n = n, m = m)
  }
  else{
    stop("Please input correct type of compression matrix.")
  }
  return(Q)
}


#Compress the seperate groups of data X1 and X2 to dimension m1 and m2, respectively
#X1 is a n1 x p matrix
#X2 is a n2 x p matrix
# m1 is a positive integer. Dimension of projection for X1
# m2 is a positive integer. Dimension of projection for X2
compressData <- function(X1, X2, m1, m2, s, type = "Rademacher"){
  
  # --- Run preliminary dimension tests ---
  if(ncol(X1) != ncol(X2)) stop("Number of features of X1 and X2 different")
  if(m1 == 0) stop("m1 is set to 0.")
  if(m2 == 0) stop("m2 is set to 0.")
  
  n1 <- nrow(X1)
  n2 <- nrow(X2)
  ncols <- ncol(X1)
  
  # --- Compute group means ----
  Mean1 <- colMeans(X1)
  Mean2 <- colMeans(X2)
  
  # --- Create compresison matricies for both groups depending on type ---
  Q1 <- createSketchMatrix(n = n1, m = m1, s = s, type = type)
  Q2 <- createSketchMatrix(n = n2, m = m2, s = s, type = type)
  
  # --- Center Data ---
  X1 <- t(t(X1)-Mean1)
  X2 <- t(t(X2) - Mean2)
  
  QX1 <- Q1 %*% as.matrix(X1)
  QX2 <- Q2 %*% as.matrix(X2)
  return(list(Cat1 = QX1, Cat2 = QX2))
}



#' @importFrom stats cov
formWithinGroupCov <- function(TrainData, TrainCat, Means1 = NULL, Means2 = NULL){
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  X1 <- TrainData[TrainCat == 1, ]
  X2 <- TrainData[TrainCat == 2, ]
  
  n <- n1 + n2
  if(is.null(Means1)){
    S1 <- cov(as.matrix(X1))
  }
  else{
    S1 <- tcrossprod(t(X1) - Means1) / (n1 - 1)
  }
  if(is.null(Means2)){
    S2 <- cov(as.matrix(X2))
  }
  else{
    S2 <- tcrossprod(t(X2) - Means2) / (n2 - 1)
  }
  W <- (n1 - 1) * S1 + (n2 - 1) * S2
  return(W/(n-2))
}




formDiscrimVector <- function(TrainData, TrainCat, W = NULL, gamma = 1E-5){
  
  # --- Form within-group covariance matrix W if NULL ---
  if(is.null(W)){
    W <- formWithinGroupCov(TrainData = TrainData,
                            TrainCat = TrainCat)
    W <- W + gamma * diag(ncol(TrainData))
  }
  
  TrainX1 <- as.matrix(TrainData[ TrainCat == 1, ])
  TrainX2 <- as.matrix(TrainData[ TrainCat == 2, ])
  
  # --- Form difference of group means vector d ---
  Means1 <- colMeans(TrainX1)
  Means2 <- colMeans(TrainX2)
  
  n1 <- nrow(TrainX1)
  n2 <- nrow(TrainX2)
  n <- n1 + n2
  
  const <- sqrt(n1) * sqrt(n2) / n
  d <- const * (Means1 - Means2)
  
  # --- Form Discriminant Vector ---
  Dvec <- solve(W) %*% d
  return(Dvec)
}


computeMahalanobis <- function(TrainData, TrainCat, TestData, W = NULL, Dvec = NULL, gamma = 1E-5){
  
  #--- Check dimensions of training and test Data---
  if(ncol(TrainData) != ncol(TestData)){
    return(warning("The Training and Test Data have different number
                   of features in computeMahalanobis function.
                   Perhaps previous operations modified data."))
  }
  
  #--- Compute W if NULL value is passed ----
  if(is.null(W)){
    W <- formWithinGroupCov(TrainData = TrainData,
                            TrainCat = TrainCat)
    W <- W + gamma * diag(ncol(W))
  }
  
  # ---- Compute Discriminant Vector if NULL value is passed ----
  if(is.null(Dvec)){
    Dvec <- formDiscrimVector(TrainData = TrainData,
                              TrainCat = TrainCat,
                              W = W,
                              gamma = 1E-5)
  }
  
  # -------- Seperate Data into Groups ------
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  n1 <- nrow(TrainX1)
  n2 <- nrow(TrainX2)
  n <- n1 + n2
  
  #----------Compute Group Means---------------
  Mean1 <- colMeans(as.matrix(TrainX1))
  Mean2 <- colMeans(as.matrix(TrainX2))
  
  #------------Compute Projected Covariance -----------------
  ProjCov <- crossprod(Dvec, W) %*% Dvec #multiply W by Dvec
  
  # ------------Compute Mahalanobis Distances for TestData------------
  Mahala1 <- ((t(t(TestData) - Mean1) %*% Dvec)^2) / as.numeric( ProjCov ) - 2 * log(n1 / n)
  Mahala2 <- ((t(t(TestData) - Mean2) %*% Dvec)^2) / as.numeric( ProjCov ) - 2 * log(n2 / n)
  
  return(list(Mahala1 = Mahala1, Mahala2 = Mahala2))
}



# title Classify
# description A function which implements full LDA on the supplied data.
# param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
# param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
# param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
# param Dvec An optinal (p x 1) vector used as the discriminant vector. Defualt value is NULL, as the function generates its own discriminant vector.
# param W An optinal (p x p) matrix used as the within-group covariance matrix. Defualt value is NULL, as the function generates its own covariance matrix.
# param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
# description Generates linear class predictions for \code{TestData}.
# details Function for full LDA.
Classify <- function(TrainData, TrainCat, TestData, Dvec = NULL, W = NULL, gamma = 1E-5){
  n1 <- nrow(TrainData[TrainCat == 1, ])
  n2 <- nrow(TrainData[TrainCat == 2, ])
  n <- n1 + n2
  
  #---------form covariance matrix ----------
  if(is.null(W)){
    W <- formWithinGroupCov(TrainData = TrainData,
                            TrainCat = TrainCat)
    W <- W + gamma*diag(ncol(W))
  }
  
  #-------- Form discrimiant vector --------------
  if(is.null(Dvec)){
    Dvec <- formDiscrimVector(TrainData = TrainData,
                              TrainCat = TrainCat,
                              W = W,
                              gamma = gamma)
  }
  #-------- Create Labels -----------
  MahalaDistances <- computeMahalanobis(TrainData = TrainData,
                                        TrainCat = TrainCat,
                                        TestData = TestData,
                                        W = W,
                                        Dvec = Dvec)
  
  Mahala1 <- MahalaDistances$Mahala1
  Mahala2 <- MahalaDistances$Mahala2
  
  #compute Mahalanobis distances
  TestPredCat <- as.numeric(Mahala1 < Mahala2)
  TestPredCat[TestPredCat == 0] <- 2
  return(list(Dvec = Dvec, Predictions = TestPredCat))
}


# title compressPredict
# description A function which implements compressed LDA on the supplied data.
# param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
# param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
# param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
# param m1 The number of class 1 compressed samples to be generated. Must be a positive integer.
# param m2 The number of class 2 compressed samples to be generated. Must be a positive integer.
# param s The sparsity level used in compression. Must satify 0 < s < 1.
# param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
# param type A string of characters determining the type of compression matrix used. The accepted values are \code{Rademacher}, \code{Gaussian}, and \code{Count}.
# description Generates the compressed class predictions for \code{TestData}.
# details Function for compressed LDA.
compressPredict <- function(TrainData, TrainCat, TestData, m1, m2, s, gamma = 1E-5, type = "Rademacher"){
  #----- Split Data into Groups ------
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  n1 <- nrow(TrainX1)
  n2 <- nrow(TrainX2)
  
  #--------- compute Group Means -------
  Means1 <- colMeans(TrainX1)
  Means2 <- colMeans(TrainX2)
  
  # -------- Compress Data ----------
  compTrain <- compressData(X1 = TrainX1,
                            X2 = TrainX2,
                            m1 = m1,
                            m2 = m2,
                            s = s,
                            type = type)
  compCat <- c(rep(1, m1), rep(2, m2))
  
  compTrainX1 <- as.matrix(compTrain$Cat1)
  compTrainX2 <- as.matrix(compTrain$Cat2)
  compTrain <- rbind(compTrainX1, compTrainX2)
  
  # -------- Form Compressed Covariance Matrix -----------
  
  Wcomp <- formWithinGroupCov(TrainData = compTrain,
                              TrainCat = compCat)
  
  Wcomp <- Wcomp + gamma*diag(ncol(Wcomp))
  
  
  # ------- Create Labels using Compressed Covariance Matrix --------
  output <- Classify(TrainData = TrainData,
                         TrainCat = TrainCat,
                         TestData = TestData,
                         W = Wcomp,
                         gamma = gamma)
  
  return(output)
}


subsampleClasses <- function(TrainData, TrainCat, m1, m2){
  stopifnot(m1 > 0 & m2 > 0)
  # ------ Split Data into Groups -----------
  TrainX1 <- TrainData[ TrainCat == 1, ]
  TrainX2 <- TrainData[ TrainCat == 2, ]
  n1 <- nrow(TrainX1)
  n2 <- nrow(TrainX2)
  n <- n1 + n2
  
  #-------- Subset Data -----------
  Index1 <- sample(1:n1, size = m1)
  Sub1 <- TrainX1[Index1, ]
  
  Index2 <- sample(1:n2, size = m2)
  Sub2 <- TrainX2[Index2, ]
  
  SubData <- rbind(Sub1, Sub2)
  SubCat <- c(rep(1, m1), rep(2, m2))
  return(list(Data = SubData, Cat = SubCat))
}


# title subsetPredict
# description A function which implements sub-sampled LDA on the supplied data.
# param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
# param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
# param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
# param m1 The number of class 1 sub-samples. Must be a positive integer.
# param m2 The number of class 2 sub-samples. Must be a positive integer.
# param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
# description Generates the sub-sampled linear discriminant vector and class predictions for \code{TestData}.
# details Function for sub-sampled LDA.
subsetPredict <- function(TrainData, TrainCat, TestData, m1, m2, gamma = 1E-5){
  stopifnot(ncol(TrainData) == ncol(TestData))
  stopifnot(m1 > 0 & m2 > 0)
  
  output <- subsampleClasses(TrainData = TrainData,
                             TrainCat = TrainCat,
                             m1 = m1,
                             m2 = m2)
  SubData <- output$Data
  SubCat <- output$Cat
  
  #------ Form Subset Labels ----------
  output <- Classify(TrainData = SubData,
                        TrainCat = SubCat,
                        TestData = TestData,
                        gamma = gamma)
  return(output)
}




# title projectPredict
# description A function which implements projected LDA on the supplied data.
# param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
# param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
# param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
# param Q1 An optional (m1 x n1) sparse matrix used for compressing class 1. The default value is NULL.
# param Q2 An optional (m2 x n2) sparse matrix used for compressing class 2. The default value is NULL.
# param m1 The number of class 1 compressed samples to be generated. Must be a positive integer.
# param m2 The number of class 2 compressed samples to be generated. Must be a positive integer.
# param s The sparsity level used in compression. Must satify 0 < s < 1.
# param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
# param type A string of characters determining the type of compression matrix used. The accepted values are \code{Rademacher}, \code{Gaussian}, and \code{Count}.
# description Generates the compressed linear discriminant vector and projected class predictions for \code{TestData}.
# details Function for projected LDA.
projectPredict <- function(TrainData, TrainCat, TestData, Q1 = NULL, Q2 = NULL, m1, m2, s = 0.01, gamma = 1E-5, type = "Rademacher"){
  #----- Split Data into Groups ------
  TrainX1 <- TrainData[TrainCat == 1, ]
  TrainX2 <- TrainData[TrainCat == 2, ]
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  n <- n1 + n2
  
  # -------- Compress Data ----------
  if(is.null(Q1)){
    compTrain <- compressData(X1 = TrainX1,
                              X2 = TrainX2,
                              m1 = m1,
                              m2 = m2,
                              s = s,
                              type = type)
    
    compTrain1 <- as.matrix(compTrain$Cat1)
  }
  else{
    Means1 <- colMeans(TrainX1)
    TrainX1 <- t(t(TrainX1) - Means1)
    compTrain1 <- as.matrix(Q1 %*% TrainX1)
    m1 <- nrow(compTrain1)
  }
  
  if(is.null(Q2)){
    compTrain <- compressData(X1 = TrainX1,
                              X2 = TrainX2,
                              m1 = m1,
                              m2 = m2,
                              s = s,
                              type = type)
    compTrain2 <- as.matrix(compTrain$Cat2)
  }
  
  else{
    Means2 <- colMeans(TrainX2)
    TrainX2 <- t(t(TrainX2) - Means2)
    compTrain2 <- as.matrix(Q2 %*% TrainX2)
    m2 <- nrow(compTrain2)
  }
  
  #--- Form Compressed Covariance Matrix
  compTrain <- as.matrix(rbind(compTrain1, compTrain2))
  compCat <- c(rep(1, m1), rep(2, m2))
  
  Wcomp <- formWithinGroupCov(TrainData = compTrain,
                              TrainCat = compCat)
  Wcomp <- Wcomp + gamma*diag(ncol(Wcomp))
  
  #----- Form Discriminant Vector----------
  Dvec <- as.numeric(formDiscrimVector(TrainData = TrainData,
                                       TrainCat = TrainCat,
                                       W = Wcomp,
                                       gamma = gamma))
  
  #------- Form Projected Training and Test Data ---------------
  
  ProjTrain <- TrainData %*% Dvec
  ProjTest <- TestData %*% Dvec
  
  
  ProjTrainX1 <- as.matrix(ProjTrain[TrainCat == 1, ])
  ProjTrainX2 <- as.matrix(ProjTrain[TrainCat == 2, ])
  
  #----- Form Covariance matrix of projected Data----------
  Wproj <- formWithinGroupCov(TrainData = ProjTrain,
                              TrainCat = TrainCat)
  
  #----------Compute Group Means---------------
  Mean1 <- colMeans(ProjTrainX1)
  Mean2 <- colMeans(as.matrix(ProjTrainX2))
  
  # ------------Compute Mahalanobis Distances for TestData------------
  Mahala1 <-  ((ProjTest - Mean1)^2) / as.numeric( Wproj ) - 2 * log(n1 / n)
  Mahala2 <- ((ProjTest - Mean2)^2) / as.numeric( Wproj ) - 2 * log(n2 / n)
  
  Labels <- as.numeric(Mahala1 < Mahala2)
  Labels[Labels == 0] <- 2
  return(list(Dvec = Dvec, Predictions = Labels))
}









testParameters <- function(m1,m2,s){
  if(is.null(m1)){
    stop("Number of compressed samples m1 is null.")
  }
  else if(is.null(m2)){
    stop("Number of compressed samples m2 is null.")
  }
  else if(is.null(s)){
    stop("Compression sparsity level s is null.")
  }
  else if(m1 <= 0){
    stop("Number of compressed samples m1 is less than or equal to 0.")
  }
  else if(m2 <= 0){
    stop("Number of compressed samples m2 is less than or equal to 0.")
  }
  else if(s <= 0){
    stop("Compression sparsity level s is less than or equal to 0.")
  }
  else if(m1 != round(m1)){
    stop("Number of compressed samples m1 is not an integer.")
  }
  else if(m2 != round(m2)){
    stop("Number of compressed samples m2 is not an integer.")
  }
  else if(s > 1){
    stop("Compression sparsity level s is greater than 1.")
  }
  else{
    return(TRUE)
  }
}



#' @title Linear Discriminant Analysis (LDA)
#' @description A wrapper function for the various LDA implementations available in this package.
#' @param TrainData A (n x p) numeric matrix without missing values consisting of n training samples each with p features.
#' @param TrainCat A vector of length n consisting of group labels of the n training samples in \code{TrainData}. Must consist of 1s and 2s.
#' @param TestData A (m x p) numeric matrix without missing values consisting of m training samples each with p features. The number of features must equal the number of features in \code{TrainData}.
#' @param Method A string of characters which determines which version of LDA to use. Must be either "Full", "Compressed", "Subsampled", "Projected", or "fastRandomFisher". Default is "Full".
#' @param Mode A string of characters which determines how the reduced sample paramters will be inputted for each method. Must be either "Research", "Interactive", or "Automatic". Default is "Automatic".
#' @param m1 The number of class 1 compressed samples to be generated. Must be a positive integer.
#' @param m2 The number of class 2 compressed samples to be generated. Must be a positive integer.
#' @param m The number of total compressed samples to be generated. Must be a positive integer.
#' @param s The sparsity level used in compression. Must satify 0 < s < 1.
#' @param gamma A numeric value for the stabilization amount gamma * I added to the covariance matrixed used in the LDA decision rule. Default amount is 1E-5. Cannot be negative.
#' @param type A string of characters determining the type of compression matrix used. The accepted values are \code{Rademacher}, \code{Gaussian}, and \code{Count}.
#' @description Generates class predictions for \code{TestData}.
#' @references Lapanowski, Alexander F., and Gaynanova, Irina. ``Compressing large sample data for discriminant analysis'' arXiv preprint arXiv:2005.03858 (2020).
#' @references Ye, Haishan, Yujun Li, Cheng Chen, and Zhihua Zhang. ``Fast Fisher discriminant analysis with randomized algorithms.'' Pattern Recognition 72 (2017): 82-92.
#' @details Function which handles all implementations of LDA. 
#' @examples 
#' TrainData <- LDA_Data$TrainData
#' TrainCat <- LDA_Data$TrainCat
#' TestData <- LDA_Data$TestData
#' plot(TrainData[,2]~TrainData[,1], col = c("blue","orange")[as.factor(TrainCat)])
#' 
#' #----- Full LDA -------
#' LDA(TrainData = TrainData,
#'     TrainCat = TrainCat,
#'     TestData = TestData,
#'     Method = "Full",
#'     gamma = 1E-5)
#'   
#' #----- Compressed LDA -------  
#'  m1 <- 700
#'  m2 <- 300
#'  s <- 0.01
#'  LDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Compressed",
#'      Mode = "Research",
#'      m1 = m1,
#'      m2 = m2,
#'      s = s,
#'      gamma = 1E-5)
#'      
#'  LDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Compressed",
#'      Mode = "Automatic",
#'      gamma = 1E-5)
#'  
#'  #----- Sub-sampled LDA ------
#'  m1 <- 700
#'  m2 <- 300
#'  LDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Subsampled",
#'      Mode = "Research",
#'      m1 = m1,
#'      m2 = m2,
#'      gamma = 1E-5)
#'  
#'   LDA(TrainData = TrainData,
#'      TrainCat = TrainCat,
#'      TestData = TestData,
#'      Method = "Subsampled",
#'      Mode = "Automatic",
#'      gamma = 1E-5)
#'      
#'  #----- Projected LDA ------
#'   m1 <- 700
#'   m2 <- 300
#'   s <- 0.01
#'   LDA(TrainData = TrainData,
#'       TrainCat = TrainCat,
#'       TestData = TestData,
#'       Method = "Projected",
#'       Mode = "Research",
#'       m1 = m1, 
#'       m2 = m2, 
#'       s = s,
#'       gamma = 1E-5)
#'       
#'    LDA(TrainData = TrainData,
#'       TrainCat = TrainCat,
#'       TestData = TestData,
#'       Method = "Projected",
#'       Mode = "Automatic",
#'       gamma = 1E-5)
#'       
#'  #----- Fast Random Fisher ------    
#'   m <- 1000 
#'   s <- 0.01
#'   LDA(TrainData = TrainData,
#'       TrainCat = TrainCat,
#'       TestData = TestData,
#'       Method = "fastRandomFisher",
#'       Mode = "Research",
#'       m = m, 
#'       s = s,
#'       gamma = 1E-5)
#'       
#'    LDA(TrainData = TrainData,
#'       TrainCat = TrainCat,
#'       TestData = TestData,
#'       Method = "fastRandomFisher",
#'       Mode = "Automatic",
#'       gamma = 1E-5)  
#' @return A list containing 
#' \item{Predictions}{(m x 1) Vector of predicted class labels for the data points in \code{TestData}.}
#' \item{Dvec}{(px1) Discriminant vector used to predict the class labels.}  
#' @export
LDA <- function(TrainData, TrainCat, TestData, Method = "Full", Mode = "Automatic", m1 = NULL, m2 = NULL, m = NULL, s = NULL, gamma = 1E-5, type = "Rademacher"){
  if(Mode == "Automatic"){
    m1 <- sum(TrainCat == 1)/10
    m2 <- sum(TrainCat == 2)/10
    m <- length(TrainCat)/10
    s <- 1/sqrt(length(TrainCat))
  }
  
  if(Method == "Full"){
    return(Classify(TrainData = TrainData, 
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
    return(compressPredict(TrainData = TrainData,
                           TrainCat = TrainCat,
                           TestData = TestData,
                           m1 = m1,
                           m2 = m2, 
                           s = s,
                           gamma = gamma,
                           type = type))
  }
  else if(Method == "Projected"){
    if(Mode == "Interactive"){
      m1 <- readline(prompt="Please enter the number of compressed group 1 samples: ")
      m2 <- readline(prompt="Please enter the number of compressed group 2 samples: ")
      s <- readline(prompt="Please enter sparsity level s used in compression: ")
      
      m1 <- as.integer(m1)
      m2 <- as.integer(m2)
      s <- as.numeric(s)
      testParameters(m1 = m1, m2 = m2, s = s)
    }
    
    return(projectPredict(TrainData = TrainData,
                          TrainCat = TrainCat,
                          TestData = TestData,
                          m1 = m1,
                          m2 = m2,
                          s = s,
                          gamma = gamma,
                          type = type))
  }
  else if(Method == "Subsampled"){
    if(Mode == "Interactive"){
      m1 <- readline(prompt="Please enter the number of group 1 sub-samples: ")
      m2 <- readline(prompt="Please enter the number of group 2 sub-samples: ")
     
      m1 <- as.integer(m1)
      m2 <- as.integer(m2)
    }
    return(subsetPredict(TrainData = TrainData,
                         TrainCat = TrainCat,
                         TestData = TestData,
                         m1 = m1,
                         m2 = m2,
                         gamma = gamma))
  }
  else if(Method == "fastRandomFisher"){
    if(Mode == "Interactive"){
      m <- readline(prompt="Please enter the number of total compressed samples m: ")
      s <- readline(prompt="Please enter sparsity level s used in compression: ")
      
      m <- as.integer(m)
      s <- as.numeric(s)
    }
    return(fastRandomFisher(TrainData = TrainData,
                            TrainCat = TrainCat,
                            TestData = TestData,
                            m = m,
                            s = s,
                            gamma = gamma,
                            type = type))
  }
  return("Error: No version of LDA was run.")
}





  