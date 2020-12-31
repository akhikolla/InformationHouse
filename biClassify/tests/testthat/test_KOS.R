
#--- Testing Kernel Matrix Functionality ---
test_that("KernelMat function equals Gaussian kernel matrix", {
 Data <- KOS_Data$TrainData
 n <- nrow(Data)
 Sigma <- 1
 K <- KernelMat(TrainData = Data, Sigma = Sigma)
 
 K_test <- matrix(0, nrow = n, ncol = n)
 for(i in 1:n){
   for(j in 1:n){
     K_test[i,j] <- exp(- sum((Data[i,] - Data[j,])^2)/Sigma)
   }
 }
 expect_equal(K, K_test)
})

test_that("Weighted Kernel Matrix equals kernel matrix when weights are 1", {
  Data <- KOS_Data$TrainData
  n <- nrow(Data)
  sigma <- 1
  K <- KernelMat(TrainData = Data, Sigma = sigma)
  Kw <- KwMat(TrainData = Data, w = rep(1, ncol(Data)), Sigma = sigma)
  expect_equal(K, Kw)
})

test_that("Weighted Kernel Matrix equals all 1s when weights are 0", {
  Data <- KOS_Data$TrainData
  n <- nrow(Data)
  sigma <- 1
  Kw <- KwMat(TrainData = Data, w = rep(0, ncol(Data)), Sigma = sigma)
  ones <- matrix(1, nrow = n, ncol = n)
  expect_equal(ones, Kw)
})


test_that("KernelCPP and Kernel function are the same.",{
  TrainData <- KOS_Data$TrainData
  x <- KOS_Data$TestData[1,]
  Sigma <- 1
  KernCPP<- KernelCPP(x = x, TrainData = TrainData, Sigma = Sigma)
  Kern <- Kernel(x = x, TrainData = TrainData, Sigma = Sigma)
  expect_equal(mean(abs(KernCPP - Kern)) < 1E-10 , TRUE)
})

#--- Testing the Optimal Scoring Function ---

test_that("Optimal Scores are +- 1 for equal class sizes", {
  TrainCat <- c(rep(1, 10), rep(2,10))
  scores <- as.numeric(OptScores(TrainCat))
  test_scores <- c(1,-1)
  expect_equal(scores, test_scores)
})

# --- Testing GetProjections ---
test_that("GetProjections produces a projection for every test sample",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  TestCat <- KOS_Data$TestCat
  
  Projections <- GetProjections(TestData = TestData, 
                                  TrainData = TrainData, 
                                  TrainCat = TrainCat,
                                  Sigma = 1,
                                  Gamma = 0.1)
  
  expect_equal(length(Projections), nrow(TestData))
})

# --- Testing FormQB ---
test_that("FormQB returns correct dimensions",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  TestCat <- KOS_Data$TestCat
  
  p <- ncol(TrainData)
    
  #Generate One-hot encoding matrix and optimal scores
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  Kw <- KwMat(TrainData, rep(1,p), 1)
  Dvec <- SolveKOSCPP(YTheta, Kw, 0.1)
  
  QBoutput <- FormQB(TrainData = TrainData, 
                     Dvec = Dvec, 
                     YTheta = YTheta, 
                     Kw = Kw,
                     w = rep(1, ncol(TrainData)), 
                     Sigma = 1, 
                     Gamma = 0.1)
  testQ <- QBoutput$Q
  testB <- QBoutput$B
  
  expect_equal(c(nrow(testQ), length(testB)), rep(p, 2))
})

# --- Testing SparseKernOptScore ---

test_that("SparseKernOptScore equals kernel ridge regression when Lambda = 0",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  
  SparseKOS_output <- SparseKernOptScore(TrainData = TrainData,
                                         TrainCat = TrainCat,
                                         Sigma = 1,
                                         Gamma = 0.1,
                                         Lambda = 0)
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  Kw <- KwMat(TrainData, rep(1,p), 1)
  Kw <- scale(Kw, center = TRUE, scale = FALSE)
  Kw <- scale(t(Kw), center = TRUE, scale = FALSE)  
  
  Dvec <- SolveKOSCPP(YTheta, Kw, 0.1)
  
  expect_equal(Dvec, SparseKOS_output$Dvec)
})

test_that("SparseKernOptScore outputs list with correct dimensions",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  n <- nrow(TrainData)
  
  Sigma <- 1
  Gamma <- 0.1
  Lambda <- 0.01
  
  output <- SparseKernOptScore(TrainData = TrainData,
                              TrainCat = TrainCat,
                              Sigma = Sigma,
                              Gamma = Gamma,
                              Lambda = Lambda)
  expect_equal(c(length(output$Dvec),length(output$Weights)), c(n,p))
})

# --- Test Select Ridge ---
test_that("Checking that the SelectRidge Parameter Gamma > 0 ",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  p <- ncol(TrainData)
  n <- nrow(TrainData)
  
  Gamma <- SelectRidge(TrainData = TrainData,
                       TrainCat = TrainCat,
                       Sigma = 1)

  expect_equal(Gamma > 0, TRUE)
})



# --- Test SelectParams ---
#test_that("Checking that all parameters are positive", {
#  TrainData <- KOS_Data$TrainData
#  TrainCat <- KOS_Data$TrainCat
#  output <- SelectParams(TrainData = TrainData,
#                         TrainCat = TrainCat)
  
#  test_sigma <- output$Sigma
#  test_ridge <- output$Gamma
#  test_lambda <- output$Lambda
  
#  expect_equal(all(test_sigma > 0, test_ridge > 0, test_lambda > 0), TRUE)
#})

# --- Testing Projection Functions Give Equal Output ---
test_that("Rcpp GetProjectionsCPP and R GetProjections functions are equal",{
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  
  n <- nrow(TrainData)
  p <- ncol(TrainData)
  Y <- IndicatMat(TrainCat)$Categorical
  Theta <- OptScores(TrainCat)
  YTheta <- Y %*% Theta
  
  Sigma <- 1
  Gamma <- 0.09853699 #result of SelectRidge
  
  K <- KernelMat(TrainData = TrainData, Sigma = Sigma)
  K <- scale(K, center = TRUE, scale = FALSE)
  K <- scale(t(K), center = TRUE, scale = FALSE)
  
  Dvec <- SolveKOSCPP(YTheta, K, Gamma)
 
  R_Projs <- GetProjections(TrainData = TrainData,
                           TrainCat = TrainCat,
                           TestData = TestData,
                           Dvec = Dvec,
                           w = rep(1, ncol(TrainData)),
                           Kw = K,
                           Sigma = Sigma, 
                           Gamma = Gamma)
  
  CPP_Projs <- GetProjectionsCPP(TrainData = TrainData,
                                 TrainCat = TrainCat,
                                 TestData = TestData,
                                 Dvec = Dvec,
                                 w = rep(1, ncol(TrainData)),
                                 Kw = K,
                                 Sigma = Sigma,
                                 Gamma = Gamma)
  
  expect_equal(sum(abs(R_Projs - CPP_Projs)) < 1E-10, TRUE)
})


#--- Testing SolveKOS --- 

test_that("SolveKOS and SolveKOSCPP are equal", {
  
  TrainData <- KOS_Data$TrainData
  TrainCat <- KOS_Data$TrainCat
  TestData <- KOS_Data$TestData
  
  
  Y<-IndicatMat(TrainCat)$Categorical
  
  #--- Get Everything Initialized ---
  n <- nrow(TrainData)
  p <- ncol(TrainData)
  error <- 1  #Initialize Error value
  niter <- 1
  Sigma <- 1
  Gamma <- 0.1
  
  #Form optimal scores and weighted kernel matrix
  Opt_Scores <- OptScores(TrainCat)
  YTheta <- Y %*% Opt_Scores
  
  Kw <- KwMat(TrainData = TrainData,
              w = rep(1, p), 
              Sigma = Sigma)
  #-- Double-center Kw --
  Kw <- scale(Kw, center = TRUE, scale = FALSE)
  Kw <- scale(t(Kw), center = TRUE, scale = FALSE)
  
  Dvec <- SolveKOS(YTheta = YTheta, 
                       K = Kw, 
                       Gamma = Gamma) 

    
  DvecCPP <- SolveKOSCPP(YTheta, Kw, Gamma)
  expect_equal(sum(abs(Dvec - DvecCPP)) < 1E-7, TRUE)
})