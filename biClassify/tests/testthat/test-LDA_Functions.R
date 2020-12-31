library(datasets)
library(DAAG)
library(mvtnorm)


test_that("Testing dimensions of compression matrix", {
  n <- 100000
  m <- 1000
  Q <- createSketchMatrix(n = n, m = m, s = 0.01)
  n_test <- ncol(Q)
  m_test <- nrow(Q)
  expect_equal(n_test, n)
  expect_equal(m_test, m)
})

test_that("Testing discrimiant vector if Null cov matrix is passed into function",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  W_test <- formWithinGroupCov( TrainData = TrainData,
                                TrainCat = TrainCat)
  DiscrimVec1 <- formDiscrimVector(TrainData = TrainData,
                                   TrainCat = TrainCat,
                                   W = W_test)
  DiscrimVec2 <- formDiscrimVector(TrainData = TrainData,
                                   TrainCat = TrainCat,
                                   gamma = 0)
  expect_equal(DiscrimVec1, DiscrimVec2)
})


test_that("Testing classification compared to MASS LDA on leaf data", {
  TrainCat <- leafshape17$arch + 1
  leaf17.LDA <- MASS::lda(arch ~ logwid+loglen, data = leafshape17)
  MASSclass <- as.numeric(predict(leaf17.LDA)$class)
  
  compLDAclass <- Classify(TrainData = data.matrix(leafshape17[, c(5, 7)]),
                           TrainCat = TrainCat,
                           TestData = data.matrix(leafshape17[, c(5, 7)]),
                           gamma = 0)
  expect_equal(compLDAclass$Predictions , MASSclass)
})

test_that("Testing classification compared to MASS LDA on large normal data.", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70000
  n2 <- 30000
  n <- n1 + n2
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  TrainData <- as.data.frame(cbind(unlist(TrainData), TrainCat))
  TestData <- as.data.frame(unlist(TestData))
  colnames(TrainData) <- c("1","2","3","4","5","6","7","8","9","10","Class")
  colnames(TestData) <- c("1","2","3","4","5","6","7","8","9","10")
  
  #--- Compute MASS Labels ---
  normal.LDA <- MASS::lda(Class ~ . , data = TrainData)
  MASSclass <- as.numeric(predict(object = normal.LDA, newdata = TestData)$class)
  
  #--- Compute compressLDA labels ----
  TrainData <- as.matrix(TrainData[,c(1:10)])
  Labels <- Classify(TrainData = TrainData,
                     TrainCat = TrainCat,
                     TestData = TestData,
                     gamma = 0)
  
  expect_equal(Labels$Predictions , MASSclass)
})




test_that("Check that Compressed Predict does not return NA",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  compLabels <- compressPredict(TrainData = TrainData,
                                TrainCat = TrainCat,
                                TestData = TrainData,
                                m1 = 10,
                                m2 = 10,
                                s = 1/2)$Predictions
  expect_equal(any(is.na(compLabels)), FALSE)
})



test_that("Check that Subset Predict does not return NA",{
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  subLabels <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = 10,
                             m2 = 10)$Predictions
  expect_equal(any(is.na(subLabels)), FALSE)
})



test_that("Check that subsetPredict equals Classify if subset entire data Leaf data", {
  TrainCat <- leafshape17$arch + 1
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  subLabels <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = n1,
                             m2 = n2)
  Labels <- Classify(TrainData = TrainData,
                     TrainCat = TrainCat,
                     TestData = TrainData)
  expect_equal(Labels , subLabels)
})

test_that("Testing if subsetPredict equals Classify if subset is entire data for large data set.", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70000
  n2 <- 30000
  n <- n1 + n2
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  subLabels <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = n1,
                             m2 = n2)
  Labels <- Classify(TrainData = TrainData,
                     TrainCat = TrainCat,
                     TestData = TrainData)
  expect_equal(Labels , subLabels)
})

test_that("Testing if subsetPredict equals Classify if subset is entire data for small data set.", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70
  n2 <- 30
  n <- n1 + n2
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  subLabels <- subsetPredict(TrainData = TrainData,
                             TrainCat = TrainCat,
                             TestData = TrainData,
                             m1 = n1,
                             m2 = n2)
  Labels <- Classify(TrainData = TrainData,
                     TrainCat = TrainCat,
                     TestData = TrainData)
  expect_equal(Labels$Predictions , subLabels$Predictions)
})



test_that("checking that projectPredict equals LDA when compression matrix is identity",{
  TrainCat <- as.numeric(leafshape17$arch + 1)
  TrainData = data.matrix(leafshape17[, c(5, 7)])
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  Q1 <- diag(n1)
  Q2 <- diag(n2)
  ProjClass <- projectPredict(TrainData = TrainData,
                              TrainCat = TrainCat,
                              TestData = TrainData,
                              Q1 = Q1,
                              Q2 = Q2,
                              m1 = n1,
                              m2 = n2,
                              s = 0.01,
                              gamma = 0)
  leaf17.LDA <- MASS::lda(arch ~ logwid+loglen, data = leafshape17)
  MASSclass <- as.numeric(predict(leaf17.LDA)$class)
  expect_equal(MASSclass, ProjClass$Predictions)
})


test_that("Testing formWithGroupCov matrix equals default when means are NULL", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  for(i in 1:p){
    for(j in 1:p){
      Sig[i,j] <- (0.5)^(abs(i-j))
    }
  }
  
  n1 <- 70000
  n2 <- 30000
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  
  Means1 <- colMeans(TrainX1)
  Means2 <- colMeans(TrainX2)
  
  TrainCat <- c(rep(1, n1), rep(2, n2))
  NullCovMat <- formWithinGroupCov(TrainData, TrainCat)
  MeansCovMat <- formWithinGroupCov(TrainData, TrainCat, Means1, Means2)
  expect_equal(NullCovMat, MeansCovMat)
})


test_that("Mahalanobis Distance equals Eucliden distance when covariance is identity and equal class sizes", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 50000
  n2 <- 50000
  n <- n1 + n2
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  Dvec_test <- rep(1/sqrt(p), p)
  W_test <- diag(p)
  
  Mahalanobis <- computeMahalanobis(TrainData = TrainData,
                                    TrainCat = TrainCat,
                                    TestData = TestData,
                                    W = W_test,
                                    Dvec = Dvec_test)
  
  TrainMean1 <- colMeans(TrainData[TrainCat == 1, ])
  TestDistances <- rowSums(t(t(TestData) - TrainMean1))/sqrt(p)
  TestDistances <- TestDistances^2 - 2* log(n1 / n)
  
  expect_equal(sum(abs(Mahalanobis$Mahala1 - TestDistances)), 0, tolerance=1e-8)
})

test_that("Testing if if formDiscrimVector gives scaling of MASS discriminant Vector for large normal data", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70000
  n2 <- 30000
  n <- n1 + n2
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  Dvec <- formDiscrimVector(TrainData, TrainCat, gamma = 0)
  MASSDvec <- normal.LDA <- MASS::lda(TrainCat ~ TrainData)$scaling
  
  expect_equal(as.numeric(Dvec/as.numeric(MASSDvec)) , rep(Dvec[1]/as.numeric(MASSDvec[1]), p ),tolerance = 1e-10)
})



#----- Testing subSampleClasses function  ------
test_that("Testing if Classify equals MASS LDA on sub-sampled data using the
subsampleClasses function to generate data. Large normal data.", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70000
  n2 <- 30000
  n <- n1 + n2
  p1 <- n1 / n
  p2 <- 1 - p1
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  TestData <- as.data.frame(unlist(TestData))
  
  #--- Generate list of subsample amounts that Classify and MASS will be compared on.
  #--- The outputs must be exactly equal across the entire sequence to pass test. ----
  subSeq <- c(50, 100, 500, 100, 500, 1000, 2000, 5000)
  
  PassTest <- rep(FALSE, length(subSeq)) #initialize tests to failure. Each trial
  # must actively pass test and change entry.
  
  for(i in 1:length(PassTest)){
    m1 <- floor(p1 * subSeq[i]) #Generate subsample amounts
    m2 <- floor(p2 * subSeq[i])
    
    subOutput <- subsampleClasses(TrainData, TrainCat, m1, m2)
    subData <- subOutput$Data
    subCat <- subOutput$Cat
    subData <- as.data.frame(cbind(unlist(subData), subCat))
    colnames(subData) <- c("1","2","3","4","5","6","7","8","9","10","Class")
    colnames(TestData) <- c("1","2","3","4","5","6","7","8","9","10")
    
    #--- Compute MASS Labels ---
    normal.LDA <- MASS::lda(Class ~ . , data = subData)
    MASSclass <- as.numeric(predict(object = normal.LDA, newdata = TestData)$class)
    
    #--- Compute compressLDA labels ----
    subData <- as.matrix(subData[,c(1:10)])
    Labels <- Classify(TrainData = subData,
                       TrainCat = subCat,
                       TestData = TestData,
                       gamma = 0)$Predictions
    
    PassTest[i] <- all(Labels == MASSclass)
    
  }
  expect_equal(PassTest, rep(TRUE, length(subSeq)))
})


#----- Testing subSampleClasses function  ------
test_that("Testing if subsammpleClasses outputs correct number of rows for each class.", {
  p <- 10
  mu1 <- rep(2, p)
  mu2 <- rep(-2, p)
  Sig <- matrix(0,nrow = p, ncol = p)
  Sig <- diag(p)
  n1 <- 70000
  n2 <- 30000
  n <- n1 + n2
  p1 <- n1 / n
  p2 <- 1 - p1
  
  TrainX1 <- rmvnorm(n = n1, mean = mu1, sigma = Sig)
  TrainX2 <- rmvnorm(n = n2, mean = mu2, sigma = Sig)
  TrainData <- rbind(TrainX1, TrainX2)
  TrainCat <- c(rep(1, n1), rep(2, n2))
  
  TestX1 <- rmvnorm(n = n1/10, mean = mu1, sigma = Sig)
  TestX2 <- rmvnorm(n = n2/10, mean = mu2, sigma = Sig)
  TestData <- rbind(TestX1, TestX2)
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))
  
  TestData <- as.data.frame(unlist(TestData))
  
  #--- Generate list of subsample amounts that Classify and MASS will be compared on.
  #--- The outputs must be exactly equal across the entire sequence to pass test. ----
  subSeq <- c(50, 100, 500, 100, 500, 1000, 2000, 5000)
  
  RowTest <- rep(FALSE, length(subSeq)) #initialize tests to failure. Each trial
  # must actively pass test and change entry.
  CatTest1 <- rep(FALSE, length(subSeq))
  CatTest2 <- rep(FALSE, length(subSeq))
  
  PassTest <- rep(FALSE, length(subSeq)) # Finaly test comparison. Will be
  # conjunction of all tests
  
  for(i in 1:length(PassTest)){
    m1 <- floor(p1 * subSeq[i]) #Generate subsample amounts
    m2 <- floor(p2 * subSeq[i])
    
    subOutput <- subsampleClasses(TrainData, TrainCat, m1, m2)
    subData <- subOutput$Data
    subCat <- subOutput$Cat
    RowTest[i] <- nrow(subData) == (m1 + m2)
    
    CatTest1[i] <- sum(subCat == 1) == m1
    CatTest2[i] <- sum(subCat == 2) == m2
    
    PassTest[i] <- (RowTest[i] & CatTest1[i] & CatTest2[i])
  }
  expect_equal(PassTest, rep(TRUE, length(subSeq)))
})





