library(DAAG)
library(datasets)
library(MASS)

test_that("Testing QDAclassify equals MASS package QDA", {
  TrainCat <- leafshape17$arch + 1
  leaf17.qda <- MASS::qda(arch ~ logwid+loglen, data = leafshape17)
  MASSclass <- as.numeric(predict(leaf17.qda)$class)

  QDAclass <- QDAclassify(TrainData = data.matrix(leafshape17[, c(5, 7)]),
                           TrainCat = TrainCat,
                           TestData = data.matrix(leafshape17[, c(5, 7)]),
                           gamma = 0)
  expect_equal(QDAclass , MASSclass)
})

test_that("Testing QDAclassify equals MASS package QDA for Iris Data", {
 TrainData <- data.matrix(iris[,-5])
 TrainCat <- as.numeric(iris$Species)
 TrainCat[TrainCat == 3] <- 2
 Data <- data.frame(TrainData, TrainCat)
 IrisMASS <- MASS::qda(TrainCat ~ Sepal.Length + Sepal.Width + Petal.Length + Petal.Width, data = Data)
 MASSclass <- as.numeric(predict(IrisMASS)$class)

  QDAclass <- QDAclassify(TrainData = TrainData,
                          TrainCat = TrainCat,
                          TestData = TrainData,
                          gamma = 0)
  expect_equal(MASSclass, QDAclass)
})

test_that("Testing subsampleQDA is equal to QDA when subsampling entire data",{
  TrainCat <- leafshape17$arch + 1
  TrainData <- data.matrix(leafshape17[, c(5, 7)])
  n1 <- nrow(TrainData[TrainCat == 1, ])
  n2 <- nrow(TrainData[TrainCat == 2, ])

  QDAclass <- QDAclassify(TrainData = TrainData,
                          TrainCat = TrainCat,
                          TestData = TrainData,
                          gamma = 1E-5)

  QDAsubclass <- subsampleQDA(TrainData = TrainData,
                              TrainCat = TrainCat,
                              TestData = TrainData,
                              m1 = n1,
                              m2 = n2,
                              gamma = 1E-5)
  expect_equal(QDAclass , QDAsubclass)
})


test_that("Testing MASS and package QDA are equal for high dimensional normal",{
  p <- 100
  n1 <- 70000
  n2 <- 30000

  Cov1 <- matrix(0, nrow = p, ncol = p)
  Cov2 <- matrix(0, nrow = p, ncol = p)

  Mean1 <- rep(1/100, p)
  Mean2 <- rep(-1/100, p)

  for(i in 1:p){
    for(j in 1:p){
      Cov1[i,j] <- (0.5)^(abs(i-j))
      Cov2[i,j] <- (-0.5)^(abs(i-j))
    }
  }

  TrainX1 <- mvrnorm(n = n1, mu = Mean1, Sigma = Cov1)
  TrainX2 <- mvrnorm(n = n2, mu = Mean2, Sigma = Cov2)
  TrainData <- rbind(TrainX1, TrainX2)

  TestX1 <- mvrnorm(n = n1/10, mu = Mean1, Sigma = Cov1)
  TestX2 <- mvrnorm(n = n2/10, mu = Mean2, Sigma = Cov2)
  TestData <- rbind(TestX1, TestX2)

  TrainCat <- c(rep(1, n1), rep(2, n2))
  TestCat <- c(rep(1, n1/10), rep(2, n2/10))

  # ------------ Build QDA Models --------------
  QDAlabels <- QDAclassify(TrainData = TrainData,
                           TrainCat = TrainCat,
                           TestData = TestData)

  MASSQDA <- MASS::qda(TrainData, as.factor(TrainCat))
  MASSlabels <- predict(MASSQDA, TestData)$class
  expect_equal(as.numeric(QDAlabels) , as.numeric(MASSlabels))
})
