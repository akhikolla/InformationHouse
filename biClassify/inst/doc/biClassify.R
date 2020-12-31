## ---- include = FALSE---------------------------------------------------------
LOCAL <- identical(Sys.getenv("LOCAL"), "true")
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
#knitr::opts_chunk$set(purl = CRAN)

## -----------------------------------------------------------------------------
library(biClassify)

## -----------------------------------------------------------------------------
data(LDA_Data)

## ---- echo = FALSE, fig.height=4, fig.width=5, fig.align = "center"-----------
plot(LDA_Data$TrainData[,2]~LDA_Data$TrainData[,1],
     col = c("orange","blue")[LDA_Data$TrainCat],
     pch = c("1","2")[LDA_Data$TrainCat],
     xlab = "Feature 1", 
     ylab = "Feature 2",
     main = "Scatter Plot of LDA Training Data")

## ---- eval = FALSE------------------------------------------------------------
#  > test_pred <- LDA(TrainData = LDA_Data$TrainData,
#                     TrainCat = LDA_Data$TrainCat,
#                     TestData = LDA_Data$TestData,
#                     Method = "Compressed")$Predictions
#  
#  > mean(test_pred != LDA_Data$TestCat)
#  [1] 0

## -----------------------------------------------------------------------------
data(QDA_Data)

## ---- echo = FALSE, fig.height=4, fig.width=5, fig.align = "center"-----------
plot(QDA_Data$TrainData[,2]~QDA_Data$TrainData[,1],
     col = c("orange","blue")[QDA_Data$TrainCat],
     pch = c("1","2")[QDA_Data$TrainCat],
     xlab = "Feature 1", 
     ylab = "Feature 2",
     main = "Scatter Plot of QDA Training Data")

## ---- eval = FALSE------------------------------------------------------------
#  > test_pred <- QDA(TrainData = QDA_Data$TrainData,
#                     TrainCat = QDA_Data$TrainCat,
#                     TestData = QDA_Data$TestData,
#                     Method = "Compressed")
#  
#  > mean(test_pred != QDA_Data$TestCat)
#  [1] 0

## -----------------------------------------------------------------------------
data(KOS_Data)

## ----echo = FALSE, fig.height=4, fig.width=8, fig.align = "center"------------
par(mfrow = c(1,2))
plot(KOS_Data$TrainData[,2]~KOS_Data$TrainData[,1], col = c("orange","blue")[KOS_Data$TrainCat],
     pch = c("1","2")[KOS_Data$TrainCat],
     xlab = "Feature 1",
     ylab = "Feature 2",
     main = "True Features")
plot(KOS_Data$TrainData[,4]~KOS_Data$TrainData[,3], col = c("orange","blue")[KOS_Data$TrainCat],
     pch = c("1","2")[KOS_Data$TrainCat],
     xlab = "Feature 3",
     ylab = "Feature 4",
     main = "Noise Features")
par(mfrow = c(1,1))

## ---- eval = FALSE------------------------------------------------------------
#  > output <- KOS(TrainData = KOS_Data$TrainData,
#                  TrainCat = KOS_Data$TrainCat,
#                  TestData = KOS_Data$TestData)
#  > output$Weight
#  [1] 1 1 0 0
#  
#  > mean(output$Predictions != KOS_Data$TestCat)
#  [1] 0
#  
#  > summary(output$Dvec)
#         V1
#   Min.   :-0.03002
#   1st Qu.:-0.01953
#   Median :-0.01445
#   Mean   : 0.00000
#   3rd Qu.: 0.03788
#   Max.   : 0.05799

## -----------------------------------------------------------------------------
TrainData <- LDA_Data$TrainData
TrainCat <- LDA_Data$TrainCat
TestData <- LDA_Data$TestData
TestCat <- LDA_Data$TestCat

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData)$Predictions
table(test_pred)
mean(test_pred != TestCat)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "Compressed", Mode = "Automatic")$Predictions
table(test_pred)
mean(test_pred != TestCat)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, Method = "Compressed")$Predictions
table(test_pred)
mean(test_pred != TestCat)

## ---- eval=FALSE--------------------------------------------------------------
#  output <- LDA(TrainData, TrainCat, TestData,
#                Method = "Compressed", Mode = "Interactive")$Predictions
#  "Please enter the number m1 of group 1 compression samples: "700
#  "Please enter the number m2 of group 2 compression samples: "300
#  "Please enter sparsity level s used in compression: "0.01

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "Compressed", Mode = "Research", 
                 m1 = 700, m2 = 300, s = 0.01)$Predictions

table(test_pred)
mean(test_pred != TestCat)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "Subsampled", Mode = "Automatic")$Predictions
table(test_pred)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "Subsampled")$Predictions
table(test_pred)

## ---- eval=FALSE--------------------------------------------------------------
#  test_pred <- LDA(TrainData, TrainCat, TestData,
#                   Method = "Subsampled", Mode = "Interactive")$Predictions
#  "Please enter the number m1 of group 1 sub-samples: "700
#  "Please enter the number m2 of group 2 sub-samples: "300

## -----------------------------------------------------------------------------
output <- LDA(TrainData, TrainCat, TestData, 
              Method = "Subsampled", Mode = "Research", 
              m1 = 700, m2 = 300)$Predictions

table(output)
mean(output != TestCat)

## -----------------------------------------------------------------------------
output <- LDA(TrainData, TrainCat, TestData, 
              Method = "Projected", Mode = "Automatic")$Predictions
table(output)
mean(output != TestCat)

## -----------------------------------------------------------------------------
output <- LDA(TrainData, TrainCat, TestData, 
              Method = "Projected")$Predictions
table(output)
mean(output != TestCat)

## ---- eval=FALSE--------------------------------------------------------------
#  output <- LDA(TrainData, TrainCat, TestData,
#                Method = "Projected", Mode = "Interactive")$Predictions
#  "Please enter the number m1 of group 1 compression samples: "700
#  "Please enter the number m2 of group 2 compression samples: "300
#  "Please enter sparsity level s used in compression: "0.01
#  

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "Projected", Mode = "Research", 
                 m1 = 700, m2 = 300, s = 0.01)$Predictions

table(test_pred)
mean(output != TestCat)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "fastRandomFisher", Mode = "Automatic")$Predictions
table(test_pred)
mean(test_pred != TestCat)

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "fastRandomFisher")$Predictions
table(test_pred)
mean(test_pred != TestCat)

## ---- eval=FALSE--------------------------------------------------------------
#  output <- LDA(TrainData, TrainCat, TestData,
#                Method = "fastRandomFisher", Mode = "Interactive")$Predictions
#  "Please enter the number m of total compressed samples: "1000
#  "Please enter sparsity level s used in compression: "0.01

## -----------------------------------------------------------------------------
test_pred <- LDA(TrainData, TrainCat, TestData, 
                 Method = "fastRandomFisher", Mode = "Research", 
                 m = 1000, s = 0.01)$Predictions

table(test_pred)
mean(test_pred != TestCat)

## -----------------------------------------------------------------------------
TrainData <- QDA_Data$TrainData
TrainCat <- QDA_Data$TrainCat
TestData <- QDA_Data$TestData
TestCat <- QDA_Data$TestCat

## -----------------------------------------------------------------------------
Predictions <- QDA(TrainData, TrainCat, TestData, Method = "Full")
table(Predictions)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Compressed", Mode = "Automatic")
table(output)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Compressed")
table(output)

## ---- eval=FALSE--------------------------------------------------------------
#  output <- QDA(TrainData, TrainCat, TestData, Method = "Compressed", Mode = "Interactive")
#  "Please enter the number m1 of group 1 compression samples: "700
#  "Please enter the number m2 of group 2 compression samples: "300
#  "Please enter sparsity level s used in compression: "0.01
#  
#  table(output)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Compressed", 
              Mode = "Research", m1 = 700, m2 = 300, s = 0.01)

summary(output)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Subsampled", Mode = "Automatic")
table(output)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Subsampled")
summary(output)

## ---- eval=FALSE--------------------------------------------------------------
#  output <- QDA(TrainData, TrainCat, TestData, Method = "Subsampled", Mode = "Interactive")
#  "Please enter the number m1 of group 1 sub-samples: "700
#  "Please enter the number m2 of group 2 sub-samples: "300
#  
#  summary(output)

## -----------------------------------------------------------------------------
output <- QDA(TrainData, TrainCat, TestData, Method = "Subsampled", 
              Mode = "Research", m1 = 700, m2 = 300)

summary(output)

## -----------------------------------------------------------------------------
TrainData <- KOS_Data$TrainData
TrainCat <- KOS_Data$TrainCat
TestData <- KOS_Data$TestData
TestCat <- KOS_Data$TestCat

## ---- eval = FALSE------------------------------------------------------------
#  > SelectParams(TrainData, TrainCat)
#  
#  $Sigma
#  [1] 0.7390306
#  
#  $Gamma
#  [1] 0.137591
#  
#  $Lambda
#  [1] 0.0401767

## ---- eval = FALSE------------------------------------------------------------
#  > SelectParams(TrainData, TrainCat, Sigma = 1, Gamma = 0.1)
#  
#  $Sigma
#  [1] 1
#  
#  $Gamma
#  [1] 0.1
#  
#  $Lambda
#  [1] 0.06186337

## ---- eval = FALSE------------------------------------------------------------
#  SelectParams(TrainData, TrainCat, Gamma = 0.1)
#  
#  Error in SelectParams(TrainData, TrainCat, Gamma = 0.1) :
#  Hierarchical order of parameters violated.
#  Please specify Sigma before Gamma, and both Sigma and Gamma before Lambda.

## ---- eval = FALSE------------------------------------------------------------
#  Sigma <- 1.325386
#  Gamma <- 0.07531579
#  Lambda <- 0.002855275
#  
#  > output <- KOS(TestData, TrainData, TrainCat, Sigma = Sigma,
#                  Gamma = Gamma, Lambda = Lambda)
#  
#  > output$Weight
#  [1] 1 1 0 0
#  
#  > table(output$Predictions)
#   1  2
#  26 68
#  
#  > summary(output$Dvec)
#         V1
#   Min.   :-0.05860
#   1st Qu.:-0.03711
#   Median :-0.02539
#   Mean   : 0.00000
#   3rd Qu.: 0.06983
#   Max.   : 0.10192

