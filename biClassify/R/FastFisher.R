# --- Functions for implementing Fast Random Fisher Discriminant Analysis ---

formIndicatorMatrix <- function(TrainCat){
  E <- matrix(0 , nrow = length(TrainCat), ncol = 2)
  E[TrainCat == 1 , 1] <- 1
  E[TrainCat == 2, 2] <- 1
  return(E)
}


formProbMat <- function(TrainCat){
  E <- formIndicatorMatrix(TrainCat)
  return(crossprod(E,E))
}


SquareRootInvPi <- function(TrainCat){
  Pi <- formProbMat(TrainCat = TrainCat)
  return(sqrt(solve(Pi)))
}



formGmat <- function(TrainData, TrainCat, m, s, Method = "Full", gamma = 1E-5, type = "Rademacher"){
  #----- form Matrices Used -----------
  sqrtInvPi <- SquareRootInvPi(TrainCat = TrainCat)
  E <- formIndicatorMatrix(TrainCat = TrainCat)
  if(Method == "Full"){
    S <- createSketchMatrix(n = nrow(TrainData),
                            m = m,
                            s = s,
                            type = type)
  }
  else{
    n1 <- nrow(TrainData[ TrainCat == 1, ])
    n2 <- nrow(TrainData[ TrainCat == 2, ])
    p1 <- n1/(n1+n2)
    p2 <- 1 - p1
    m1 <- floor(p1 * m)
    m2 <- floor(p2 * m)

    S1 <- createSketchMatrix(n = n1, m = m1, s = s, type = type)
    S2 <- createSketchMatrix(n = n2, m = m2, s = s, type = type)
    S <- Matrix::bdiag(S1, S2)
  }

  #----- column-center TrainData ----------
  TrainData <- t(t(TrainData) - colMeans(TrainData))

  #----- Multiply Sketch to matrices-----------
  SketchData <- S %*% TrainData

  SketchResponse <- (S %*% E) %*% sqrtInvPi

  #----- Form G -----
  Mat1 <- crossprod(as.matrix(SketchData), as.matrix(SketchResponse))
  Mat2 <- crossprod(as.matrix(SketchData), as.matrix(SketchData)) + gamma * diag(ncol(TrainData))

  return(solve(Mat2)%*% Mat1)
}



fastRandomFisher <- function(TrainData, TrainCat, TestData, G = NULL, m, s = 0.01, gamma = 1E-5, type = "Rademacher"){
  #------- Form Projection Matrix and Projected Data -----------
  if(is.null(G)){
    G <- formGmat(TrainData = TrainData,
                  TrainCat = TrainCat,
                  m = m,
                  s = s,
                  gamma = gamma, 
                  type = type)
  }
  # --- isolate first column as discriminant vector ---
  G <- as.matrix(G)
  G <- G[,1]
  
  
  n1 <- sum(TrainCat == 1)
  n2 <- sum(TrainCat == 2)
  n <- n1 + n2

  ProjTrain <- as.matrix(TrainData %*% G)
  ProjTest <- as.matrix(TestData %*% G)
  ProjTrainX1 <- as.matrix(ProjTrain[TrainCat == 1, ])
  ProjTrainX2 <- as.matrix(ProjTrain[TrainCat == 2, ])

  #----- Form Covariance matrix of projected Data----------
  Wproj <- formWithinGroupCov(TrainData = ProjTrain,
                              TrainCat = TrainCat)

  #----------Compute Group Means---------------
  Mean1 <- colMeans(ProjTrainX1)
  Mean2 <- colMeans(as.matrix(ProjTrainX2))

  # ------------Compute Mahalanobis Distances for Projected Test Data------------
  ProjTest1 <- t(t(ProjTest) - Mean1)
  Mahala1 <-  t(solve(Wproj) %*% t(ProjTest1))
  Mahala1 <- rowSums(ProjTest1 * Mahala1) - 2 * log(n1 / n)

  ProjTest2 <- t(t(ProjTest) - Mean2)
  Mahala2 <- t(solve(Wproj) %*% t(ProjTest2))
  Mahala2 <- rowSums(ProjTest2 * Mahala2) - 2 * log(n2 / n)

  Labels <- as.numeric(Mahala1 < Mahala2)
  Labels[Labels == 0] <- 2
  return(list(Dvec = G, Predictions = Labels))
}
