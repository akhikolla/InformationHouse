library(hetGP)
context("IMSE")


test_that("IMSE",{
  
  set.seed(1)
  nvar <- 2
  
  ## Branin redefined in [0,1]^2
  branin <- function(x){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    x1 <- x[,1] * 15 - 5
    x2 <- x[,2] * 15
    (x2 - 5/(4 * pi^2) * (x1^2) + 5/pi * x1 - 6)^2 + 10 * (1 - 1/(8 * pi)) * cos(x1) + 10
  }
  
  ## Noise field via standard deviation
  noiseFun <- function(x){
    if(is.null(nrow(x)))
      x <- matrix(x, nrow = 1)
    return(1/20*(3*(2 + 2*sin(x[,1]*pi)*cos(x[,2]*3*pi) + 5*rowSums(x^2))))
  }
  
  ## data generating function combining mean and noise fields
  ftest <- function(x){
    return(branin(x) + rnorm(nrow(x), mean = 0, sd = noiseFun(x)))
  }
  
  ## Predictive locations for comparisons
  ntest <- 1e5
  Xgrid <- matrix(runif(ntest*2), ntest, 2)
  
  ## Unique (randomly chosen) design locations
  n <- 20
  Xu <- matrix(runif(n * 2), n)
  
  ## Select replication sites randomly
  X <- Xu[sample(1:n, 10*n, replace = TRUE),]
  
  ## obtain training data response at design locations X
  Z <- ftest(X)
  
  ## Formating of data for model creation (find replicated observations) 
  prdata <- find_reps(X, Z, rescale = FALSE, normalize = FALSE)
  
  for(covtype in c("Gaussian", "Matern5_2", "Matern3_2")){
    for(trend in c(NA, 0)){
      
      if(is.na(trend)) trend <- NULL
      ## Model fitting
      model <- mleHetGP(X = list(X0 = prdata$X0, Z0 = prdata$Z0, mult = prdata$mult), Z = prdata$Z,
                        lower = rep(0.01, nvar), upper = rep(1, nvar), known = list(beta0 = trend),
                        covtype = covtype)
      IMSE <- IMSPE(model)
      
      # Finite sum estimation:
      preds <- predict(model, Xgrid)
      IMSE_est <- mean(preds$sd2)
      
      expect_equal(IMSE, IMSE_est, tol = 1e-2)
    }
  }
  
})