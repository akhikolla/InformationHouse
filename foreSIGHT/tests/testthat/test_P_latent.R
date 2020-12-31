context("P_latent model check")


test_that("P_latent output matches saved reference", {
  
  randomVector <- c(0.5349710, 0.5463762, 0.8611103, 0.6066115, 0.6701731)
  parSigma <- rep_len(5, length.out = 5)
  parMu <- rep_len(-3, length.out = 5)
  parLambda <- rep_len(1.2, length.out = 5)
  parAlpha <- rep_len(0.8, length.out = 5)
  
  parTS <- list(sigma = parSigma, mu = parMu, lambda = parLambda, alpha = parAlpha)
  
  file_name = "../P_latent_output1.rds"
  expect_equal_to_reference(P_latent(parTS, randomVector), file = file_name)
  
})

test_that("P_latent_master output matches saved reference", {
  
  modelInfo <- list()
  modelInfo$simVar="P"
  modelInfo$simPriority=1
  modelInfo$nperiod=1
  modelInfo$fixedPars=NA
  modelInfo$ncycle=NA
  modelInfo$npars=modelInfo$nperiod*4
  modelInfo$parNam=c("alpha", "sigma", "mu", "lambda")
  modelInfo$minBound=c(0, 0.001, -15, 1)
  modelInfo$maxBound=c(0.999, 10, 0, 2)
  
  modelTag <- "P-ann-latent"
  
  datInd <- list()
  datInd$ndays <- 5
  
  test_modelEnv <- new.env(parent = emptyenv())
  test_modelEnv$modelInfo <- modelInfo
  test_modelEnv$modelTag <- modelTag
  test_modelEnv$datInd <- datInd
  
  parS <- c(0.8, 5, -3, 1.2)
  randomVector <- c(0.5349710, 0.5463762, 0.8611103, 0.6066115, 0.6701731)
  
  file_name = "../P_latent_output1.rds"
  expect_equal_to_reference(P_latent_master(parS=parS,
                                            modelEnv = test_modelEnv,
                                            randomVector = randomVector), file = file_name)
  
})


