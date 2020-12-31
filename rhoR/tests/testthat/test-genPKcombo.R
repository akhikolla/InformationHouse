test_that("verify precision kappa combinations", {
  kappa_dist <- seq(0.01, 0.9, length = 100)
  prec_dist <- 0.21
  
  pk <- rhoR:::genPKcombo(kappa_dist, NULL, prec_dist, NULL, 0.2)
  testthat::expect_lte(pk$kappa, 0.02)
  testthat::expect_equal(pk$precision, 0.21)
  
  
  kappa_dist <- seq(0.01, 0.9, length = 100)
  prec_dist <- 0.21
  
  br = 0.8
  kappaMin = 0.1
  kappaMax = 0.1
  precisionMin = 0.2
  precisionMax = 0.2
  # kps.r = rhoR::generateKPs(numNeeded = 100, br, kappaMin, kappaMax, precisionMin, precisionMax)
  # kps.df = do.call("rbind", kps.r)
  # # kps.df[,1] = 0.5
  # rhoR:::find_valid_pk(kappaDistribution = as.numeric(kps.df[,2]), kappaProbability = rep(1, nrow(kps.df)), precisionDistribution = as.numeric(kps.df[,1]), precisionProbability = rep(1, nrow(kps.df)), baserate = br)
  # rhoR:::find_valid_pk(c(0.1, 0.1), c(1, 1), 0.2, 1, 0.6)
  
  ## Thorws error in endless loop
  pk <- rhoR:::find_valid_pk(c(rep(0.6, 1000), 0.4), rep(1, 1001), 0.5, 1, 0.2)
  testthat::expect_equal(pk[1], 0.5)
  testthat::expect_equal(pk[2], 0.4)
  
  ## Force find new
  pk <- rhoR:::find_valid_pk(0.2, 1, c(rep(0.1, 10000), 0.9), c(rep(1, 10000), 0.1), 0.2)
  testthat::expect_equal(pk[1], 0.9)
  testthat::expect_equal(pk[2], 0.2)
})

test_that("create KP values", {
  kps_with_bell <- generateKPs(10, 0.2, 0.7, 0.9, 0.6, 0.9, distributionType = 'BELL')
  kp_kappas <- sapply(kps_with_bell, `[[`, "kappa")
  kp_precs <- sapply(kps_with_bell, `[[`, "precision")
  
  expect_false(any(is.na(c(kp_kappas, kp_precs))))
  expect_true(all(kp_kappas < 0.9 & kp_kappas > 0.7))
})

test_that("force KP error", {
  expect_error(generateKPs(10, 0.2, 0.7, 0.7, 0.1, 0.1), regexp = "Could not.*?ranges of kappa and precision")
})

test_that("verify p boots", {
  rand_ks <- runif(1000, min = 0.75, max = 0.8)
  actual_k <- 0.78
  
  p_r <- rhoR:::getBootPvalue(rand_ks, actual_k)
  p_c <- rhoR:::getBootPvalue_c(rand_ks, actual_k)
  expect_equal(p_r, p_c)
})

# test_that("calc rho tests", {
#   observedKappa = 0.67
#   fullSetBaserate = 0.2
#   handSetLength = 20
#   handSetBaserateInflation = 0.2
#   
#   kps <- generateKPs(100, 0.2, 0.6, 0.8, 0.6, 0.9,)
#   rhoR:::calcRho(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation, KPs = kps)
#   # Rs = replicate(10, rhoR:::calcRho(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation))
#   # Cs = replicate(10, rhoR:::calcRho_c(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation))
# })