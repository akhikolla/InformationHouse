library(testthat)
context("Testing utility methods within package")

test_that("p value", {
  dist <- seq(60, 80) / 100
  dist <- dist[dist != 0.7]

  expect_equal(getBootPvalue(dist, 0.6), 1)
  expect_equal(getBootPvalue(dist, 0.7), 0.5)
  
  expect_equal(getBootPvalue_c(rep(0.8, 10), 0.7), 1)
})

test_that("brpk values", {
  expect_false(checkBRPKcombo(0.2, 0.2, 0.8))
  expect_true(checkBRPKcombo(0.2, 0.8, 0.8))
})

# test_that("p combos", {
#   p_min <- 0.80;
#   p_max <- 0.85;
#   
#   combo <- rhoR:::genPcombo(0.8, p_min, p_max, 0.5)
#   expect_true(p_min < combo && combo < p_max)
#   
#   expect_match(tryCatch(genPcombo(0.9, p_min, p_max, 0.5), error = function(e) { e$message }), regexp = "nested too deeply")
# })

test_that("create random sets", {
  len <- 200
  br <- 0.2
  set <- rhoR:::createRandomSet(len, br, 0.4, 0.65, 0, 1)
  
  expect_gte(sum(set[,1]), len * br)
  expect_gte(round(kappa(set), 2), 0.38)
  expect_lte(round(kappa(set), 2), 0.67)
  
  ct <- rhoR:::createRandomSet(len, br, 0.4, 0.65, 0, 1, type = "ct")
  expect_gte(sum(ct[1, ]), len * br)
  expect_gte(round(kappa(ct), 2), 0.38)
  expect_lte(round(kappa(ct), 2), 0.67)
  
  ct <- rhoR:::random_contingency_table(len, br, 0.5, 0.75, 0, 1)
  expect_gte(sum(ct[1, ]), len * br)
  expect_gte(round(kappa(ct), 2), 0.48)
  expect_lte(round(kappa(ct), 2), 0.77)
})

test_that("rho min", {
  br <- 0.2
  rho_min <- rhoMin(br)
  expect_lte(rho_min, 50)
  
  expect_error(rhoMin(br, inc = 12.5), regexp = "Inc value must")
  expect_error(rhoMin(br, alpha = 2), regexp = "Alpha must")
  
  expect_output(rhoR:::rhoMin(br, printInc = TRUE), regexp = "10\\s0\\.")
})

test_that("rho given kappa errors", {
  expect_error(rhoK(2, OcSBaserate = 0.5, testSetLength = 10), regexp = "between 0 and 1")
  expect_error(rhoK(-2, OcSBaserate = 0.5, testSetLength = 10), regexp = "between 0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = -1), regexp = "OcSBaserate.*?0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 1.5), regexp = "testSetLength.*?positive integer")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, OcSLength = 2.5), regexp = "OcSLength.*?positive integer")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, replicates = 1.5), regexp = "replicates.*?positive integer")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSKappaThreshold = -1), regexp = "ScSKappaTh.*?positive")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSKappaMin = -1), regexp = "ScSKappaMin.*?positive")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSKappaMin = 0.9, ScSKappaThreshold = 0.8), regexp = "ScSKappaThresh.*?greater.*?Min")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSPrecisionMin = -2), regexp = "ScSPrecisionMin.*?0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSPrecisionMin =  2), regexp = "ScSPrecisionMin.*?0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSPrecisionMax = -2), regexp = "ScSPrecisionMax.*?0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSPrecisionMax =  2), regexp = "ScSPrecisionMax.*?0 and 1")
  expect_error(rhoK(0.9, OcSBaserate = 0.2, testSetLength = 200, ScSPrecisionMin = 0.6, ScSPrecisionMax = 0.5), regexp = "ScSPrecisionMax.*?greater.*?Min")
})

test_that("create from PR", {
  prec <- 0.8
  rec <- 0.6
  hsl <- 40
  br <- 0.2
  
  set <- rhoR:::prset(prec, rec, hsl, br)
  ct  <-   contingency_table(prec, rec, hsl, br)

  set_tp <- sum(set[, 1] == 1 & set[, 2] == 1)
  set_tn <- sum(set[, 1] == 0 & set[, 2] == 0)
  set_fn <- sum(set[, 1] == 1 & set[, 2] == 0)
  set_fp <- sum(set[, 1] == 0 & set[, 2] == 1)
  
  expect_equal(set_tp, ct[1, 1])
  expect_equal(set_tn, ct[2, 2])
  expect_equal(set_fn, ct[1, 2])
  expect_equal(set_fp, ct[2, 1])
})

test_that("rhoK on set and ct", { 
  len = 1000
  br = 0.2
  set <- rhoR:::createRandomSet(len, br, 0.9, 1, 0.8, 1)

  set_k <- kappa(set)
  ocs_br <- 0.2
  tsl <- 80

  set_k_std <- replicate(10, rhoK(set_k, OcSBaserate = ocs_br, testSetLength = tsl))
  set_k_c   <- replicate(10, rhoK(set_k, OcSBaserate = ocs_br, testSetLength = tsl, method = "c"))
  testthat::expect_equivalent(mean(set_k_std), mean(set_k_c), tolerance = 0.01)
})

test_that("Generate KPs", {
  num <- 10000
  br  <- 0.2
  kmin <- 0.60
  kmax <- 0.65
  pmin <- 0.8
  pmax <- 0.85
  
  kps_list <- generateKPs(num, br, kmin, kmax, pmin, pmax)
  
  kps_r <- matrix(unlist(kps_list), byrow = TRUE, ncol = 2, dimnames = list(NULL, names(kps_list[[1]])))
  
  kps_c <- generate_kp_list(num, br, kmin, kmax, pmin, pmax)
  
  r_means <- colMeans(kps_r)
  c_means <- colMeans(kps_c)
  testthat::expect_true(rhoR:::all_equal(r_means[1], c_means[1], tolerance = 0.001))
  testthat::expect_true(rhoR:::all_equal(r_means[2], c_means[2], tolerance = 0.001))
})

test_that("CT kappa comparisons", {
  prec <- 0.8
  rec <- 0.6
  hsl <- 40
  br <- 0.2

  ct <- contingency_table(prec, rec, hsl, br)
  k_r <- kappaCT(ct)
  k_c <- kappa_ct(ct)
  
  testthat::expect_equivalent(k_r, k_c)
})

test_that("compare calcRhos", {
  num <- 100
  br  <- 0.2
  kmin <- 0.85
  kmax <- 0.95
  pmin <- 0.8
  pmax <- 1
  
  # kps_list <- rhoR:::generateKPs(num, br, kmin, kmax, pmin, pmax)
  # kps_r <- matrix(unlist(kps_list), byrow = TRUE, ncol = 2, dimnames = list(NULL, names(kps_list[[1]])))
  
  # col_1 <- round(kps_r[, 1], 4)
  # col_2 <- round(kps_r[, 2], 4)
  # col_1_c <- paste0("c(", paste(shQuote(col_1), collapse = ", "), ")")
  # col_2_c <- paste0("c(", paste(shQuote(col_2), collapse = ", "), ")")
  # col_1_n <- as.numeric(eval(parse(text = col_1_c)))
  # col_2_n <- as.numeric(eval(parse(text = col_2_c)))
  
  col_1_n <- c('0.9477', '0.8846', '0.995', '0.9396', '0.9107', '0.9608', '0.9244', '0.9663', '0.9119', '0.9362', '0.951', '0.9921', '0.9482', '0.9848', '0.8578', '0.8077', '0.8813', '0.8965', '0.9972', '0.9852', '0.8374', '0.8802', '0.9167', '0.9315', '0.9851', '0.8017', '0.9858', '0.9597', '0.9911', '0.9765', '0.9331', '0.9192', '0.9422', '0.9097', '0.9745', '0.9111', '0.9888', '0.9791', '0.9715', '0.9704', '0.87', '0.9081', '0.9594', '0.9945', '0.8935', '0.8569', '0.9383', '0.8428', '0.9888', '0.8261', '0.9201', '0.8785', '0.9254', '0.9884', '0.8937', '0.9242', '0.9903', '0.9901', '0.9476', '0.8804', '0.9043', '0.9409', '0.8611', '0.9427', '0.8889', '0.9737', '0.8799', '0.9152', '0.9924', '0.8264', '0.9396', '0.888', '0.9308', '0.947', '0.8757', '0.9344', '0.9977', '0.9203', '0.9397', '0.993', '0.993', '0.9421', '0.9979', '0.9192', '0.8453', '0.89', '0.9293', '0.8432', '0.8989', '0.88', '0.929', '0.9512', '0.9832', '0.8459', '0.9716', '0.8141', '0.9671', '0.9729', '0.9396', '0.8969')
  col_2_n <- c('0.9112', '0.8881', '0.8552', '0.9354', '0.8877', '0.872', '0.8874', '0.9495', '0.9168', '0.9312', '0.9483', '0.9083', '0.9246', '0.8948', '0.8749', '0.8502', '0.8654', '0.9071', '0.8625', '0.9197', '0.8561', '0.86', '0.919', '0.931', '0.9106', '0.8557', '0.8745', '0.9066', '0.9257', '0.944', '0.9105', '0.8832', '0.9038', '0.9081', '0.938', '0.8537', '0.9294', '0.9156', '0.878', '0.9453', '0.8798', '0.9372', '0.8897', '0.9069', '0.8939', '0.8985', '0.9465', '0.8858', '0.8757', '0.8584', '0.9423', '0.857', '0.8533', '0.9303', '0.8974', '0.8847', '0.9486', '0.9432', '0.9316', '0.9058', '0.9326', '0.8833', '0.8728', '0.8769', '0.9217', '0.936', '0.8633', '0.8762', '0.9494', '0.8729', '0.9094', '0.8722', '0.9401', '0.9283', '0.8552', '0.9378', '0.8546', '0.8658', '0.8843', '0.9097', '0.8675', '0.8783', '0.8855', '0.8707', '0.8832', '0.8795', '0.883', '0.8515', '0.8611', '0.9091', '0.8809', '0.9363', '0.9112', '0.8765', '0.9125', '0.8524', '0.8681', '0.8538', '0.9349', '0.8765')
  kps_r <- matrix(as.numeric(c(col_1_n, col_2_n)), ncol = 2)
  kps_list <- apply(kps_r, 1, function(r) list(precision = r[1], kappa = r[2]))

  observedKappa = 0.91
  fullSetBaserate = 0.2
  handSetLength = 100
  handSetBaserateInflation = 0.2
  
  calc_rs <- replicate(50, rhoR:::calcRho(  observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation, KPs = kps_list))
  calc_cs <- replicate(50, rhoR:::calcRho_c(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation, KPs = kps_r))
  expect_true(rhoR:::all_equal(mean(calc_rs), mean(calc_cs), tolerance = 0.1))
  
  calc_rs <- replicate(50, rhoR:::calcRho(  observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation, KPs = kps_list, replicates = 50))
  calc_cs <- replicate(50, rhoR:::calcRho_c(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation, KPs = kps_r, replicates = 50))
  expect_true(rhoR:::all_equal(mean(calc_rs), mean(calc_cs), tolerance = 0.3))
  
  calc_rs <- replicate(10, rhoR:::calcRho(  observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation))
  calc_cs <- replicate(10, rhoR:::calcRho_c(observedKappa, fullSetBaserate, handSetLength, handSetBaserateInflation))
  expect_true(rhoR:::all_equal(mean(calc_rs), mean(calc_cs), tolerance = 0.1))
})

test_that("kps from c", {
  num <- 50000
  br  <- 0.2
  kmin <- 0
  kmax <- 1
  pmin <- 0
  pmax <- 1

  kps_c <- generate_kp_list(num, br, kmin, kmax, pmin, pmax, distributionType = 1)
  expect_equal(mean(kps_c[,2]), 0.9, tolerance = 0.04)
})

test_that("ct equals set", {
  ct_as_set <- as.code.set(contingencyTable)
  expect_equal(ct_as_set[,1], codeSet[,1])
  expect_equal(ct_as_set[,2], codeSet[,2])
})
test_that("set equals ct", {
  set_as_ct <- as.contingency.table(codeSet)
  expect_equivalent(unclass(set_as_ct), contingencyTable)
})

test_that("verify rating set printing", {
  # tmp <- tempfile()
  # on.exit(unlink(tmp), add = TRUE)
  
  ct_as_set <- as.code.set(contingencyTable)
  # sink(tmp)
  # print(ct_as_set)

  # output <- readLines(tmp)
  # unlink(tmp)
  # sink(file = NULL)
  # testthat::expect_true(grepl(x = output, pattern = "\\s+\\[,1\\]\\s\\[,2\\]")[1])
  testthat::expect_output(print(ct_as_set), regexp = "first_rater second_rater" )
})

test_that("test dollar completion", {
  set_as_ct <- as.contingency.table(codeSet)
  expect_null(rhoR:::`$.rating.set`(set_as_ct, ""))
  expect_equal(rhoR:::`$.rating.set`(set_as_ct, "kappa"), 0.625)
  expect_equal(rhoR:::`$.rating.set`(set_as_ct, "second_negative"), 35)
  expect_equal(rhoR:::`$.rating.set`(set_as_ct, "first_positive"), 4)
  
  ct_as_set <- as.code.set(contingencyTable)
  expect_null(rhoR:::`$.rating.set`(ct_as_set, ""))
  expect_equal(rhoR:::`$.rating.set`(ct_as_set, "kappa"), 0.625)
  expect_equal(sum(rhoR:::`$.rating.set`(ct_as_set, "second_rater") == 0), 35)
  expect_true(all(rhoR:::.codegen_dollar(ct_as_set) %in% rhoR:::additional_attributes))
  
  expect_true(all(`.DollarNames.rating.set`(ct_as_set) %in% c("first_rater", "second_rater", "baserate", "kappa")))
  expect_true(all(`.DollarNames.rating.set`(set_as_ct) %in% c("first_positive", "first_negative",
                                                                   "second_positive", "second_negative",
                                                                   "agreement", "baserate", "kappa")))
})

# test_that("package loading", {
#   testthat::expect_error(rhoR:::.onLoad("", "rhoR"))
# })
