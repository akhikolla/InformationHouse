
# library(testthat)

###############################################################################@
############################# item_analysis ####################################
###############################################################################@

test_that("Test item_analysis function", {
  n_item <- 10 # sample(8:12, 1)
  n_theta <- 100 # sample(100:200, 1)
  ip <- itempool(b = rnorm(n_item))
  theta <- setNames(rnorm(n_theta), paste0("Ex-", 1:n_theta))
  resp <- sim_resp(ip = ip, theta = theta)
  # Add some missing responses
  resp[sample(1:length(resp), round(length(resp)*.1))] <- NA
  # Run item analysis:
  ia <- item_analysis(resp = resp)
  # Calculate the p-value
  i <- sample(1:n_item, 1)
  pval <- sum(resp[, i], na.rm = TRUE)/sum(!is.na(resp[, i]))
  expect_equal(ia$pval[i], pval)

})

###############################################################################@
############################# distractor_analysis ##############################
###############################################################################@
test_that("Test distractor_analysis function", {
  n_item <- 10 # sample(8:12, 1)
  n_theta <- 50 # sample(100:200, 1)
  raw_resp <- matrix(sample(LETTERS[1:4], n_item * n_theta, replace = TRUE),
                     nrow = n_theta, ncol = n_item,
                     dimnames = list(paste0("Examinee-", 1:n_theta),
                                     paste0("Item-", 1:n_item)))
  # Add some missing responses
  raw_resp[sample(1:length(raw_resp), round(length(raw_resp)*.1))] <- NA
  key <- sample(LETTERS[1:4], n_item, replace = TRUE)
  da <- distractor_analysis(raw_resp = raw_resp, key = key,
                            suppress_output = TRUE)
  expect_equal(key, da$biserial$key) # Check the key column

  # Check whether responses scored correctly
  for (i in 1:n_item) { # i is for item number
    j <- sample(1:n_theta, 1) # j is for examinee number
    expect_equal(1 * (raw_resp[j, i] == key[i]), da$scores[j, i])
  }

  for (item_col in sample(1:n_item, 3)) {
    choice <- sample(key, 1)

    # Check the number of occurences
    expect_equal(sum(!is.na(raw_resp[, item_col])), da$prop$n[item_col])
    # Check choice proportion
    expect_equal(sum(raw_resp[, item_col] == choice, na.rm = TRUE) /
                   sum(!is.na(raw_resp[, item_col])),
                 da$prop[[choice]][item_col])
    # Check biserial correlation
    expect_equal(da$biserial[item_col, choice],
                 biserial(score = 1 * (raw_resp[, item_col] == choice),
                          total_score =  rowSums(da$scores, na.rm = TRUE)))
  }

  # -------------------------------------------------------------------------- #
  # Check distractor analysis with adjustment
  da <- distractor_analysis(raw_resp = raw_resp, key = key, adjusted = TRUE,
                            suppress_output = TRUE)

  for (item_col in sample(1:n_item, 3)) {
    choice <- sample(colnames(da$biserial)[-1], 1)
    ts <- rowSums(da$scores[, -item_col], na.rm = TRUE)
    score <- (raw_resp[, item_col] == choice) * 1
    expect_equal(biserial(score = score, total_score = ts),
                 da$biserial[item_col, choice])
  }
})



###############################################################################@
############################# biserial #########################################
###############################################################################@
test_that("Test biserial function", {
  # -------------------------------------------------------------------------- #
  # Check biserial correlation
  # The example is from page 72 of: Lewis-Beck, Michael S., Alan E. Bryman, and
  # Tim Futing Liao. 2004. The SAGE Encyclopedia of Social Science Research
  # Methods. Thousand Oaks, CA: SAGE Publications.
  score <- c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 0, 1,
             1, 1, 1, 0, 1, 1, 1, 1)
  total_score <- c(8, 12, 6, 12, 8, 8, 8, 11, 13, 4, 14, 13, 10, 9, 8, 33, 28,
                   29, 30, 29, 28, 33, 32, 32, 33, 34, 35, 34, 38, 37)
  expect_equal(biserial(score = score, total_score = total_score), 0.8449,
               tol = 0.001)

  # -------------------------------------------------------------------------- #
  # Check biserial and point-biserial correlation
  # The example is from Salkind, Rasmussen (2007) Encyclopedia of measurement
  # and statistics, pages 94-97
  score <- c(rep(0, 16), rep(1, 22))
  total_score <- c(87, 90, 94, 94, 97, 103, 103, 104, 106, 108, 109, 109, 109,
                   112, 119, 132, 100, 103, 103, 106, 112, 113, 114, 114, 118,
                   119, 120, 120, 124, 133, 135, 135, 136, 141, 155, 157, 159,
                   162)
  expect_equal(biserial(score = score, total_score = total_score), 0.703,
               tol = 0.01)
  expect_equal(biserial(score = score, total_score = total_score,
                        method = "point-biserial"), 0.557, tol = 0.001)

  # -------------------------------------------------------------------------- #
  # Check biserial and point-biserial correlation, biserial_brogden
  # The example is from Kotz (2005) Encyclopedia of Statistical Sciences,
  # pages 577-580
  # This exact same example is also on SAS:
  #   https://support.sas.com/kb/24/991.html
  score <- c(1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
             0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1)
  total_score <- c(14.8, 13.8, 12.4, 10.1, 7.1, 6.1, 5.8, 4.6, 4.3, 3.5, 3.3,
                   3.2, 3, 2.8, 2.8, 2.5, 2.4, 2.3, 2.1, 1.7, 1.7, 1.5, 1.3,
                   1.3, 1.2, 1.2, 1.1, 0.8, 0.7, 0.6, 0.5, 0.2, 0.2, 0.1)
  expect_equal(biserial(score = score, total_score = total_score), 0.45369,
               tol = 0.0001)
  expect_equal(biserial(score = score, total_score = total_score,
                        method = "point-biserial"), 0.36149, tol = 0.0001)
  expect_equal(biserial(score = score, total_score = total_score,
                        method = "brogden"), 0.50, tol = 0.01)
  expect_equal(biserial(score = score, total_score = total_score,
                        method = "rank"), 0.4359, tol = 0.0001)

  # biserial_cpp(score = score, total_score = total_score,
  #              method = "clemans-lord")
  # biserial_cpp(score = score, total_score = total_score,
  #              method = "brogden")
  # biserial_cpp(score = score, total_score = total_score,
  #              method = "rank")
  # biserial_cpp(score = score, total_score = total_score,
  #              method = "default")
  # n <- 1e6
  # score = sample(0:1, n, replace = TRUE)
  # total_score = sample(10:31, n, replace = TRUE)
  # microbenchmark::microbenchmark(
  #   R = biserial(score = score, total_score = total_score),
  #   Cpp = biserial(score = score, total_score = total_score)
  #   )



})



