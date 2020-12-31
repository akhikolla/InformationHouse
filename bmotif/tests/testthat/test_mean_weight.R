context("mean weight")

test_that("Testing for binary matrices", {
  for (i in 1:10) {
    W <- rbm(20,20)
    mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
    mw <- mean_weight(W, six_node = TRUE, mc = mc)
    expect(all(mw[which(mc > 0)] == 1), failure_message = "failed")
    expect(all(is.na(mw[which(mc == 0)])), failure_message = "failed")
  }
})

test_that("Testing for matrices with equal weights", {
  W <- 3 * rbm(20,20)
  mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  mw <- mean_weight(W, mc, TRUE)
  expect(all(mw[which(mc > 0)] == 3), failure_message = "failed")
  expect(all(is.na(mw[which(mc == 0)])), failure_message = "failed")

  W <- 0.75 * rbm(20,20)
  mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  mw <- mean_weight(W, mc, TRUE)
  expect(all(mw[which(mc > 0)] == 0.75), failure_message = "failed")
  expect(all(is.na(mw[which(mc == 0)])), failure_message = "failed")
})

test_that("Testing matrix with single link", {
  W <- matrix(0, 3, 3)
  W[1,1] <- 1
  mw <- mean_weight(W)
  expect_equal(mw[1], 1, failure_message = "failed")
  expect(all(is.na(mw[2:17])), failure_message = "failed")
})

test_that("Comparing with motif extraction algorithm", {
  for (i in 1:10) {
    W <- rwm(6,6)
    mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
    mw <- mean_weight(W, mc, TRUE)
    mw_minors <- mean_weight_minors(W, TRUE)
    expect_equal(mw, mw_minors)
  }
})

test_that("Testing simple version of M4", {
  W <- matrix(c(1,2,3), nrow = 3)
  mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  mw <- mean_weight(W, mc, TRUE)
  expect_equal(mw[1], 2)
  expect(is.na(mw[2]), failure_message = "failed")
  expect_equal(mw[3], 2)
  expect_equal(mw[4], 2)
  expect(all(is.na(mw[5:44])), failure_message = "failed")
})

test_that("Test: If mcount is zero, mean should be NA", {
  for (i in 1:10) {
    W <- rwm(10,10, 0.1)
    if(sum(W) != 0) {
      mc <- mcount (W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
      mw <- mean_weight(W, mc, TRUE)
      expect_identical(which(is.na(mw)), which(mc == 0))
    }
  }
})
