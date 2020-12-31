context("sd")

test_that("Testing for binary matrices", {
  W <- rbm(20,20)
  mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  msd <- motif_sd(W)
  expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
  expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")

  W <- rbm(5,5)
  mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  msd <- motif_sd_minors(W)
  expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
  expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")

  for (i in 1:5) {
    W <- rbm(3,3, 0.2) # make it sparse
    if (any(W != 0)) {
      mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
      msd <- motif_sd_minors(W)
      expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
      expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")
    }
  }
})

test_that("Testing for matrices with equal weights", {
  W <- 2 * rbm(20,20)
  msd <- motif_sd(W)
  mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
  expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")

  W <- 2 * rbm(5,5)
  msd <- motif_sd_minors(W)
  mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
  expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")

  W <- 3.75 * rbm(20,20)
  msd <- motif_sd(W)
  mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
  expect(all(msd[which(mc > 0)] == 0), failure_message = "failed")
  expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")
})

test_that("Testing matrix with single link", {
  W <- matrix(0, 3, 3)
  W[1,1] <- 1
  msd <- motif_sd(W)
  expect_equal(msd[1], 0)
  expect(all(is.na(msd[2:17])), failure_message = "failed")
})

test_that("Testing simple version of M2", {
  W <- matrix(c(1,2), nrow = 1)
  msd <- motif_sd(W)
  expect_equal(msd[1], pop_sd(c(1,2)))
  expect_equal(msd[2], 0)
  expect(all(is.na(msd[3:17])), failure_message = "failed")
})


test_that("Testing simple version of M3", {
  W <- matrix(c(2.34,7.21), nrow = 2)
  msd <- motif_sd(W)
  expect_equal(msd[1], pop_sd(c(2.34,7.21)))
  expect(is.na(msd[2]), failure_message = "failed")
  expect_equal(msd[3], 0)
  expect(all(is.na(msd[4:17])), failure_message = "failed")

  msd <- motif_sd_minors(W)
  expect_equal(msd[1], pop_sd(c(2.34,7.21)))
  expect(is.na(msd[2]), failure_message = "failed")
  expect_equal(msd[3], 0)
  expect(all(is.na(msd[4:17])), failure_message = "failed")
})

test_that("Testing simple version of M4", {
  W <- matrix(c(5,7, 2), nrow = 3)
  msd <- motif_sd(W)
  expect_equal(msd[1], pop_sd(c(5,2,7)))
  expect(is.na(msd[2]), failure_message = "failed")
  m3_weights <- c(mean(c(5,7)), mean(c(5,2)), mean(c(7,2)))
  expect_equal(msd[3], pop_sd(m3_weights))
  expect_equal(msd[4], 0)
  expect(all(is.na(msd[5:17])), failure_message = "failed")

  msd <- motif_sd_minors(W)
  expect_equal(msd[1], pop_sd(c(5,2,7)))
  expect(is.na(msd[2]), failure_message = "failed")
  m3_weights <- c(mean(c(5,7)), mean(c(5,2)), mean(c(7,2)))
  expect_equal(msd[3], pop_sd(m3_weights))
  expect_equal(msd[4], 0)
  expect(all(is.na(msd[5:17])), failure_message = "failed")
})

# test_that("Testing simple version of M5", {
#   W <- matrix(c(5,7, 0, 2), nrow = 2)
#   msd <- motif_sd(W)
#   expect_equal(msd[1], pop_sd(c(5,2,7)))
#   expect(all(msd[c(2,3,5)] == 0))
#   expect(all(is.na(msd[c(4,6:17)])))
# })
test_that("Testing simple version of M5_2",{
  W <- matrix(c(5,7, 0, 2), nrow = 2)
  msd <- motif_sd_minors(W)
  expect_equal(msd[1], pop_sd(c(5,2,7)))
  expect(all(msd[c(2,3,5)] == 0), failure_message = "failed")
  expect(all(is.na(msd[c(4,6:17)])), failure_message = "failed")
})

test_that("Testing simple version of M9", {
  W <- matrix(c(5, 0, 7, 0, 2, 3), nrow = 3, byrow = 3)
  msd <- motif_sd(W)
  expect_equal(msd[1], pop_sd(c(5,2,7, 3)))
  expect_equal(msd[2], 0)
  m3_weights <- c(mean(c(5,2)), mean(c(7,2)), mean(c(5,7)))
  expect_equal(msd[3], pop_sd(m3_weights))
  expect_equal(msd[4], 0)
  m5_weights <- c(mean(c(5,2,3)), mean(c(7,2,3)))
  expect_equal(msd[5], pop_sd(m5_weights))
  expect_equal(msd[9], 0)
  expect(all(is.na(msd[6:8])), failure_message = "failed")
  expect(all(is.na(msd[10:17])), failure_message = "failed")

  W <- matrix(c(5, 0, 7, 0, 2, 3), nrow = 3, byrow = TRUE)
  msd <- motif_sd_minors(W)
  expect_equal(msd[1], pop_sd(c(5,2,7, 3)))
  expect_equal(msd[2], 0)
  expect_equal(msd[3], pop_sd(m3_weights))
  expect_equal(msd[4], 0)
  expect_equal(msd[5], pop_sd(m5_weights))
  expect_equal(msd[9], 0)
  expect(all(is.na(msd[6:8])), failure_message = "failed")
  expect(all(is.na(msd[10:17])), failure_message = "failed")
})


test_that("Compare two versions", {
  for (i in 1:5) {
    W <- rwm(5,5, 0.6)
    msd_minors <- motif_sd_minors(W)
    msd <- motif_sd(W)
    expect_equal(msd_minors, msd)
  }
  for (i in 1:5) {
    W <- rwm(5,5, 0.7)
    msd_minors <- motif_sd_minors(W)
    msd <- motif_sd(W)
    expect_equal(msd_minors, msd)
  }
  for (i in 1:5) {
    W <- rwm(5,5, 0.8)
    msd_minors <- motif_sd_minors(W)
    msd <- motif_sd(W)
    expect_equal(msd_minors, msd)
  }
})

test_that("Test: If mcount is zero, sd should be NA", {
  for (i in 1:10) {
    W <- rwm(10,10, 0.1)
    if(sum(W) != 0) {
      mc <- mcount(W, six_node = FALSE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
      msd <- motif_sd(W)
      expect(all(is.na(msd[which(mc == 0)])), failure_message = "failed")
    }
  }
})
