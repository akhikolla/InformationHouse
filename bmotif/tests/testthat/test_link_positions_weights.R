context("Weighted link positions")

test_that("Check binary matrices", {
  for (i in 1:10) {
    W <- rbm(20,20)
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    expect(all(lw == lp), failure_message = "failed")
  }
})

test_that("Check matrices with equal weights", {
  for (i in 1:5) {
    W <- 3 * rbm(20,20)
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    expect(all(lw == 3 * lp), failure_message = "failed")
  }
  for (i in 1:5) {
    W <- 0.5 * rbm(20,20)
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    expect(all(lw == 0.5 * lp), failure_message = "failed")
  }
})


test_that("Check matrices with one weighted link", {
  for (i in 1:5) {
    W <- rbm(20,20)
    W[1,1] <- 10
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    expect(all(lw[1,] == 10 * lp[1,]), failure_message = "failed")
    expect(all(lw[-1,] == lp[-1,]), failure_message = "failed")
  }
  for (i in 1:5) {
    W <- rbm(20,20)
    W[3,3] <- 0.5
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    ind <- which(rownames(lp) == "r3 -- c3")
    expect(all(lw[ind,] == 0.5 * lp[ind,]), failure_message = "failed")
    expect(all(lw[-ind,] == lp[-ind,]), failure_message = "failed")
  }
})

test_that("Check matrices with several weighted links", {
  for (i in 1:5) {
    W <- rbm(20,20)
    W[7,8] <- 1.5
    W[4,6] <- 0.7
    lp <- link_positions(W, six_node = TRUE, weights = FALSE, normalisation = "none")
    lw <- link_positions(W, six_node = TRUE, weights = TRUE, normalisation = "none")
    ind <- which(rownames(lp) == "r7 -- c8")
    ind2 <- which(rownames(lp) == "r4 -- c6")
    expect(all(lw[ind,] == 1.5 * lp[ind,]), failure_message = "failed")
    expect(all(lw[ind2,] == 0.7 * lp[ind2,]), failure_message = "failed")
    expect(all(lw[-c(ind, ind2),] == lp[-c(ind, ind2),]), failure_message = "failed")
  }
})
