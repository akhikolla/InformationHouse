context("mcount")

# test against subgraph extraction method
test_that("Comparing with motif extraction algorithm", {
  for (i in 1:10) {
    W <- rbm(6,6)
    mc <- mcount(W, six_node = TRUE, mean_weight = FALSE, standard_dev = FALSE, normalisation = FALSE)$frequency
    mc_minors <- mcount_minors(W = W, six_node = TRUE)
    expect_equal(mc, mc_minors)
  }
})

# tests with block matrices
mlist2 <- motifs[1]
bm2 <- block_matrix(mlist2)

test_that("Test counts for the motifs with 2 nodes", {
  mc <- mcount(M = bm2, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  expect_identical(all(mc[mc$nodes == 2,"frequency"] == 1), TRUE)
})

mlist3 <- motifs[2:3]
bm3 <- block_matrix(mlist3)

test_that("Test counts for the motifs with 3 nodes", {
  mc <- mcount(M = bm3, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  expect_identical(all(mc[mc$nodes == 3,"frequency"] == 1), TRUE)
})

mlist4 <- motifs[4:7]
bm4 <- block_matrix(mlist4)

test_that("Test counts for the motifs with 4 nodes", {
  mc <- mcount(M = bm4, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  expect_identical(all(mc[mc$nodes == 4,"frequency"] == 1), TRUE)
})

mlist5 <- motifs[8:17]
bm5 <- block_matrix(mlist5)

test_that("Test counts for the motifs with 5 nodes", {
  mc <- mcount(M = bm5, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  expect_identical(all(mc[mc$nodes == 5,"frequency"] == 1), TRUE)
})


mlist6 <- motifs[18:44]
bm6 <- block_matrix(mlist6)

test_that("Test counts for the motifs with 6 nodes", {
  mc <- mcount(M = bm6, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  expect_identical(all(mc[mc$nodes == 6,"frequency"] == 1), TRUE)
})

# error messages
test_that("Error messages work",{
  expect_error(object = mcount("a"), "'M' must be an object of class 'matrix'")
  expect_error(object = mcount(matrix("a",3,3)), "Elements of 'M' must be numeric")
  expect_error(object = mcount(matrix(-1,3,3)), "Elements of 'M' must be greater than or equal to zero")
  expect_error(object = mcount(matrix(1,3,3), six_node = TRUE, normalisation = 7), "'normalisation' must be of class 'logical' i.e. TRUE or FALSE")
  expect_error(object = mcount(matrix(1,3,3), six_node = 7, normalisation = TRUE), "'six_node' must be of class 'logical' i.e. TRUE or FALSE")
})

# validated results
test_that("Matches validated results",{
  expect_equal(object = mcount(M = mat, six_node = FALSE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_SF_NF)
  expect_equal(object = mcount(M = mat, six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_ST_NF)
  expect_equal(object = mcount(M = mat, six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)[,-6], mcount_SF_NT)
  expect_equal(object = mcount(M = mat, six_node = TRUE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)[,-6], mcount_ST_NT)
  expect_equal(mcount(M = t(mat), six_node = TRUE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)[,-6], mcount_transpose_ST_NT)
  expect_equal(mcount(M = t(mat), six_node = FALSE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)[,-6], mcount_transpose_SF_NT)
  expect_equal(mcount(M = t(mat), six_node = TRUE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_ST_NF)
  expect_equal(mcount(M = t(mat), six_node = FALSE, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE), mcount_transpose_SF_NF)
})
