library(testthat)
context("Testset generation")

test_that("default inflation to 0", {
  testSetLength = 20
  
  testSet_matrix = getTestSet(
    set = codeSet, 
    testSetLength = testSetLength
  )
  
  testSetBR = baserate(testSet_matrix)
  
  expect_equal(colnames(testSet_matrix), c("first_rater", "second_rater"))
  expect_equal(nrow(testSet_matrix), testSetLength)
  expect_gte(testSetBR$firstBaserate, 0)
})

test_that("custom inflation", {
  testSetLength = 20
  testSetInflation = 0.2
  
  testSet_matrix = getTestSet(
    set = codeSet, 
    testSetLength = testSetLength,
    testSetBaserateInflation = testSetInflation
  )
  
  testSetBR = baserate(testSet_matrix)
  
  expect_equal(colnames(testSet_matrix), c("first_rater", "second_rater"))
  expect_equal(nrow(testSet_matrix), testSetLength)
  expect_gte(testSetBR$firstBaserate, testSetInflation)
})

test_that("verify error checks", {
  expect_error(getTestSet(codeSet, 200), regexp = "testSetLength must be less than the length of set")
  expect_error(getTestSet(codeSet, -1), regexp = "testSetLength value must be a positive integer.")
  expect_error(getTestSet(codeSet, 20, -1), regexp = "must be positive")
  expect_error(getTestSet(codeSet, 20, 1.1), regexp = "must be below 1")
})

test_that("check short-circuits", {
  sample_set <- matrix(1, nrow=20, ncol=2)
  set_1 <- getTestSet(sample_set, 10)
  expect_equal(sum(set_1), 20)
  
  sample_set <- matrix(0, nrow=20, ncol=2)
  set_2 <- getTestSet(sample_set, 20)
  expect_equal(sum(set_2), 0)
  expect_error(getTestSet(sample_set, 10, 0.2), regexp = "Not enough positives")
})

test_that("Verify set kappas", {
  sample_set <- matrix(1, nrow=20, ncol=2)
  sample_k <- kappaSet(sample_set)
  expect_equal(sample_k, 1)
  
  sample_set_bad <- sample_set
  sample_set_bad[[1]] = -1
  expect_error(kappaSet(sample_set_bad), regexp = "Each set must only consist of zeros and ones.")
})