context("PCR")

test_that("PCR works", {
  expect_is(PCR_reconstruction(NPannual, NPpc, start.year = 1200),
            "list")
})

test_that("PCR cross-validation works", {
  expect_is(cvPCR(NPannual, NPpc, start.year = 1200), "list")
})

test_that("PCR ensemble works", {
  expect_is(PCR_reconstruction(NPannual, list(NPpc, NPpc[, 1:2]), start.year = 1200), "list")
})

test_that("PCR ensemble cross-validation works", {
  expect_is(cvPCR(NPannual, list(NPpc, NPpc[, 1:2]), start.year = 1200), "list")
})
