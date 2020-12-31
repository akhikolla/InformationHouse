
if (startsWith(Sys.getenv("USERNAME"), "hugh")) {
  options(hutilscpp.nThread = 10)
}

test_that("start", {
  expect_true(TRUE)
})
