context("valgrind stressor")

test_that("pminC large", {
  x <- rep_len(c(0, -0.5, 1.2), 7e6)
  for (i in 1:128) {
    z <- pminC(x, 0)
    expect_equal(range_rcpp(z)[2], 0)
    expect_equal(range_rcpp(x)[2], 1.2)
  }
  expect_true(is.double(x))
})

test_that("bench_system_time without bench", {
  expr <- x <- rep_len(c(0, -0.5, 1.2), 5e6)
  bench_x <- (pminC(x, 0))
  expect_true(is.numeric(bench_x))
})

