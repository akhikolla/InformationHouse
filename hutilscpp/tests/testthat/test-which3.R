test_that("which3 works", {
  skip_if_not_installed("data.table")
  expect_equal(which3(c(TRUE, FALSE, FALSE),
                      c(FALSE, TRUE, FALSE),
                      c(FALSE, FALSE, TRUE)),
               integer(0))
  expect_equal(which3(c(TRUE, FALSE, FALSE),
                      c(FALSE, TRUE, FALSE),
                      c(FALSE, FALSE, TRUE),
                      And = FALSE),
               1:3)
  expect_equal(which3(c(TRUE, NA, FALSE),
                      c(FALSE, TRUE, FALSE),
                      c(FALSE, FALSE, TRUE),
                      And = FALSE),
               1:3)

  x <- c(TRUE, FALSE, NA)
  expect_equal(which3(x, x, x),
               which(x & x & x))
  y <- c(TRUE, TRUE, TRUE)
  expect_equal(which3(x, y, x),
               which(x & y & x))

  library(data.table)
  TFN <- CJ(a = c(TRUE, FALSE, NA),
            b = c(TRUE, NA, FALSE),
            d = c(FALSE, TRUE, NA),
            e = c(FALSE, FALSE, FALSE),
            f = c(TRUE, TRUE, TRUE),
            g = c(NA, NA, NA))

  base1 <- TFN[, which(a & b & d)]
  huti1 <- TFN[, which3(a, b, d)]
  expect_equal(base1, huti1)

  base2 <- TFN[, which(b & a & d)]
  huti2 <- TFN[, which3(b, a, d)]
  expect_equal(base2, huti2)

  base3 <- TFN[, which(a & d & e)]
  huti3 <- TFN[, which3(a, d, e)]
  expect_equal(base3, huti3)

  base4 <- TFN[, which(d & e & f)]
  huti4 <- TFN[, which3(d, e, f)]
  expect_equal(base4, huti4)

  base5 <- TFN[, which(e & f & g)]
  huti5 <- TFN[, which3(e, f, g)]
  expect_equal(base5, huti5)

  base6 <- TFN[, which(a | b | d)]
  huti6 <- TFN[, which3(a, b, d, And = FALSE)]
  expect_equal(base6, huti6)





})


test_that("long failing instance", {
  skip_on_cran()
  N <- 1e8L
  if (nzchar(Sys.getenv("TRAVIS"))) {
    N <- 1e6L
  }

  x <- rep_len(1:20 %in% c(1, 5, 8, 9), N)
  y <- rep_len(1:10 %in% c(1, 5, 8, 9), N)
  z <- rep_len(TRUE, N)

  expect_identical(which3(x, y, z),
                   which(x & y & z))

  x <- rep_len(1:20 %in% c(2, 5, 8, 9), N)
  y <- rep_len(1:10 %in% c(1, 5, 8, 9), N)
  z <- rep_len(TRUE, N)

  expect_identical(which3(x, y, z),
                   which(x & y & z))

  z <- logical(N)

  expect_identical(which3(x, y, z),
                   which(x & y & z))

})


