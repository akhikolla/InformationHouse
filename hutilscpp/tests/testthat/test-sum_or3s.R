test_that("sum_or3s works", {
  sum_bor3 <- function(exprA, exprB = TRUE, exprC = TRUE, ..., .parent_nframes = 1L, nThread = 1L) {
    if (missing(..1)) {
      sum(exprA | exprB | exprC)
    } else if (missing(..2)) {
      return(exprA | exprB)
    } else {
      sum_bor3(exprA | exprB, exprC, ...)
    }
  }
  abc <- 1:100
  def <- 1:100
  ghi <- 1:100
  x <- abc == 1
  expect_equal(sum_or3s(x, def >= 1, ghi <= 2), 100)
  expect_equal(sum_or3s(x, x, x), 1)

  expect_equal(sum_or3s(abc == 1, def >= 1, ghi <= 2), 100)

  S3 <- 100:1001
  B3 <- c(99L, 100:1001)

  DT <- data.table(A = 1:100, B = 2L, Z = 3L)
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), B == 2L, Z >= 1L)],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), B == 2L, Z >= 1L)])
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), B != 2L, Z >= 1L)],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), B != 2L, Z >= 1L)])
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), B != 2L, Z %between% c(1L, 4L))],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), B != 2L, Z %between% c(1L, 4L))])
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), Z %between% c(1L, 4L))],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), Z %between% c(1L, 4L))])
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), B != 2L, Z %in% 1:4)],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), B != 2L, Z %in% 1:4)])
  expect_equal(DT[, sum_or3s(A %in% c(5L, 5L, 6L), B != 2L, Z %in% c(2L, 1:4))],
               DT[, sum_bor3(A %in% c(5L, 5L, 6L), B != 2L, Z %in% c(2L, 1:4))])
  expect_equal(DT[, sum_or3s(A %in% c(5L, integer(101)), B != 2L, Z %in% c(2L, 1:4))],
               DT[, sum_bor3(A %in% c(5L, integer(101)), B != 2L, Z %in% c(2L, 1:4))])
  expect_equal(DT[, sum_or3s(A %in% c(5L, integer(101)), B %in% S3, Z %in% c(integer(101)))],
               DT[, sum_bor3(A %in% c(5L, integer(101)), B %in% S3, Z %in% c(integer(101)))])
  expect_equal(DT[, sum_or3s(A %in% c(5L, integer(101)), B %in% B3, Z %in% c(integer(101)))],
               DT[, sum_bor3(A %in% c(5L, integer(101)), B %in% B3, Z %in% c(integer(101)))])
})


test_that("internal-betweens", {
  x <- 1:10 + 0L
  expect_equal(sum_or3s(x %between% c(1L, 3L)), 3L)
  expect_equal(sum_or3s(x %(between)% c(1L, 3L)), 1L)
  expect_equal(sum_or3s(x %]between[% c(1L, 10L)), 2L)
  x <- as.double(x)
  expect_equal(sum_or3s(x %between% c(1, 3)), 3L)
  expect_equal(sum_or3s(x %(between)% c(1, 3)), 1L)
  expect_equal(sum_or3s(x %]between[% c(1, 10)), 2L)

})

test_that("mixture double ints sum_or3s", {
  sum_bor3 <- function(exprA, exprB = TRUE, exprC = TRUE, ..., .parent_nframes = 1L, nThread = 1L) {
    if (missing(..1)) {
      sum(exprA | exprB | exprC)
    } else if (missing(..2)) {
      return(exprA | exprB)
    } else {
      sum_bor3(exprA | exprB, exprC, ...)
    }
  }

  y <- c(101L, 103L, 110L, 999L)
  expect_equal(sum_or3s(y == 103.5,
                        y <= 103.5,
                        y %in% c(998:999)),
               3)
  expect_equal(sum_or3s(y <= 103.5,
                        y %in% c(998:999, 1000:9999)),
               3)
  y <- as.double(y)
  expect_equal(sum_or3s(y == 103.5,
                        y <= 103.5,
                        y %in% c(998:999)),
               3)


})


