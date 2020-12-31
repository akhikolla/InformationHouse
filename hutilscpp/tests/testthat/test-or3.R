context("test-or3")

test_that("or3 error handling", {
  expect_error(or3(logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(z = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(y = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(y = TRUE, logical(2), x = logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(x = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(x = logical(3), logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(or3(x = logical(1), logical(4), logical(3)), regexp = "permissible.*lengths")
})

test_that("logical3 works", {
  x <- c(TRUE, FALSE, TRUE)
  y <- logical(3)
  z <- c(TRUE, TRUE, FALSE)
  expect_equal(or3(x, y, z), x | y | z)
  z <- NULL
  expect_equal(or3(x, y, z), x | y)
})

test_that("or3 length-1", {
  expect_equal(or3(TRUE, FALSE, TRUE),
               TRUE)
})

test_that("length-1s 2/3", {
  expect_equal(or3(FALSE, FALSE, c(TRUE, FALSE, FALSE)),
               c(TRUE, FALSE, FALSE))
  expect_equal(or3(logical(3), FALSE, c(TRUE, FALSE, FALSE)),
               c(TRUE, FALSE, FALSE))
  expect_equal(or3(FALSE, logical(3), c(TRUE, FALSE, FALSE)),
               c(TRUE, FALSE, FALSE))
})

test_that("or3 works with NAs", {
  skip_if_not_installed("data.table")
  library(data.table)
  DT <- CJ(xx = c(NA, FALSE, TRUE),
           yy = c(NA, FALSE, TRUE),
           zz = c(NA, FALSE, TRUE))
  DT[, res1 := or3(xx, yy, zz)]
  DT[, res2 := xx | yy | zz]
  expect_equal(DT[["res1"]], DT[["res2"]])

  DT[, res3 := or3(xx, yy, TRUE)]
  DT[, res4 := xx | yy | TRUE]
  expect_equal(DT[["res3"]], DT[["res4"]])
  DT[, res3 := or3(xx, yy, FALSE)]
  DT[, res4 := xx | yy | FALSE]
  expect_equal(DT[["res3"]], DT[["res4"]])

  for (pos in 1:3) {
    for (lgvl in c(TRUE, FALSE, NA)) {
      switch(pos,
             {
               DT[, res5 := or3(lgvl, xx, yy)]
             },
             {
               DT[, res5 := or3(xx, lgvl, yy)]
             },
             {
               DT[, res5 := or3(xx, yy, lgvl)]
             })
      DT[, res6 := xx | yy | lgvl]
      expect_equal(DT[["res5"]], DT[["res6"]],
                   info = paste("single", pos, lgvl))
    }
  }
  for (x in c(TRUE, FALSE, NA)) {
    for (y in c(TRUE, FALSE, NA)) {
      for (z in c(TRUE, FALSE, NA)) {
        expect_equal(or3(c(x, z), y, x),
                     c(x, z) | y | x,
                     info = paste("1", x, y, z))
        expect_equal(or3(y = c(x, y), y, x),
                     c(x, y) | y | x,
                     info = paste("2", x, y, z))
      }
    }
  }
})

test_that("C++", {
  expect_error(do_or3(logical(3), logical(2), TRUE),
               regexp = "different lengths")
  expect_equal(do_or3(c(TRUE, FALSE), c(TRUE, FALSE), TRUE), c(TRUE, TRUE))
  expect_equal(do_or3(c(TRUE, FALSE), c(TRUE, FALSE), FALSE), c(TRUE, FALSE))
  expect_equal(do_or3(c(TRUE, FALSE), c(TRUE, FALSE), logical(0)), c(TRUE, FALSE))
  expect_error(do_or3(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE, FALSE)),
               regexp = "wrong length")
})
