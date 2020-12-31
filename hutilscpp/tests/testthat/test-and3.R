context("test-and3")

test_that("and3 error handling", {
  expect_error(and3(logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(z = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(y = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(y = TRUE, logical(2), x = logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(x = TRUE, logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(x = logical(3), logical(2), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(x = logical(4), logical(1), logical(3)), regexp = "permissible.*lengths")
  expect_error(and3(y = logical(4), logical(1), logical(3)), regexp = "permissible.*lengths")
})

test_that("and3 works", {
  x <- c(TRUE, FALSE, TRUE)
  y <- logical(3)
  z <- c(TRUE, TRUE, FALSE)
  expect_equal(and3(x, y, z),
               x & y & z)
  expect_equal(and3(x, y), x & y)
})

test_that("and3 length-1s", {
  expect_equal(and3(TRUE, TRUE, TRUE), TRUE)
  expect_equal(and3(TRUE, TRUE, FALSE), FALSE)
  expect_equal(and3(TRUE, TRUE, NA), NA)
})

test_that("length-1s 2/3", {
  expect_equal(and3(TRUE, TRUE, c(TRUE, FALSE, TRUE)),
               c(TRUE, FALSE, TRUE))
  expect_equal(and3(TRUE, logical(3), c(TRUE, FALSE, TRUE)),
               c(FALSE, FALSE, FALSE))
})

test_that("and3 works with NAs", {
  skip_if_not_installed("data.table")
  library(data.table)
  DT <- CJ(xx = c(NA, FALSE, TRUE),
           yy = c(NA, FALSE, TRUE),
           zz = c(NA, FALSE, TRUE))
  DT[, res1 := and3(xx, yy, zz)]
  DT[, res2 := xx & yy & zz]
  expect_equal(DT[["res1"]], DT[["res2"]])

  DT[, res3 := and3(xx, yy, TRUE)]
  DT[, res4 := xx & yy & TRUE]
  expect_equal(DT[["res3"]], DT[["res4"]])

  for (pos in 1:3) {
    for (lgvl in c(TRUE, FALSE, NA)) {
      switch(pos,
             {
               DT[, res5 := and3(lgvl, xx, yy)]
             },
             {
               DT[, res5 := and3(xx, lgvl, yy)]
             },
             {
               DT[, res5 := and3(xx, yy, lgvl)]
             })
      DT[, res6 := xx & yy & lgvl]
      expect_equal(DT[["res5"]], DT[["res6"]],
                   info = paste("single", pos, lgvl))
    }
  }
  for (x in c(TRUE, FALSE, NA)) {
    for (y in c(TRUE, FALSE, NA)) {
      for (z in c(TRUE, FALSE, NA)) {
        expect_equal(and3(c(x, z), y, x),
                     c(x, z) & y & x,
                     info = paste("1", x, y, z))
        expect_equal(and3(y = c(x, y), y, x),
                     c(x, y) & y & x,
                     info = paste("2", x, y, z))
      }
    }
  }
})

test_that("nas_absent ok", {
  expect_equal(and3(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE), nas_absent = TRUE),
               c(TRUE, FALSE))
})

test_that("C++", {
  expect_error(do_and3(logical(3), logical(2), TRUE),
               regexp = "lengths")
  expect_error(do_and3(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE, FALSE)),
               regexp = "wrong length")
  expect_equal(do_and3(c(TRUE, FALSE), c(TRUE, FALSE), TRUE), c(TRUE, FALSE))
  expect_equal(do_and3(c(TRUE, FALSE), c(TRUE, FALSE), FALSE), c(FALSE, FALSE))
  expect_equal(do_and3(c(TRUE, FALSE), c(TRUE, FALSE), logical(0)), c(TRUE, FALSE))
  expect_equal(na_and(c(TRUE, FALSE, NA)), c(TRUE, FALSE, NA) & NA)
  expect_error(do_and3(c(TRUE, FALSE), c(TRUE, FALSE), c(TRUE, FALSE, FALSE)),
               regexp = "wrong length")
})

test_that("do_and3_na", {
  skip_if_not_installed("data.table")
  skip_if_not_installed("hutils")
  library(data.table)
  library(hutils)
  DT <- CJ(x = c(TRUE, FALSE, NA), y = c(TRUE, FALSE, NA), z = c(TRUE, FALSE, NA))
  DT[, xf := coalesce(x, FALSE)]
  DT[, xt := coalesce(x, TRUE)]
  DT[, yf := coalesce(y, FALSE)]
  DT[, yt := coalesce(y, TRUE)]
  DT[, zf := coalesce(z, FALSE)]
  DT[, zt := coalesce(z, TRUE)]
  DT[, ans0 := do_and3_na(x, y, z)]
  DT[, ans1 := do_and3_na(x, y, z, na_value = -1L)]
  DT[, ans2 := do_and3_na(x, y, z, na_value = +1L)]

  DT[, ans0x := do_and3_na(.BY[[1]], y, z), by = "x"]
  DT[, ans1x := do_and3_na(.BY[[1]], y, z, na_value = -1L), by = "x"]
  DT[, ans2x := do_and3_na(.BY[[1]], y, z, na_value = +1L), by = "x"]
  DT[, expect_identical(ans0, ans0x)]
  DT[, expect_identical(ans1, ans1x)]
  DT[, expect_identical(ans2, ans2x)]

  DT[, ans0y := do_and3_na(x, .BY[[1]], z), by = "y"]
  DT[, ans1y := do_and3_na(x, .BY[[1]], z, na_value = -1L), by = "y"]
  DT[, ans2y := do_and3_na(x, .BY[[1]], z, na_value = +1L), by = "y"]
  DT[, expect_identical(ans0, ans0y)]
  DT[, expect_identical(ans1, ans1y)]
  DT[, expect_identical(ans2, ans2y)]

  DT[, ans0z := do_and3_na(x, y, .BY[[1]]), by = "z"]
  DT[, ans1z := do_and3_na(x, y, .BY[[1]], na_value = -1L), by = "z"]
  DT[, ans2z := do_and3_na(x, y, .BY[[1]], na_value = +1L), by = "z"]
  DT[, expect_identical(ans0, ans0z)]
  DT[, expect_identical(ans1, ans1z)]
  DT[, expect_identical(ans2, ans2z)]

  DT[, ans0xyz := do_and3_na(.BY[[1]], .BY[[2]], .BY[["z"]]), by = c("x", "y", "z")]
  DT[, expect_identical(ans0, ans0xyz)]

  # default is any NA => FALSE
  DT[seq_len(.N - 1L), expect_false(any(ans1))]

  expect_equal(do_and3_na(logical(0), logical(0), logical(0)), logical(0))
})

test_that("prev segfault", {
  a <- do_and3_na(TRUE, c(TRUE, TRUE), TRUE, maxCall = 5)
  expect_true(all(a))
})

test_that("and3_na errors", {
  expect_error(do_and3_na(TRUE, TRUE, FALSE, maxCall = -1), "maxCall")
  expect_error(do_and3_na(TRUE, c(TRUE, FALSE), c(TRUE, FALSE, NA), maxCall = 1), "lengths")
})









