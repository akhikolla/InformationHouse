context("which_first")

test_that("which_first works", {
  x <- runif(10, 1, 2)
  expr <- c(0, 5, 3, 2)
  expect_identical(which_first(x < 0), 0L)
  expect_identical(which_first(x > 0), 1L)
  expect_identical(which_first(x == 0), 0L)
  expect_identical(which_first(x >= 0), 1L)
  expect_identical(which_first(x <= 0), 0L)
  expect_identical(which_first(x != 0), 1L)
  x30 <- rep(3, 4)
  expect_identical(which_first(x30 != 3.0), 0L)

  x <- c(-2, -1.5)
  expect_identical(which_first(x == -1.5), 2L)
  expect_identical(which_first(x == -2), 1L)
  expect_identical(which_first(x == -2L), 1L)
  expect_identical(which_first(x == 2L), 0L)
  expect_identical(which_first(x == 2), 0L)

  expect_identical(which_first(x >= -1.5), 2L)
  expect_identical(which_first(x >= -2), 1L)
  expect_identical(which_first(x >= -2L), 1L)
  expect_identical(which_first(x >= 2L), 0L)
  expect_identical(which_first(x >= 2), 0L)

  expect_identical(which_first(x > -1.5), 0L)
  expect_identical(which_first(x > -1.6), 2L)
  expect_identical(which_first(x > -2), 2L)
  expect_identical(which_first(x > -2L), 2L)
  expect_identical(which_first(x > 2L), 0L)
  expect_identical(which_first(x > 2), 0L)

  expect_identical(which_first(x <= -1.5), 1L)
  expect_identical(which_first(x <= -2), 1L)
  expect_identical(which_first(x <= -2L), 1L)
  expect_identical(which_first(x <= 2L), 1L)
  expect_identical(which_first(x <= -32L), 0L)
  expect_identical(which_first(x <= -32), 0L)
  expect_identical(which_first(x <= 2), 1L)
  expect_identical(which_first(x <= -2.2), 0L)

  expect_identical(which_first(x < -1.5), 1L)
  expect_identical(which_first(x < -2), 0L)
  expect_identical(which_first(x < -2L), 0L)
  expect_identical(which_first(x < 2L), 1L)
  expect_identical(which_first(x < -1L), 1L)
  expect_identical(which_first(x < -32), 0L)
  expect_identical(which_first(x < 2), 1L)

  expect_identical(which_first(x != -1.5), 1L)
  expect_identical(which_first(x != -2), 2L)
  expect_identical(which_first(x != -2L), 2L)
  expect_identical(which_first(x != 2L), 1L)
  x <- -1
  expect_identical(which_first(x != -1), 0L)
  expect_identical(which_first(x != -1L), 0L)


  y <- -5:5
  y0 <- integer(0)
  yr <- rev(y)
  expect_identical(which_first(y == -4), 2L)
  expect_identical(which_first(y == -4L), 2L)
  expect_identical(which_first(y == -6L), 0L)
  expect_identical(which_first(y == -5.5), 0L)
  expect_identical(which_first(y == -0.01), 0L)
  expect_identical(which_first(y == 5), length(y))
  expect_identical(which_first(y == 5L), length(y))
  expect_identical(which_first(y == 2.5), 0L)

  expect_identical(which_first(y != -4), 1L)
  expect_identical(which_first(y != -5), 2L)
  expect_identical(which_first(y != -5L), 2L)
  expect_identical(which_first(y != 2.5), 1L)
  # corner case
  expect_identical(which_first(y0 != 0), 0L)
  expect_identical(which_first(y0 != 0L), 0L)
  expect_identical(which_first(y0 != 0.5), 0L)


  expect_identical(which_first(y >= 5), length(y))
  expect_identical(which_first(y >= 5L), length(y))
  expect_identical(which_first(y >= 55), 0L)
  expect_identical(which_first(y >= 55L), 0L)
  expect_identical(which_first(y >= -5L), 1L)
  expect_identical(which_first(y >= 4.5), length(y))
  expect_identical(which_first(y >= -4.5), 2L)

  expect_identical(which_first(y <= 5), 1L)
  expect_identical(which_first(y <= 5.5), 1L)
  expect_identical(which_first(y <= 1), 1L)
  expect_identical(which_first(y <= -5.5), 0L)
  expect_identical(which_first(yr <= -5), length(yr))

  expect_identical(which_first(y < 5), 1L)
  expect_identical(which_first(y < 5L), 1L)
  expect_identical(which_first(y < 5.5), 1L)
  expect_identical(which_first(y < -4.5), 1L)
  expect_identical(which_first(y < -5.5), 0L)
  expect_identical(which_first(y < 1), 1L)
  expect_identical(which_first(y < 1L), 1L)
  expect_identical(which_first(y < -11L), 0L)
  expect_identical(which_first(yr < 1.5), 5L)
  expect_identical(which_first(yr < -4), length(yr))

  expect_identical(which_first(y > 0), 7L)
  expect_identical(which_first(y > -1), 6L)
  expect_identical(which_first(y > -2.5),
                   any(y > -2.5) * which.max(y > -2.5))

  y2 <- as.integer(c(0, -1, -2, -1, 0))
  y3 <- as.integer(c(2, 1, 0, 3))
  expect_identical(which_first(y2 <= -2.5),
                   Position(f = function(k) k <= -2.5,
                            x = y2,
                            nomatch = 0L))
  expect_identical(which_first(y3 <= 0.5),
                   Position(f = function(k) k <= 0.5,
                            x = y3,
                            nomatch = 0L))

  expr <- c(0, 5, 3, 2)
  expect_identical(which_first(expr %in% c(22, 32)), 0L)
  expect_identical(which_first(expr %in% c(2, 3)), 3L)
  expr <- as.integer(expr)
  expect_identical(which_first(expr %in% c(2, 3)), 3L)
  expect_identical(which_first(expr %in% c(7, 9)), 0L)
  expect_identical(which_first(expr %in% c(2L, 3L)), 3L)
  expect_identical(which_first(expr %in% c(8L, 9L)), 0L)
})

test_that("Conflicts with expressions", {
  expr <- c(0, 5, 3, 2)
  expect_identical(which_first(expr == 5), 2L)
  lhs <- rhs <- 3
  expect_identical(which_first(expr == lhs), 3L)
})

test_that("Fall through", {
  y <- logical(5)
  expect_identical(which_first(y), 0L)
  expect_identical(which_first(!y), 1L)
  expect_message(which_first(!y,
                             verbose = TRUE),
                 regexp = "which.max",
                 fixed = TRUE)

  x <- c(letters, 1)
  expect_identical(which_first(x == 1), 27L)
  expect_identical(which_first(x == 155), 0L)
  expect_identical(which_first(x < 1), 0L)
})

test_that("match", {
  int_m <- 1:100
  expect_identical(which_first(int_m == 2L), 2L)
})



test_that(".which_first_logical all NA", {
  skip("For later consideration.")
  x <- c(NA, NA, NA)
  expect_equal(.which_first_logical(x, NA), 1)
  expect_equal(.which_first_logical(x, TRUE), 0)
  expect_equal(.which_first_logical(x, FALSE), 0)
})
test_that(".which_first_logical one NA", {
  skip("For later consideration.")
  x <- c(NA, NA, TRUE)
  expect_equal(.which_first_logical(x, NA), 1)
  expect_equal(.which_first_logical(x, TRUE), 3)
  expect_equal(.which_first_logical(x, FALSE), 0)
  y <- c(NA, NA, FALSE)
  expect_equal(.which_first_logical(y, NA), 1)
  expect_equal(.which_first_logical(y, TRUE), 0)
  expect_equal(.which_first_logical(y, FALSE), 3)
  z <- c(TRUE, NA, NA)
  expect_equal(.which_first_logical(z, NA), 2)
  expect_equal(.which_first_logical(z, TRUE), 1)
  expect_equal(.which_first_logical(z, FALSE), 0)
})
test_that(".which_first_logical every", {
  skip("For later consideration.")
  x <- c(NA, TRUE, FALSE)
  expect_equal(.which_first_logical(x, NA), 1)
  expect_equal(.which_first_logical(x, TRUE), 2)
  expect_equal(.which_first_logical(x, FALSE), 3)
  y <- c(FALSE, NA, TRUE)
  expect_equal(.which_first_logical(y, NA), 2)
  expect_equal(.which_first_logical(y, TRUE), 3)
  expect_equal(.which_first_logical(y, FALSE), 1)
})
test_that(".which_first_logical no NA", {
  x <- c(TRUE, TRUE, TRUE)

  expect_equal(.which_first_logical(x, TRUE), 1)
  expect_equal(.which_first_logical(x, FALSE), 0)

  y <- c(FALSE, FALSE, FALSE)

  expect_equal(.which_first_logical(y, TRUE), 0)
  expect_equal(.which_first_logical(y, FALSE), 1)

  z <- c(FALSE, TRUE, FALSE)

  expect_equal(.which_first_logical(z, TRUE), 2)
  expect_equal(.which_first_logical(z, FALSE), 1)
})

test_that("unexpected o", {
  x <- c(NA, TRUE, TRUE)
  expect_equal(which_first(x == TRUE), 2)
  expect_equal(which_first(x == 1L), 2)
  expect_equal(which_first(x == 1), 2)
  expect_equal(which_first(x != 1), 0)

  x.raw <- c(raw(5), charToRaw("A"))
  expect_equal(which_first(x.raw == 0), 1L)
  expect_equal(which_first(x.raw != 0), 6L)
  expect_equal(which_first(x.raw < 0), 0L)
  expect_equal(which_first(x.raw <= 0), 1L)
  expect_equal(which_first(x.raw > 0), 6L)
  expect_equal(which_first(x.raw >= 0), 1L)

})

test_that("LHS logical length-one", {
  true <- TRUE
  expect_equal(which_first(true == 1L), 1L)
  expect_equal(which_first(true == 0L), 0L)
  expect_equal(suppressWarnings(which_first(true == NA)), 0L)
  false <- FALSE
  expect_equal(which_first(false == 1), 0L)
  expect_equal(which_first(false == 0), 1L)
  expect_equal(suppressWarnings(which_first(false == NA)), 0L)
  missy <- NA
  expect_equal(which_first(missy == 1), 0L)
  expect_equal(which_first(missy == 1), 0L)
})

test_that("RHS NA", {
  x <- c(NA, NA)
  expect_error(which_first(x > NA),
               regexp = "This is not supported for operator '>'.",
               fixed = TRUE)
  expect_warning(wf_xisna <- which_first(x == NA),
                 regexp = "`rhs` appears to be logical NA.",
                 fixed = TRUE)
  expect_equal(wf_xisna, 1L)
  wf_xisfalse <- which_first(x == 0)
  expect_equal(wf_xisfalse, 0L)

  y <- c(TRUE, FALSE, NA)
  expect_warning(wf_yisna <- which_first(y == NA),
                 regexp = "which_first(is.na",
                 fixed = TRUE)
  expect_equal(wf_yisna, 3L)
  expect_warning(wf_yisntna <- which_first(y != NA),
                 regexp = "which_first(!is.na",
                 fixed = TRUE)
  expect_equal(wf_yisntna, 1L)

  z <- c(NA, FALSE)
  expect_warning(wf_zisntna <- which_first(z != NA),
                 regexp = "`rhs` appears to be logical NA.",
                 fixed = TRUE)
  expect_equal(wf_zisntna, 2L)

})

test_that("lhs_eval length 0", {
  x <- integer(0)
  expect_equal(which_first(x == 0.5), 0L)
  expect_equal(which_first(x != 0.5), 0L)
  expect_equal(which_first(x >= 0.5), 0L)
  expect_equal(which_first(x <= 0.5), 0L)
  expect_equal(which_first(x < 0.5), 0L)
  expect_equal(which_first(x > 0.5), 0L)
  expect_equal(.which_first(c(NA, NA)), 0L)
})


test_that("which_first logical", {
  #10
  x <- y <- c(TRUE, FALSE, NA)
  expect_equal(which_first(x != y), 0)
  expect_equal(which_first(x == y), 1)
  expect_equal(which_first(x > y), 0)
  expect_equal(which_first(x < y), 0)

  y <- as.integer(x)
  expect_equal(which_first(x != y), 0)

  y <- c(2L, 0L, NA)
  expect_equal(which_first(x != y), 1)


})

test_that("which_first_logical %in%", {

  abc <- logical(0)
  lgl0 <- logical(0)
  ttt <- c(TRUE, TRUE, TRUE)
  ttf <- c(TRUE, TRUE, FALSE)
  ttn <- c(TRUE, TRUE, NA)
  tfn <- c(TRUE, FALSE, NA)
  fff <- c(FALSE, FALSE, FALSE)
  nana <- c(NA, NA, NA)



  expect_equal(which_first(abc %in% lgl0),
               first_which(abc %in% lgl0))
  expect_equal(which_first(ttt %in% fff),
               first_which(ttt %in% fff))
  expect_equal(which_first(ttt %in% ttt),
               first_which(ttt %in% ttt))
  expect_equal(which_first(fff %in% fff),
               first_which(fff %in% fff))

  x <- c(TRUE, FALSE, TRUE, FALSE)
  y <- c(TRUE, FALSE)
  expect_equal(which_first(x %in% y), 1)

  x <- rev(x)
  expect_equal(which_first(x %in% y), 1)

  xx <- logical(10)
  yy <- logical(4)
  expect_equal(which_first(xx %in% yy), 1)
  yy <- !yy
  expect_equal(which_first(xx %in% yy), 0)

  yy <- c(TRUE, NA)
  expect_equal(which_first(xx %in% yy), 0)
  xx[3] <- TRUE
  expect_equal(which_first(xx %in% yy), 3)
  xx[3] <- NA
  expect_equal(which_first(xx %in% yy), 3)



  yyy <- TRUE
  expect_equal(which_first(xx %in% yyy), 0)
  yyy <- FALSE
  expect_equal(which_first(xx %in% yyy), 1)
  yyy <- NA
  expect_equal(which_first(xx %in% yyy), 3)

  xxx <- c(TRUE, FALSE)
  yyy <- NA
  expect_equal(which_first(xxx %in% yyy), 0)
  xxx <- c(TRUE, FALSE)
  yyy <- TRUE
  expect_equal(which_first(xxx %in% yyy), 1)
  yyy <- FALSE
  expect_equal(which_first(xxx %in% yyy), 2)

  yyy <- c(TRUE, FALSE, NA)
  expect_equal(which_first(xx %in% yyy), 1)


  a <- c(rep(TRUE, 100), NA)
  b <- c(TRUE, FALSE, TRUE)
  expect_equal(which_first(a %in% b), 1)
  a <- c(rep(FALSE, 100), NA)
  expect_equal(which_first(a %in% b), 1)
  b <- c(TRUE, TRUE)
  expect_equal(which_first(a %in% b), 0)
  b <- c(TRUE, FALSE)
  expect_equal(which_first(a %in% b), 1)
  b <- c(TRUE, NA)
  expect_equal(which_first(a %in% b), 101)

  a <- c(TRUE, TRUE, rep(NA, 100), FALSE, TRUE)
  b <- c(FALSE, NA)
  expect_equal(which_first(a %in% b), 3)
  b <- c(TRUE, NA)
  expect_equal(which_first(a %in% b), 1)

})




test_that("Error handling", {
  x <- c(TRUE, FALSE)
  y <- c(TRUE, FALSE, TRUE)
  expect_error(which_first(x == y), "length")
  x <- as.integer(x)
  y <- as.integer(y)
  expect_error(which_first(x >= y), "length")
  expect_error(which_first(x %between% y), "length")

})


test_that("which_first fall through when no name", {
  expect_message(which_first(1 == 1, verbose = TRUE),
                 regexp = "which.max")
})

test_that("which_first(<lgl> <operator> <TRUE/FALSE>)", {
  lhs <- (c(1:10) %% 2L) | (c(1:10) == 7L)

  expect_equal(which_first(lhs < TRUE), 2L)
  expect_equal(which_first(lhs <= TRUE), 1L)
  expect_equal(which_first(lhs > TRUE), 0L)
  expect_equal(which_first(lhs >= TRUE), 1L)

  expect_equal(which_first(lhs < FALSE), 0L)
  expect_equal(which_first(lhs <= FALSE), 2L)
  expect_equal(which_first(lhs > FALSE), 1L)
  expect_equal(which_first(lhs >= FALSE), 1L)

  expect_equal(which_last(lhs < TRUE), 10L)
  expect_equal(which_last(lhs <= TRUE), 10L)
  expect_equal(which_last(lhs > TRUE), 0L)
  expect_equal(which_last(lhs >= TRUE), 9L)

  expect_equal(which_last(lhs < FALSE), 0L)
  expect_equal(which_last(lhs <= FALSE), 10L)
  expect_equal(which_last(lhs > FALSE), 9L)
  expect_equal(which_last(lhs >= FALSE), 10L)


})

test_that("which_first() == notTRUE", {
  allTRUE3 <- rep(TRUE, 3)
  expect_equal(which_first(allTRUE3 != TRUE), 0L)
  allFALSE3 <- logical(3)
  expect_equal(which_first(allFALSE3 != FALSE), 0L)
})

test_that("which_first() N == N, N != N", {
  skip_if_not_installed("data.table")
  library(data.table)
  o <- c(-4:4, 4:0)
  x <- c(-4:4, 5:1)
  y <- o + 1L
  expect_equal(which_first(o == y), 0L)
  expect_equal(which_first(o == x), 1L)
  expect_equal(which_first(o != x), first(which(o != x)))
  expect_equal(which_first(o == x, reverse = TRUE), last(which(o == x)))
  expect_equal(which_first(o != x, reverse = TRUE), last(which(o != x)))

  od <- as.double(o)
  xd <- as.double(x)
  yd <- as.double(y)
  expect_equal(which_first(od == yd), 0L)
  expect_equal(which_first(od == xd), 1L)
  expect_equal(which_first(od != xd), first(which(od != xd)))
  expect_equal(which_first(od == xd, reverse = TRUE), last(which(od == xd)))
  expect_equal(which_first(od != xd, reverse = TRUE), last(which(od != xd)))
})

test_that("which_first(NA . NA)", {
  x <- c(NA, "a", "abc", "Def")
  y <- "abc"
  z <- NA_character_
  n <- ""
  expect_equal(which_first(x == y), 3L)
  expect_equal(which_first(x == z), 1L)
  expect_equal(which_first(x == ""), 0L)
  expect_equal(which_first(x == n), 0L)
  expect_equal(which_first(x != y), 1L)
  expect_equal(which_first(x != z), 2L)
  expect_equal(which_first(x != ""), 1L)
})

test_that("which_first not trues", {
  expect_equal(which_first(c(NA, NA)), 0L)
  expect_equal(which_last(c(NA, NA)), 0L)
})

test_that("which_first_quick", {
  x <- 11:20
  expect_equal(which_first(x == 13L, use.which.max = TRUE), 3L)
})

test_that("which_first_lgl_NA", {
  x <- c(TRUE, NA, TRUE, FALSE)
  expect_equal(which_first(x %in% c(FALSE, NA)), 2L)
  expect_equal(which_first(x %in% c(TRUE, NA)), 1L)
  x2 <- c(NA, x)
  expect_equal(which_first(x2 %in% c(TRUE, NA)), 1L)
})

test_that("which_first(x > 1) (len = 1)", {
  x <- 1
  expect_equal(which_first(x > 1), 0)
  expect_equal(which_first(x <= 1), 1)

  expect_equal(which_first(x != 1), 0)
  expect_equal(which_first(x == 1), 1)
  expect_equal(which_first(x >= 1), 1)
  expect_equal(which_first(x <= 1), 1)
  expect_equal(which_first(x > 1), 0)
  expect_equal(which_first(x < 1), 0)
  expect_equal(which_first(x %in% 1), 1)
  expect_equal(which_first(x %between% c(1L, 1L)), 1)
  expect_equal(which_first(x %(between)% c(1, 1)), 0)
  expect_equal(which_first(x %]between[% c(1, 1)), 1)

  expect_equal(which_first(x != 0.9), 1)
  expect_equal(which_first(x == 0.9), 0)
  expect_equal(which_first(x >= 0.9), 1)
  expect_equal(which_first(x <= 0.9), 0)
  expect_equal(which_first(x > 0.9), 1)
  expect_equal(which_first(x < 0.9), 0)
  expect_equal(which_first(x %in% 0.9), 0)
  expect_equal(which_first(x %between% c(0.9, 0.95)), 0)
  expect_equal(which_first(x %between% c(0.9, 1.95)), 1)
  expect_equal(which_first(x %(between)% c(0.9, 1)), 0)
  expect_equal(which_first(x %(between)% c(0.9, 1.2)), 1)
  expect_equal(which_first(x %]between[% c(0.9, 0.95)), 1)
  expect_equal(which_first(x %]between[% c(0.9, 1.95)), 0)

  x <- 1L
  expect_equal(which_first(x > 1), 0)
  expect_equal(which_first(x <= 1), 1)

  expect_equal(which_first(x != 1), 0)
  expect_equal(which_first(x == 1), 1)
  expect_equal(which_first(x >= 1), 1)
  expect_equal(which_first(x <= 1), 1)
  expect_equal(which_first(x > 1), 0)
  expect_equal(which_first(x < 1), 0)
  expect_equal(which_first(x %in% 1), 1)
  expect_equal(which_first(x %between% c(1, 2)), 1)
  expect_equal(which_first(x %(between)% c(1, 2)), 0)
  expect_equal(which_first(x %]between[% c(1, 2)), 1)

  x <- 0L
  expect_equal(which_first(x != 0.9), 1)
  expect_equal(which_first(x == 0.9), 0)
  expect_equal(which_first(x >= 0.9), 0)
  expect_equal(which_first(x <= 0.9), 1)
  expect_equal(which_first(x > 0.9), 0)
  expect_equal(which_first(x < 0.9), 1)
  expect_equal(which_first(x %in% 0.9), 0)
  expect_equal(which_first(x %in% c(0, 0.9)), 1)
  expect_equal(which_first(x %between% c(0, 0.9)), 1)
  expect_equal(which_first(x %(between)% c(0.9, 1)), 0)
  expect_equal(which_first(x %(between)% c(0.9, 0)), 0)
  expect_equal(which_first(x %]between[% c(0.9, 0.95)), 1)
  expect_equal(which_first(x %]between[% c(0.9, 1.95)), 1)
})

test_that("first_which", {
  expect_equal(first_which(c(FALSE, FALSE)), 0L)
  expect_equal(first_which(c(FALSE, TRUE)), 2L)
  expect_equal(first_which(c(TRUE, TRUE)), 1L)
  expect_equal(first_which(c(TRUE, NA)), 1L)
  expect_equal(first_which(c(NA, NA)), 0L)
})


test_that("which_first(<x> <o> <y>) lens equal", {
  x <- c(113L, 102L, 106L, 100L, 114L)
  y <- c(108L, 106L, 114L, 100L, 109L)
  y100 <- c(rep_len(y, length(-100:0)) + -100:0)

  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))
  z <- c(min(x), min(y))
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))
  z <- c(0L, 0L)
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))
  z <- c(1L, 0L)
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))


  x <- as.double(x)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x != x),
               first_which(x != x))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x >= x),
               first_which(x >= x))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x > x),
               first_which(x > x))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))
  z <- c(min(x), min(y))
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))

  y <- as.double(y)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))
  z <- c(min(x), min(y))
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))

  x <- as.integer(x)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))
  z <- c(min(x), min(y))
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))

  y <- y - 0.5
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))
  z <- c(min(x), min(y))
  expect_equal(which_first(x %(between)% z),
               first_which(x %(between)% z))
  expect_equal(which_first(x %]between[% z),
               first_which(x %]between[% z))

  x <- integer(0)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %in% y100),
               first_which(x %in% y100))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))

  x <- double(0)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))
  expect_equal(which_first(x %between% y[1:2]),
               first_which(x %between% y[1:2]))
  expect_equal(which_first(x %(between)% y[2:3]),
               first_which(x %(between)% y[2:3]))

  y <- integer(0)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))

  y <- double(0)
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))
  expect_equal(which_first(x %in% y),
               first_which(x %in% y))

  y <- Inf
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))

  y <- -Inf
  expect_equal(which_first(x != y),
               first_which(x != y))
  expect_equal(which_first(x == y),
               first_which(x == y))
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
  expect_equal(which_first(x > y),
               first_which(x > y))
  expect_equal(which_first(x < y),
               first_which(x < y))

  x <- c(5, -5, 2)
  y <- c(4, -4, 2)
  expect_equal(which_first(y != y), 0L)
  expect_equal(which_first(y >= x),
               first_which(y >= x))
  expect_equal(which_first(y > x),
               first_which(y > x))
  expect_equal(which_first(x < x),
               first_which(x < x))
  expect_equal(which_first(x <= y),
               first_which(x <= y))

  x <- c(5L, -5L, 2L)
  y <- c(4L, -4L, 2L)
  expect_equal(which_first(y != y), 0L)
  expect_equal(which_first(y >= x),
               first_which(y >= x))
  expect_equal(which_first(y > x),
               first_which(y > x))
  expect_equal(which_first(x < x),
               first_which(x < x))
  expect_equal(which_first(x <= y),
               first_which(x <= y))
})

test_that("which_first(lgl lgl)", {
  x <- c(TRUE, FALSE, TRUE)
  y <- c(FALSE, FALSE, TRUE)
  expect_equal(which_first(x == y), 2)
  expect_equal(which_first(x != y), 1)
  expect_equal(which_first(x >= y), 1)
  expect_equal(which_first(x <= y), 2)
  expect_equal(which_first(x > y), 1)
  expect_equal(which_first(x < y), 0)
  expect_equal(which_first(y >= x), 2)
  expect_equal(which_first(y < x), 1)

  B1 <- c(NA, TRUE)
  B2 <- c(TRUE, FALSE)
  B3 <- c(FALSE, TRUE)
  B4 <- c(NA, FALSE)
  T2 <- c(TRUE, TRUE)
  expect_equal(which_first(x %in% B1),
               first_which(x %in% B1))
  expect_equal(which_first(x %between% B2), 0)
  expect_equal(which_first(x %between% B3), 1)
  expect_equal(which_first(x %between% B4), 2)
  expect_equal(which_first(x %(between)% B3), 0)
  expect_equal(which_first(x %in% c(NA, TRUE, FALSE, TRUE)), 1)
  expect_equal(which_first(x %between% T2),
               first_which(and(x >= T2[1], x <= T2[2])))
  expect_equal(which_first(x %in% T2), which_first(x))
  T0 <- logical(0)
  expect_equal(which_first(x %in% T0),
               first_which(x %in% T0))
})

test_that("which_first internals", {
  expect_equal(do_which_first_notTRUE(c(FALSE)), 1L)
  expect_equal(do_which_first_notTRUE(c(TRUE)), 0L)
  expect_equal(do_which_first_notTRUE(c(TRUE, TRUE)), 0L)
  expect_equal(do_which_first_notTRUE(c(NA)), 1L)

  all_lgls <- c(TRUE, FALSE, NA)
  expect_equal(do_which_first_lgl_lgl_op(logical(0), logical(0), 1, TRUE), 0)
  expect_error(do_which_first_lgl_lgl_op(logical(5), logical(3), do_op2M("%between%")))
               # regex = "length.(2|two)")
  expect_equal(do_which_first_lgl_lgl_op(logical(11), all_lgls, do_op2M("%in%")), 1)
  expect_equal(do_which_first_lgl_lgl_op(logical(11), all_lgls, do_op2M("%in%"), TRUE), 11)
  expect_equal(do_which_first_lgl_lgl_op(logical(11), TRUE, do_op2M("%in%")), 0)
  expect_equal(do_which_first_lgl_lgl_op(logical(11), FALSE, do_op2M("%in%")), 1)
  expect_equal(do_which_first_lgl_lgl_op(logical(11), c(TRUE, FALSE), do_op2M("%between%")), 0)
  expect_equal(do_which_first_lgl_lgl_op(c(TRUE, NA), c(FALSE, NA), do_op2M("%in%")), 2)
  expect_equal(do_which_first_lgl_lgl_op(c(TRUE, TRUE), c(TRUE, FALSE), do_op2M("%between%")), 0)
  expect_equal(do_which_first_lgl_lgl_op(c(TRUE, TRUE), c(TRUE, TRUE), do_op2M("%in%")), 1)
  expect_equal(do_which_first_lgl_lgl_op(c(FALSE, TRUE), c(TRUE, TRUE), do_op2M("%in%")), 2)
  expect_equal(do_which_first_lgl_lgl_op(c(FALSE, FALSE), c(TRUE, TRUE), do_op2M("%between%")), 0)

})



test_that("which_firstNA", {
  expect_equal(which_firstNA(1:10), 0)
  expect_equal(which_firstNA(c(1:10, NA)), 11)
  expect_equal(which_firstNA(NA), 1L)
  expect_equal(which_firstNA(TRUE), 0)
  expect_equal(which_firstNA(c(NA, " ")), 1)
  expect_equal(which_firstNA(letters), 0)
  expect_equal(which_firstNA(raw(5)), 0)
  expect_equal(which_firstNA(NULL), 0)
  expect_equal(which_firstNA(c(NA, 0)), 1)
  expect_equal(which_firstNA(c(0, 0, NA, 0)), 3)
  expect_equal(which_firstNA(c(0, 0, NA, NA)), 3)
  expect_equal(which_firstNA(c(0, 0, 0)), 0)
  expect_equal(which_firstNA(c(NA_real_)), 1)
})

test_that("which_first(x %between% c(NA, 1))", {
  x <- c(1L, 5L, 3L, 10L, -1L)
  expect_equal(which_first(x %between% c(NA_integer_, NA_integer_)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 1L)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 0L)), length(x))
  expect_equal(which_first(x %between% c(NA_integer_, -2L)), 0)
  expect_equal(which_first(x %between% c(2L, NA_integer_)), 2)
  expect_equal(which_first(x %between% c(10L, NA_integer_)), 4)
  expect_equal(which_first(x %between% c(12L, NA_integer_)), 0)

  expect_equal(which_first(x %between% c(NA_integer_, NaN)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 1)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 0)), length(x))
  expect_equal(which_first(x %between% c(NA_integer_, -2)), 0)
  expect_equal(which_first(x %between% c(2, NA_integer_)), 2)
  expect_equal(which_first(x %between% c(10, NA_integer_)), 4)
  expect_equal(which_first(x %between% c(12, NA_integer_)), 0)

  x <- as.double(x)
  expect_equal(which_first(x %between% c(NA_integer_, NA_integer_)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 1L)), 1)
  expect_equal(which_first(x %between% c(NA_integer_, 0L)), length(x))
  expect_equal(which_first(x %between% c(NA_integer_, -2L)), 0)
  expect_equal(which_first(x %between% c(2L, NA_integer_)), 2)
  expect_equal(which_first(x %between% c(10L, NA_integer_)), 4)
  expect_equal(which_first(x %between% c(12L, NA_integer_)), 0)

  x <- c(0.1, 0.01, 0.25, 0.3, 0.9, 0.8)
  expect_equal(which_first(x %between% c(NA, NaN)), 1)
  expect_equal(which_first(x %between% c(NA, 1)), 1)
  expect_equal(which_first(x %between% c(NA, -1)), 0)
  expect_equal(which_first(x %between% c(NA, 0.09)), 2)
  expect_equal(which_first(x %between% c(0.5, NA)),
               first_which(x %between% c(0.5, NA)))
  expect_equal(which_first(x %between% c(0.9, NA)),
               first_which(x %between% c(0.9, NA)))
})

test_that("anyNA(x) implies which_first(x %in% c(NA, <rtype>))", {
  x_w_na <- c(1L, 5L, 3L, 10L, -1L, NA)
  expect_equal(which_first(x_w_na %in% c(NA, 11L)),
               first_which(x_w_na %in% c(NA, 11L)))
  expect_equal(which_first(x_w_na %in% c(NA, 10L)),
               first_which(x_w_na %in% c(NA, 10L)))
  expect_equal(which_first(x_w_na %in% c(NA, 11)),
               first_which(x_w_na %in% c(NA, 11)))
  expect_equal(which_first(x_w_na %in% c(NA, 10)),
               first_which(x_w_na %in% c(NA, 10)))
  x_w_na <- as.double(x_w_na)
  expect_equal(which_first(x_w_na %in% c(NA, 11L)),
               first_which(x_w_na %in% c(NA, 11L)))
  expect_equal(which_first(x_w_na %in% c(NA, 10L)),
               first_which(x_w_na %in% c(NA, 10L)))
  expect_equal(which_first(x_w_na %in% c(NA, 11)),
               first_which(x_w_na %in% c(NA, 11)))
  expect_equal(which_first(x_w_na %in% c(NA, 10)),
               first_which(x_w_na %in% c(NA, 10)))
})

test_that("which_first %between% error", {
  x <- 1:10
  expect_error(which_first(x %between% x))
  xd <- as.double(x)
  expect_error(which_first(x %between% xd))
  expect_error(which_first(xd %between% x))
  expect_error(which_first(xd %between% xd))
})

test_that("(between) and ]between[ with NA", {
  x <- c(11L, 1:5, 10L)
  expect_equal(which_first(x %(between)% c(NA, 2L)),
               which_first(x < 2L))
  expect_equal(which_first(x %(between)% c(5L, NA)),
               which_first(x > 5L))
  expect_equal(which_first(x %]between[% c(NA, 2L)),
               which_first(x >= 2L))
  expect_equal(which_first(x %]between[% c(5L, NA)),
               which_first(x <= 5L))

  expect_equal(first_which(x %(between)% c(NA, 2L)),
               which_first(x < 2L))
  expect_equal(first_which(x %(between)% c(5L, NA)),
               which_first(x > 5L))
  expect_equal(first_which(x %]between[% c(NA, 2L)),
               which_first(x >= 2L))
  expect_equal(first_which(x %]between[% c(5L, NA)),
               which_first(x <= 5L))

  expect_equal(which_first(x %(between)% c(NA, 2)),
               which_first(x < 2))
  expect_equal(which_first(x %(between)% c(5, NA)),
               which_first(x > 5))
  expect_equal(which_first(x %]between[% c(NA, 2)),
               which_first(x >= 2))
  expect_equal(which_first(x %]between[% c(5, NA)),
               which_first(x <= 5))

  expect_equal(first_which(x %(between)% c(NA, 2)),
               which_first(x < 2))
  expect_equal(first_which(x %(between)% c(5, NA)),
               which_first(x > 5))
  expect_equal(first_which(x %]between[% c(NA, 2)),
               which_first(x >= 2))
  expect_equal(first_which(x %]between[% c(5, NA)),
               which_first(x <= 5))


  x <- as.double(x)

  expect_equal(which_first(x %(between)% c(NA, 2L)),
               which_first(x < 2L))
  expect_equal(which_first(x %(between)% c(5L, NA)),
               which_first(x > 5L))
  expect_equal(which_first(x %]between[% c(NA, 2L)),
               which_first(x >= 2L))
  expect_equal(which_first(x %]between[% c(5L, NA)),
               which_first(x <= 5L))

  expect_equal(first_which(x %(between)% c(NA, 2L)),
               which_first(x < 2L))
  expect_equal(first_which(x %(between)% c(5L, NA)),
               which_first(x > 5L))
  expect_equal(first_which(x %]between[% c(NA, 2L)),
               which_first(x >= 2L))
  expect_equal(first_which(x %]between[% c(5L, NA)),
               which_first(x <= 5L))

  expect_equal(which_first(x %(between)% c(NA, 2)),
               which_first(x < 2))
  expect_equal(which_first(x %(between)% c(5, NA)),
               which_first(x > 5))
  expect_equal(which_first(x %]between[% c(NA, 2)),
               which_first(x >= 2))
  expect_equal(which_first(x %]between[% c(5, NA)),
               which_first(x <= 5))

  expect_equal(first_which(x %(between)% c(NA, 2)),
               which_first(x < 2))
  expect_equal(first_which(x %(between)% c(5, NA)),
               which_first(x > 5))
  expect_equal(first_which(x %]between[% c(NA, 2)),
               which_first(x >= 2))
  expect_equal(first_which(x %]between[% c(5, NA)),
               which_first(x <= 5))
})

test_that("lens 0", {
  xi <- c(1L, 2L)
  xd <- c(10, 20)
  x0 <- integer(0)
  x0d <- double(0)

  yi <- c(1L, 2L)
  yd <- c(10, 20)
  y0 <- integer(0)
  y0d <- double(0)

  expect_equal(which_first(xi == yi), first_which(xi == yi))
  expect_equal(which_first(xi == yd), first_which(xi == yd))
  expect_equal(which_first(xi == y0), first_which(xi == y0))
  expect_equal(which_first(xi == y0d), first_which(xi == y0d))

  expect_equal(which_first(xd == yi),
               first_which(xd == yi))
  expect_equal(which_first(xd == yd),
               first_which(xd == yd))
  expect_equal(which_first(xd == y0),
               first_which(xd == y0))
  expect_equal(which_first(xd == y0d),
               first_which(xd == y0d))

  expect_equal(which_first(x0 == yi),
               first_which(x0 == yi))
  expect_equal(which_first(x0 == yd),
               first_which(x0 == yd))
  expect_equal(which_first(x0 == y0),
               first_which(x0 == y0))
  expect_equal(which_first(x0 == y0d),
               first_which(x0 == y0d))

  expect_equal(which_first(x0d == yi),
               first_which(x0d == yi))
  expect_equal(which_first(x0d == yd),
               first_which(x0d == yd))
  expect_equal(which_first(x0d == y0),
               first_which(x0d == y0))
  expect_equal(which_first(x0d == y0d),
               first_which(x0d == y0d))


})

test_that("%in% with integers outside integer range", {
  s <- c(1L, 0L, -2L, 3L, NA_integer_)
  t <- c(-1e10, 1e10, 3)
  expect_equal(which_first(s %in% t), 4L)
  s <- as.double(s)
  expect_equal(which_first(s %in% t), 4L)
})

test_that("which_first bench mark", {
  skip_on_cran()
  skip_on_travis()
  skip_on_appveyor()
  skip_if_not_installed("bench")
  skip_on_ci()
  x <- double(1e8)
  which_first_time <- bench_system_time(which_first(x > 0))
  first_which_time <- bench_system_time(first_which(x > 0))
  expect_lt(which_first_time[2], 0.5 * first_which_time[2])
})


test_that("do_which_first_xi_ad", {
  x <- c(-.Machine$integer.max, -1L, -1L, .Machine$integer.max)
  w <- c(-3L, 1L, -2L)
  expect_equal(which_first(x >= Inf),
               first_which(x >= Inf))
  expect_equal(which_first(x <= Inf),
               first_which(x <= Inf))
  expect_equal(which_first(x >  Inf),
               first_which(x >  Inf))
  expect_equal(which_first(x <  Inf),
               first_which(x <  Inf))

  expect_equal(which_first(x >= -Inf),
               first_which(x >= -Inf))
  expect_equal(which_first(x <= -Inf),
               first_which(x <= -Inf))
  expect_equal(which_first(x >  -Inf),
               first_which(x >  -Inf))
  expect_equal(which_first(x <  -Inf),
               first_which(x <  -Inf))

  expect_equal(which_first(x > NA_real_),
               first_which(x > NA_real_))

  expect_equal(which_first(x != -.Machine$integer.max),
               first_which(x != -.Machine$integer.max))
  expect_equal(which_first(x == 0),
               first_which(x == 0))
  expect_equal(which_first(w > -10),
               first_which(w > -10))
  z3 <- integer(3)
  expect_equal(which_first(z3 != 0),
               first_which(z3 != 0))
  expect_equal(which_first(z3 > 0),
               first_which(z3 > 0))

  # Just outside int range
  expect_equal(which_first(z3 >= 2147483647.1),
               first_which(z3 >= 2147483647.1))
  expect_equal(which_first(z3 <= 2147483647.1),
               first_which(z3 <= 2147483647.1))
  expect_equal(which_first(z3 >= -2147483647.1),
               first_which(z3 >= -2147483647.1))
  expect_equal(which_first(z3 <= -2147483647.1),
               first_which(z3 <= -2147483647.1))

})

test_that("do_which_first_xd_yd", {
  x <- c(0, 2.9, 3.5)
  y <- x - 1.1
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))

  y <- x + 1
  expect_equal(which_first(x >= y),
               first_which(x >= y))
  expect_equal(which_first(x <= y),
               first_which(x <= y))

})

test_that("do_which_first_xi_yi", {
  x <- c(-1L, 2L, -3L)
  expect_equal(which_first(x >= -4L),
               first_which(x >= -4L))
  expect_equal(which_first(x > -4L),
               first_which(x > -4L))
  y <- as.integer(x) + 1L
  expect_equal(which_first(x >= y),
               which_first(x >= y))
  expect_equal(which_first(x > y),
               which_first(x > y))

  y <- as.integer(x) - 1L
  expect_equal(which_first(x <= y),
               which_first(x <= y))
  expect_equal(which_first(x < y),
               which_first(x < y))
})

test_that("do_which_first_xi_yd", {
  x <- rep(1L, 3)
  z <- as.double(x)
  expect_equal(which_first(x != z),
               first_which(x != z))
  z1.5 <- rep(1.5, 3)
  expect_equal(which_first(x >= z1.5),
               first_which(x >= z1.5))
  expect_equal(which_first(x > z1.5),
               first_which(x > z1.5))
  n1.5 <- rep(-1.5, 3)
  expect_equal(which_first(x <= n1.5),
               first_which(x <= n1.5))
  expect_equal(which_first(x < n1.5),
               first_which(x < n1.5))
})

test_that("do_which_first_xi_aii", {
  x <- rep(2L, 3)
  expect_equal(which_first(x %]between[% c(1L, 3L)),
               first_which(x %]between[% c(1L, 3L)))
})

test_that("do_which_first_xi_add", {
  x <- rep(2L, 3)
  expect_equal(which_first(x %]between[% c(1, 3)),
               first_which(x %]between[% c(1, 3)))
})

test_that("do_which_first_xd_add", {
  x <- rep(2, 3)
  expect_equal(which_first(x %]between[% c(1, 3)),
               first_which(x %]between[% c(1, 3)))
})

test_that("do_which_first_xi_ai", {
  x <- rep(7L, 4)
  expect_equal(which_first(x != 7L),
               first_which(x != 7L))
  expect_equal(which_first(x > 7L),
               first_which(x > 7L))
  expect_equal(which_first(x <= 6L),
               first_which(x <= 6L))
})

test_that("do_which_first_xi_ind", {
  x <- c(-400L, 4L, 5L, 2L)
  d0 <- double(0)
  tl <- seq(-50, 50, by = 0.5)
  t2 <- c(tl, NA)
  expect_equal(which_first(x %in% d0),
               first_which(x %in% d0))
  expect_equal(which_first(x %in% tl),
               first_which(x %in% tl))
  expect_equal(which_first(x %in% t2),
               first_which(x %in% t2))
  x <- c(-400.5, 400)
  expect_equal(which_first(x %in% t2),
               first_which(x %in% t2))
})

test_that("do_which_first_xd_ind", {
  x <- c(4, 5L, 2L)
  d0 <- double(0)
  expect_equal(which_first(x %in% d0),
               first_which(x %in% d0))
})







