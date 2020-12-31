context("which_last")

test_that("which_last works on integers", {
  x <- c(2L, 1L, 0L, 1L, 2L)
  expect_equal(which_last(x == 2), length(x))
  expect_equal(which_last(x != 2), length(x) - 1L)
  expect_equal(which_last(x == 1), length(x) - 1L)
  expect_equal(which_last(x != 1), length(x))
  expect_equal(which_last(x >= 1), length(x))
  expect_equal(which_last(x <= 1), length(x) - 1L)
  expect_equal(which_last(x > 1), length(x))
  expect_equal(which_last(x < 1), 3L)

  expect_equal(which_last(x == 2L), length(x))
  expect_equal(which_last(x != 2L), length(x) - 1L)
  expect_equal(which_last(x == 1L), length(x) - 1L)
  expect_equal(which_last(x != 1L), length(x))
  expect_equal(which_last(x >= 1L), length(x))
  expect_equal(which_last(x <= 1L), length(x) - 1L)
  expect_equal(which_last(x > 1L), length(x))
  expect_equal(which_last(x < 1L), 3L)

  y <- x + 0L
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- x + 1L
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- x - 1L
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- Inf
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y <- -Inf
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  # y <- NaN
  # expect_equal(which_last(x != y),
  #              last_which(x != y))
  # expect_equal(which_last(x == y),
  #              last_which(x == y))
  # expect_equal(which_last(x >= y),
  #              last_which(x >= y))
  # expect_equal(which_last(x <= y),
  #              last_which(x <= y))
  # expect_equal(which_last(x > y),
  #              last_which(x > y))
  # expect_equal(which_last(x < y),
  #              last_which(x < y))

  y <- max(x) + 1
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y <- min(x) - 1
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y <- as.double(x)
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y <- as.double(x) + 1
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y <- as.double(x) - 1
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))

  y[2] <- NA
  expect_equal(which_last(x != y),
               last_which(x != y))
  expect_equal(which_last(x == y),
               last_which(x == y))
  expect_equal(which_last(x >= y),
               last_which(x >= y))
  expect_equal(which_last(x <= y),
               last_which(x <= y))
  expect_equal(which_last(x > y),
               last_which(x > y))
  expect_equal(which_last(x < y),
               last_which(x < y))



})

test_that("which_last works on doubles", {
  x <- as.double(c(2L, 1L, 0L, 1L, 2L))
  expect_equal(which_last(x == 2), length(x))
  expect_equal(which_last(x != 2), length(x) - 1L)
  expect_equal(which_last(x == 1), length(x) - 1L)
  expect_equal(which_last(x != 1), length(x))
  expect_equal(which_last(x >= 1), length(x))
  expect_equal(which_last(x <= 1), length(x) - 1L)
  expect_equal(which_last(x > 1), length(x))
  expect_equal(which_last(x < 1), 3L)

  expect_equal(which_last(x == 2L), length(x))
  expect_equal(which_last(x != 2L), length(x) - 1L)
  expect_equal(which_last(x == 1L), length(x) - 1L)
  expect_equal(which_last(x != 1L), length(x))
  expect_equal(which_last(x >= 1L), length(x))
  expect_equal(which_last(x <= 1L), length(x) - 1L)
  expect_equal(which_last(x > 1L), length(x))
  expect_equal(which_last(x < 1L), 3L)

  y <- as.integer(x + 1L)
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- as.integer(x - 1L)
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- x + 1
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- x - 1
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  x0 <- double(9)
  expect_equal(which_last(x0 != 0),
               last_which(x0 != 0))
  expect_equal(which_last(x0 == 0.5),
               last_which(x0 == 0.5))
  expect_equal(which_last(x0 >= 0.5),
               last_which(x0 >= 0.5))
  expect_equal(which_last(x0 > 0.5),
               last_which(x0 > 0.5))
  expect_equal(which_last(x0 <= -0.5),
               last_which(x0 <= -0.5))
  expect_equal(which_last(x0 < -0.5),
               last_which(x0 < -0.5))

  x0 <- integer(9)
  expect_equal(which_last(x0 != 0),
               last_which(x0 != 0))
  expect_equal(which_last(x0 == 0.5),
               last_which(x0 == 0.5))
  expect_equal(which_last(x0 >= 0.5),
               last_which(x0 >= 0.5))
  expect_equal(which_last(x0 > 0.5),
               last_which(x0 > 0.5))
  expect_equal(which_last(x0 <= -0.5),
               last_which(x0 <= -0.5))
  expect_equal(which_last(x0 < -0.5),
               last_which(x0 < -0.5))

  expect_equal(which_last(x0 != 0L),
               last_which(x0 != 0L))
  expect_equal(which_last(x0 == 1L),
               last_which(x0 == 1L))
  expect_equal(which_last(x0 >= 1L),
               last_which(x0 >= 1L))
  expect_equal(which_last(x0 > 1L),
               last_which(x0 > 1L))
  expect_equal(which_last(x0 <= -1L),
               last_which(x0 <= -1L))
  expect_equal(which_last(x0 < -1L),
               last_which(x0 < -1L))


})

test_that("which_last on logicals", {
  x <- c(FALSE, TRUE, FALSE, TRUE)
  current_opt <- getOption("hutilscpp_suppressWarning")
  skip_if_not(is.null(current_opt))
  withr::with_options(list(hutilscpp_suppressWarning = FALSE),
                      {
                        expect_warning(which_last(x == NA),
                                       regexp = "appears to be logical NA")
                      })

  options(hutilscpp_suppressWarning = TRUE)

  expect_equal(which_last(x), 4L)
  expect_equal(which_last(x == TRUE), 4L)
  expect_equal(which_last(x == FALSE), 3L)
  expect_equal(which_last(x != TRUE), 3L)
  expect_equal(which_last(x != FALSE), 4L)
  expect_equal(which_last(x %in% c(TRUE, FALSE)), length(x))
  expect_equal(which_last(x %in% c(NA)), 0L)
  expect_equal(which_last(x %in% c(NA, FALSE)), 3L)
  expect_equal(which_last(x == NA), 0L)
  expect_equal(which_last(x != NA), length(x))

  axna <- c(FALSE, TRUE, FALSE, NA, FALSE)
  expect_equal(which_last(axna == TRUE), 2L)
  expect_equal(which_last(axna == FALSE), length(axna))
  expect_equal(which_last(axna == NA), length(axna) - 1L)
  expect_equal(which_last(axna != NA), length(axna))
  expect_equal(which_last(axna != TRUE), length(axna))
  options(hutilscpp_suppressWarning = NULL)
})

test_that("which_last is.na(lhs_eval)", {
  x <- c(TRUE, TRUE, FALSE)
  expect_equal(which_last(x == NA, suppressWarning = TRUE), 0L)
  expect_equal(which_last(x %in% c(NA, TRUE)), 2L)

  xna <- c(TRUE, FALSE, NA, TRUE, FALSE, NA)
  expect_equal(which_last(xna == NA, suppressWarning = TRUE), length(xna))
  expect_equal(which_last(xna != NA, suppressWarning = TRUE), length(xna) - 1L)
  expect_equal(which_last(xna != TRUE), length(xna) - 1L)
  expect_equal(which_last(xna == TRUE), 4L)
  expect_equal(which_last(xna %in% NA), length(xna))
  expect_equal(which_last(xna != FALSE, reverse = TRUE), 1L)


  noTRUES <- logical(10)
  expect_equal(which_last(noTRUES == TRUE), 0L)
})

test_that("which_last(x == y)", {
  x <- c(5:1, 0:4)
  y <- c(5:1, 0:4)
  w <- y + 1L
  expect_equal(which_last(x == y), length(x))
  expect_equal(which_last(x == w), 0L)
  expect_equal(which_last(x != y), 0L)
  expect_equal(which_last(x != w), length(x))
  z <- c(5:2, 0:5)
  expect_equal(which_last(x == z), 4L)
  expect_equal(which_last(x != z), length(x))
  expect_equal(which_last(x > z), 5L)

  xd <- as.double(x)
  yd <- as.double(y)
  wd <- as.double(w)
  expect_equal(which_last(xd == yd), length(x))
  expect_equal(which_last(xd == wd), 0L)
  expect_equal(which_last(xd != yd), 0L)
  expect_equal(which_last(xd != wd), length(x))
  zd <- as.double(z)
  expect_equal(which_last(xd == zd), 4L)
  expect_equal(which_last(xd != zd), length(x))
  expect_equal(which_last(xd > zd), 5L)

  i <- c(0.5, 0.4, 0.3, 0.5)
  j <- c(0.5, 0.4, 0.2 + 0.1, 0)
  expect_equal(which_last(i == j), 2)
})

test_that("which_last(all FALSE)", {
  x <- logical(10)
  expect_equal(which_last(x), 0L)
  expect_equal(which_last(x != FALSE), 0L)
  allTRUE <- rep(TRUE, 8)
  expect_equal(which_last(allTRUE != TRUE), 0L)
})

test_that("which_first long", {
  skip_on_travis()
  skip_on_appveyor()
  skip_if_covr()
  skip_if(.Machine$sizeof.pointer != 8)
  x <-
    tryCatch(allocate0_except(25e8, c(2e9, 23e8), c(1L, -1L)),
             error = function(e) {
               NULL
             })
  skip_if(is.null(x), message = "error during allocation")
  expect_equal(which_first(x > 2L), 0L)
  expect_equal(which_last(x == 0), length(x))
  expect_equal(which_first(x == 1), 2e9)
  expect_equal(which_first(x < 0), 23e8)
  expect_equal(which_first(x > 0), 2e9)
  expect_equal(which_first(x < -0.5), 23e8)
})

test_that("which_last internals", {
  expect_equal(do_which_last_notTRUE(c(FALSE)), 1L)
  expect_equal(do_which_last_notTRUE(c(TRUE)), 0L)
  expect_equal(do_which_last_notTRUE(c(TRUE, TRUE)), 0L)
  expect_equal(do_which_last_notTRUE(c(NA)), 1L)
  expect_equal(do_which_last_notFALSE(c(FALSE)), 0L)
  expect_equal(do_which_last_notFALSE(c(FALSE, FALSE)), 0L)
  expect_equal(do_which_last_notFALSE(c(TRUE)), 1L)
  expect_equal(do_which_last_notFALSE(c(NA)), 1L)
})

test_that("which_last( %between% )", {
  x <- c(1:10, -1L)
  expect_equal(which_last(x %between% c(1L, 10L)), 10L)
  expect_equal(which_last(x %between% c(-1L, -10L)), 0L)
  expect_equal(which_last(x %(between)% c(1L, 10L)),
               last_which(x %(between)% c(1L, 10L)))
  expect_equal(which_last(x %(between)% c(-1L, -10L)),
               last_which(x %(between)% c(-1L, -10L)))
  expect_equal(which_last(x %]between[% c(1L, 10L)), 11L)

  x <- as.double(x)
  expect_equal(which_last(x %between% c(1L, 10L)), 10L)
  expect_equal(which_last(x %]between[% c(1L, 10L)), 11L)

  x <- c(1:10, -1L)
  expect_equal(which_last(x %between% c(1, 10)), 10L)
  expect_equal(which_last(x %(between)% c(1, 10)),
               last_which(x %(between)% c(1, 10)))
  expect_equal(which_last(x %]between[% c(1, 10)), 11L)

  x <- as.double(x)
  expect_equal(which_last(x %between% c(1, 10)), 10L)
  expect_equal(which_last(x %(between)% c(1L, 10L)),
               last_which(x %(between)% c(1L, 10L)))
  expect_equal(which_last(x %(between)% c(1, 10)),
               last_which(x %(between)% c(1, 10)))
  expect_equal(which_last(x %]between[% c(1, 10)), 11L)

  x <- rep(5L, 11)
  expect_equal(which_last(x %]between[% c(4L, 6L)),
               last_which(x %]between[% c(4L, 6L)))
  expect_equal(which_last(x %]between[% c(4, 6)),
               last_which(x %]between[% c(4, 6)))
  x <- as.double(x)
  expect_equal(which_last(x %]between[% c(4L, 6L)),
               last_which(x %]between[% c(4L, 6L)))
  expect_equal(which_last(x %]between[% c(4, 6)),
               last_which(x %]between[% c(4, 6)))

})

test_that("last(lgl lgl)", {
  x <- c(TRUE, FALSE, TRUE)
  y <- c(FALSE, FALSE, TRUE)
  ftn <- c(FALSE, TRUE, NA)
  expect_equal(which_last(x == y), 3)
  expect_equal(which_last(x != y), 1)
  expect_equal(which_last(x >= y), 3)
  expect_equal(which_last(x <= y), 3)
  expect_equal(which_last(x > y), 1)
  expect_equal(which_last(x < y), 0)
  expect_equal(which_last(x %in% ftn), 3)
})

test_that("which_firstNA", {
  expect_equal(which_lastNA(1:10), 0)
  expect_equal(which_lastNA(c(1:10, NA)), 11)
  expect_equal(which_lastNA(NA), 1L)
  expect_equal(which_lastNA(TRUE), 0)
  expect_equal(which_lastNA(c(NA, " ")), 1)
  expect_equal(which_lastNA(letters), 0)
  expect_equal(which_lastNA(raw(5)), 0)
  expect_equal(which_lastNA(NULL), 0)
  expect_equal(which_lastNA(c(NA, 0)), 1)
  expect_equal(which_lastNA(c(0, 0, NA, 0)), 3)
  expect_equal(which_lastNA(c(0, 0, NA, NA)), 4)
  expect_equal(which_lastNA(c(0, 0, 0)), 0)
  expect_equal(which_lastNA(c(NA_real_)), 1)
})

test_that("last_which", {
  expect_equal(last_which(FALSE), 0L)
  expect_equal(last_which(logical(2)), 0L)
})

test_that("which_last(x %between% c(NA, 1))", {
  x <- c(1L, 5L, 3L, 10L, -1L)
  expect_equal(which_last(x %between% c(NA_integer_, NA_integer_)), length(x))
  expect_equal(which_last(x %between% c(NA_integer_, NaN)), length(x))

  x <- c(0.1, 0.01, 0.25, 0.3, 0.9, 0.8)
  expect_equal(which_last(x %between% c(NA, NA_integer_)), length(x))
  expect_equal(which_last(x %between% c(NA, NaN)), length(x))
  expect_equal(which_last(x %between% c(NA, 1)), length(x))
  expect_equal(which_last(x %between% c(NA, -1)), 0)
  expect_equal(which_last(x %between% c(NA, 0.09)), 2)
  expect_equal(which_last(x %between% c(0.5, NA)),
               last_which(x %between% c(0.5, NA)))
  expect_equal(which_last(x %between% c(0.9, NA)),
               last_which(x %between% c(0.9, NA)))
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

  expect_equal(which_last(xi == yi), last_which(xi == yi))
  expect_equal(which_last(xi == yd), last_which(xi == yd))
  expect_equal(which_last(xi == y0), last_which(xi == y0))
  expect_equal(which_last(xi == y0d), last_which(xi == y0d))

  expect_equal(which_last(xd == yi),
               last_which(xd == yi))
  expect_equal(which_last(xd == yd),
               last_which(xd == yd))
  expect_equal(which_last(xd == y0),
               last_which(xd == y0))
  expect_equal(which_last(xd == y0d),
               last_which(xd == y0d))

  expect_equal(which_last(x0 == yi),
               last_which(x0 == yi))
  expect_equal(which_last(x0 == yd),
               last_which(x0 == yd))
  expect_equal(which_last(x0 == y0),
               last_which(x0 == y0))
  expect_equal(which_last(x0 == y0d),
               last_which(x0 == y0d))

  expect_equal(which_last(x0d == yi),
               last_which(x0d == yi))
  expect_equal(which_last(x0d == yd),
               last_which(x0d == yd))
  expect_equal(which_last(x0d == y0),
               last_which(x0d == y0))
  expect_equal(which_last(x0d == y0d),
               last_which(x0d == y0d))
})


test_that("%in% with integers outside integer range", {
  s <- c(1L, 0L, -2L, 3L, NA_integer_)
  t <- c(-1e10, 1e10, 3)
  expect_equal(which_last(s %in% t), 4L)
  s <- as.double(s)
  expect_equal(which_last(s %in% t), 4L)
})

test_that("which_last(x %in% lgl)", {
  abc <- logical(0)
  lgl0 <- logical(0)
  ttt <- c(TRUE, TRUE, TRUE)
  ttf <- c(TRUE, TRUE, FALSE)
  ttn <- c(TRUE, TRUE, NA)
  tfn <- c(TRUE, FALSE, NA)
  fff <- c(FALSE, FALSE, FALSE)
  nana <- c(NA, NA, NA)



  expect_equal(which_last(abc %in% lgl0),
               last_which(abc %in% lgl0))
  expect_equal(which_last(ttt %in% fff),
               last_which(ttt %in% fff))
  expect_equal(which_last(ttt %in% ttt),
               last_which(ttt %in% ttt))
  expect_equal(which_last(fff %in% fff),
               last_which(fff %in% fff))
  expect_equal(which_last(ttf %in% lgl0),
               last_which(ttf %in% lgl0))
})

test_that("do_which_last_xi_ini", {
  x50 <- c(0L, 40L, 53:99)
  x100 <- c(0L, 100:400L)
  x <- c(1:40, 50:550, 40L)
  expect_equal(which_last(x %in% x50),
               last_which(x %in% x50))
  expect_equal(which_last(x %in% x100),
               last_which(x %in% x100))

  x <- -x
  expect_equal(which_last(x %in% x50),
               last_which(x %in% x50))
  expect_equal(which_last(x %in% x100),
               last_which(x %in% x100))


  x50 <- as.double(x50)
  x100 <- as.double(x100)
  expect_equal(which_last(x %in% x50),
               last_which(x %in% x50))
  expect_equal(which_last(x %in% x100),
               last_which(x %in% x100))

  x <- -x
  expect_equal(which_last(x %in% x50),
               last_which(x %in% x50))
  expect_equal(which_last(x %in% x100),
               last_which(x %in% x100))

  x50 <- c(NA, x50)
  x100 <- c(NA, x100)
  z <- c(NA, x)
  expect_equal(which_last(z %in% x50),
               last_which(z %in% x50))
  expect_equal(which_last(z %in% x100),
               last_which(z %in% x100))


})

test_that("do_which_last_xd_yi", {
  x <- double(3)
  expect_equal(which_last(x != 0L),
               last_which(x != 0L))
  expect_equal(which_last(x == 0L),
               last_which(x == 0L))

  x <- c(double(10), 101, 88)
  y <- as.integer(x)

  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- c(-5:4, 101, 89)

  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- c(-5:4, 100, 88)
  expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

  y <- c(-5:4, 100, NA)
  # expect_equal(which_last(x != y), last_which(x != y))
  expect_equal(which_last(x == y), last_which(x == y))
  expect_equal(which_last(x >= y), last_which(x >= y))
  expect_equal(which_last(x <= y), last_which(x <= y))
  expect_equal(which_last(x > y), last_which(x > y))
  expect_equal(which_last(x < y), last_which(x < y))

})

test_that("do_which_last_xd_ind", {
  x <- c(-400, 4, 5, NA, 2)
  d0 <- double(0)
  tl <- seq(-50, 50, by = 0.5)
  t2 <- c(tl, NA)
  expect_equal(which_last(x %in% d0),
               last_which(x %in% d0))
  expect_equal(which_last(x %in% tl),
               last_which(x %in% tl))
  expect_equal(which_last(x %in% t2),
               last_which(x %in% t2))
  x <- c(-400.5, 400)
  expect_equal(which_last(x %in% t2),
               last_which(x %in% t2))
  x <- c(x, NA, x)
  expect_equal(which_last(x %in% t2),
               last_which(x %in% t2))
  x <- c(NA, x, x)
  expect_equal(which_last(x %in% t2),
               last_which(x %in% t2))
  t2 <- c(head(t2, 98), NA)
  expect_equal(which_last(x %in% t2),
               last_which(x %in% t2))


})

test_that("do_which_last_xi_add", {
  x <- rep(5L, 3)
  expect_equal(which_last(x %between% c(6, 7)),
               last_which(x %between% c(6, 7)))
  expect_equal(which_last(x %(between)% c(6, 8)),
               last_which(x %(between)% c(6, 8)))
})

test_that("do_which_last_xi_add", {
  x <- rep(5, 3)
  expect_equal(which_last(x %between% c(6, 7)),
               last_which(x %between% c(6, 7)))
  expect_equal(which_last(x %(between)% c(6, 8)),
               last_which(x %(between)% c(6, 8)))
})

test_that("do_which_last_xi_ad", {

  z3 <- integer(3)
  # Just outside int range
  expect_equal(which_last(z3 >= 2147483647.1),
               last_which(z3 >= 2147483647.1))
  expect_equal(which_last(z3 <= 2147483647.1),
               last_which(z3 <= 2147483647.1))
  expect_equal(which_last(z3 >= -2147483647.1),
               last_which(z3 >= -2147483647.1))
  expect_equal(which_last(z3 <= -2147483647.1),
               last_which(z3 <= -2147483647.1))
})




