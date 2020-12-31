test_that("which2 works", {
  frog <- 1:5
  pear <- 101:105

  which2 <- which2_Year

  expect_identical(which2(frog >= 4L, pear <= 105L),
                   which(frog >= 4L & pear <= 105L))
  expect_identical(which2(frog == 4L, pear <= 105L),
                   which(frog == 4L & pear <= 105L))
  expect_identical(which2(frog == 4L, pear >= 105L),
                   which(frog == 4L & pear >= 105L))
  expect_identical(which2(frog == 4L, pear < 105L),
                   which(frog == 4L & pear < 105L))

  coffee <- c(integer(10), c(1L, 1L, 2L, 2L, 8L, 3L))
  wrench <- c(integer(10), c(3L, -2L, 4L, 9L, 6L, 3L))

  expect_identical(which2(coffee == 0L, wrench != 1L),
                   which(coffee == 0L & wrench != 1L))
  expect_identical(which2(coffee == 0L, wrench == 1L),
                   which(coffee == 0L & wrench == 1L))
  expect_identical(which2(coffee == 0L, wrench >= 1L),
                   which(coffee == 0L & wrench >= 1L))
  expect_identical(which2(coffee == 0L, wrench <= 1L),
                   which(coffee == 0L & wrench <= 1L))
  expect_identical(which2(coffee == 0L, wrench < 1L),
                   which(coffee == 0L & wrench < 1L))
  expect_identical(which2(coffee != 0L, wrench != 1L),
                   which(coffee != 0L & wrench != 1L))
  expect_identical(which2(coffee != 0L, wrench != 1L),
                   which(coffee != 0L & wrench != 1L))

  expect_identical(which2(coffee != 0L , wrench != 1L), which(coffee != 0L & wrench != 1L))
  expect_identical(which2(coffee != 0L , wrench == 1L), which(coffee != 0L & wrench == 1L))
  expect_identical(which2(coffee != 0L , wrench >= 1L), which(coffee != 0L & wrench >= 1L))
  expect_identical(which2(coffee != 0L , wrench <= 1L), which(coffee != 0L & wrench <= 1L))
  expect_identical(which2(coffee != 0L , wrench > 1L), which(coffee != 0L & wrench > 1L))
  expect_identical(which2(coffee != 0L , wrench < 1L), which(coffee != 0L & wrench < 1L))

  expect_identical(which2(coffee == 0L , wrench != 1L), which(coffee == 0L & wrench != 1L))
  expect_identical(which2(coffee == 0L , wrench == 1L), which(coffee == 0L & wrench == 1L))
  expect_identical(which2(coffee == 0L , wrench >= 1L), which(coffee == 0L & wrench >= 1L))
  expect_identical(which2(coffee == 0L , wrench <= 1L), which(coffee == 0L & wrench <= 1L))
  expect_identical(which2(coffee == 0L , wrench > 1L), which(coffee == 0L & wrench > 1L))
  expect_identical(which2(coffee == 0L , wrench < 1L), which(coffee == 0L & wrench < 1L))

  expect_identical(which2(coffee >= 0L , wrench != 1L), which(coffee >= 0L & wrench != 1L))
  expect_identical(which2(coffee >= 0L , wrench == 1L), which(coffee >= 0L & wrench == 1L))
  expect_identical(which2(coffee >= 0L , wrench >= 1L), which(coffee >= 0L & wrench >= 1L))
  expect_identical(which2(coffee >= 0L , wrench <= 1L), which(coffee >= 0L & wrench <= 1L))
  expect_identical(which2(coffee >= 0L , wrench > 1L), which(coffee >= 0L & wrench > 1L))
  expect_identical(which2(coffee >= 0L , wrench < 1L), which(coffee >= 0L & wrench < 1L))

  expect_identical(which2(coffee <= 0L , wrench != 1L), which(coffee <= 0L & wrench != 1L))
  expect_identical(which2(coffee <= 0L , wrench == 1L), which(coffee <= 0L & wrench == 1L))
  expect_identical(which2(coffee <= 0L , wrench >= 1L), which(coffee <= 0L & wrench >= 1L))
  expect_identical(which2(coffee <= 0L , wrench <= 1L), which(coffee <= 0L & wrench <= 1L))
  expect_identical(which2(coffee <= 0L , wrench > 1L), which(coffee <= 0L & wrench > 1L))
  expect_identical(which2(coffee <= 0L , wrench < 1L), which(coffee <= 0L & wrench < 1L))

  expect_identical(which2(coffee > 0L , wrench != 1L), which(coffee > 0L & wrench != 1L))
  expect_identical(which2(coffee > 0L , wrench == 1L), which(coffee > 0L & wrench == 1L))
  expect_identical(which2(coffee > 0L , wrench >= 1L), which(coffee > 0L & wrench >= 1L))
  expect_identical(which2(coffee > 0L , wrench <= 1L), which(coffee > 0L & wrench <= 1L))
  expect_identical(which2(coffee > 0L , wrench > 1L), which(coffee > 0L & wrench > 1L))
  expect_identical(which2(coffee > 0L , wrench < 1L), which(coffee > 0L & wrench < 1L))

  expect_identical(which2(coffee < 0L , wrench != 1L), which(coffee < 0L & wrench != 1L))
  expect_identical(which2(coffee < 0L , wrench == 1L), which(coffee < 0L & wrench == 1L))
  expect_identical(which2(coffee < 0L , wrench >= 1L), which(coffee < 0L & wrench >= 1L))
  expect_identical(which2(coffee < 0L , wrench <= 1L), which(coffee < 0L & wrench <= 1L))
  expect_identical(which2(coffee < 0L , wrench > 1L), which(coffee < 0L & wrench > 1L))
  expect_identical(which2(coffee < 0L , wrench < 1L), which(coffee < 0L & wrench < 1L))

  expect_identical(which2(coffee != 8L , wrench != 1L), which(coffee != 8L & wrench != 1L))
  expect_identical(which2(coffee != 8L , wrench == 1L), which(coffee != 8L & wrench == 1L))
  expect_identical(which2(coffee != 8L , wrench >= 1L), which(coffee != 8L & wrench >= 1L))
  expect_identical(which2(coffee != 8L , wrench <= 1L), which(coffee != 8L & wrench <= 1L))
  expect_identical(which2(coffee != 8L , wrench > 1L), which(coffee != 8L & wrench > 1L))
  expect_identical(which2(coffee != 8L , wrench < 1L), which(coffee != 8L & wrench < 1L))

  expect_identical(which2(coffee == 8L , wrench != 1L), which(coffee == 8L & wrench != 1L))
  expect_identical(which2(coffee == 8L , wrench == 1L), which(coffee == 8L & wrench == 1L))
  expect_identical(which2(coffee == 8L , wrench >= 1L), which(coffee == 8L & wrench >= 1L))
  expect_identical(which2(coffee == 8L , wrench <= 1L), which(coffee == 8L & wrench <= 1L))
  expect_identical(which2(coffee == 8L , wrench > 1L), which(coffee == 8L & wrench > 1L))
  expect_identical(which2(coffee == 8L , wrench < 1L), which(coffee == 8L & wrench < 1L))

  expect_identical(which2(coffee >= 8L , wrench != 1L), which(coffee >= 8L & wrench != 1L))
  expect_identical(which2(coffee >= 8L , wrench == 1L), which(coffee >= 8L & wrench == 1L))
  expect_identical(which2(coffee >= 8L , wrench >= 1L), which(coffee >= 8L & wrench >= 1L))
  expect_identical(which2(coffee >= 8L , wrench <= 1L), which(coffee >= 8L & wrench <= 1L))
  expect_identical(which2(coffee >= 8L , wrench > 1L), which(coffee >= 8L & wrench > 1L))
  expect_identical(which2(coffee >= 8L , wrench < 1L), which(coffee >= 8L & wrench < 1L))

  expect_identical(which2(coffee <= 8L , wrench != 1L), which(coffee <= 8L & wrench != 1L))
  expect_identical(which2(coffee <= 8L , wrench == 1L), which(coffee <= 8L & wrench == 1L))
  expect_identical(which2(coffee <= 8L , wrench >= 1L), which(coffee <= 8L & wrench >= 1L))
  expect_identical(which2(coffee <= 8L , wrench <= 1L), which(coffee <= 8L & wrench <= 1L))
  expect_identical(which2(coffee <= 8L , wrench > 1L), which(coffee <= 8L & wrench > 1L))
  expect_identical(which2(coffee <= 8L , wrench < 1L), which(coffee <= 8L & wrench < 1L))

  expect_identical(which2(coffee > 8L , wrench != 1L), which(coffee > 8L & wrench != 1L))
  expect_identical(which2(coffee > 8L , wrench == 1L), which(coffee > 8L & wrench == 1L))
  expect_identical(which2(coffee > 8L , wrench >= 1L), which(coffee > 8L & wrench >= 1L))
  expect_identical(which2(coffee > 8L , wrench <= 1L), which(coffee > 8L & wrench <= 1L))
  expect_identical(which2(coffee > 8L , wrench > 1L), which(coffee > 8L & wrench > 1L))
  expect_identical(which2(coffee > 8L , wrench < 1L), which(coffee > 8L & wrench < 1L))

  expect_identical(which2(coffee < 8L , wrench != 1L), which(coffee < 8L & wrench != 1L))
  expect_identical(which2(coffee < 8L , wrench == 1L), which(coffee < 8L & wrench == 1L))
  expect_identical(which2(coffee < 8L , wrench >= 1L), which(coffee < 8L & wrench >= 1L))
  expect_identical(which2(coffee < 8L , wrench <= 1L), which(coffee < 8L & wrench <= 1L))
  expect_identical(which2(coffee < 8L , wrench > 1L), which(coffee < 8L & wrench > 1L))
  expect_identical(which2(coffee < 8L , wrench < 1L), which(coffee < 8L & wrench < 1L))
})

test_that("which2_year", {
  xya <- xyb <- xyc <- 1:100
  expect_identical(which2_Year(xya >= 10L, xyb >= 20L, xyc, 100L), 100L)
})

test_that("Warnings errors", {
  ab <- c(1:10, 11:15)
  bc <- c(11:15, 1:10)
  expect_warning(which2_Year(bc <= 10L, ab >= 3L, yr = 2004L),
                 regexp = "has been provided")
  expect_error(which2_Year(bc <= 10L, ab >= 3L, Year = ab, yr = NULL),
               regexp = "provided but `yr` remained NULL",
               fixed = TRUE)
  expect_error(which2_Year(bc <= 10L, ab >= 3L, Year = ab + 0, yr = 0L),
               regexp = "(strictly) type integer",
               fixed = TRUE)
  expect_error(which2_Year(bc <= 10L, ab >= 3L, Year = 1:10, yr = 0L),
               regexp = "evaluates to")

  expect_error(do_which2_yr(Year = integer(10), yr = 0L, consider_yr = FALSE,
                            x = integer(7), xa = 1L, eqx = FALSE, gtx = FALSE, ltx = FALSE,
                            y = integer(11), ya = 2L, eqy = FALSE, gty = FALSE, lty = FALSE),
               "lengths differ")
  expect_error(do_which2_yr(Year = integer(11), yr = 0L, consider_yr = FALSE,
                            x = integer(11), xa = 1L, eqx = FALSE, gtx = TRUE, ltx = TRUE,
                            y = integer(11), ya = 2L, eqy = FALSE, gty = FALSE, lty = FALSE),
               "both > and < requested")
  expect_error(do_which2_yr(Year = integer(11), yr = 0L, consider_yr = FALSE,
                            x = integer(11), xa = 1L, eqx = FALSE, gtx = FALSE, ltx = FALSE,
                            y = integer(11), ya = 2L, eqy = FALSE, gty = TRUE, lty = TRUE),
               "both > and < requested")
  expect_error(do_which2_yr(Year = integer(10), yr = 0L, consider_yr = TRUE,
                            x = integer(11), xa = 1L, eqx = FALSE, gtx = FALSE, ltx = FALSE,
                            y = integer(11), ya = 2L, eqy = FALSE, gty = FALSE, lty = FALSE),
               "lengths differ (Year)",
               fixed = TRUE)
})

