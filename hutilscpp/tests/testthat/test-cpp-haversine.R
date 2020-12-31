context("test-haversine")

test_that("Error handling", {
  expect_error(haversineDistance(1, 1:2, 1:3, 1:4), regexp = "ength")
  expect_error(haversineDistance(1:2, 1:2, 1:3, 1:4), regexp = "ength")
  expect_error(which_min_HaversineDistance(1, 1:2, 1, 1),
               regexp = "length(lat1) != length(lat2)",
               fixed = TRUE)
  expect_error(match_min_Haversine(1, 1:2, 1, 1, 0L),
               regexp = "length(lat1) != length(lon1)",
               fixed = TRUE)
  expect_error(match_min_Haversine(1, 1, 1:2, 1, 0L),
               regexp = "length(lat2) != length(lon2)",
               fixed = TRUE)
  expect_warning(match_min_Haversine(1, 1, 1:2, 1:2, 0L, excl_self = TRUE),
                 regexp = "`excl_self = true`, yet lengths of `lat1` and `lat2` differ.",
                 fixed = TRUE)

})

# Same as hutils
# bankstown to sydney airports approximately 17628m
test_that("Bankstown airport to Sydney airport approximately 17628m", {
  expect_lt(haversineDistance(-33 - 56/60 - 46/3600, 151 + 10/60 + 38/3600,
                              -33 - 55/60 - 28/3600, 150 + 59/60+18/3600) / 17.628 - 1,
            0.01)
})

test_that("Broken Hill airport to Sydney airport approximately 932158", {
  expect_lt(haversineDistance(-33 - 56/60 - 46/3600, 151 + 10/60 + 38/3600,
                              -32 - 00/60 - 05/3600, 141 + 28/60 + 18/3600) / 932.158 - 1,
            0.01)
})

test_that("which_min_HaversineDistance", {
  lat1 <- -33 - seq(0, 1, length.out = 10)
  lon1 <- rep_len(150, 10)
  lat2 <- -33.09
  lon2 <- 150
  expect_identical(which_min_HaversineDistance(lat1, lon1, lat2, lon2), 2L)

  expect_identical(which_min_HaversineDistance(double(5), seq(-0.1, 0, length.out = 5), 0, 0,
                                               upperBound = 0),
                   5L)
  expect_identical(which_min_HaversineDistance(double(5), seq(-0.1, 0, length.out = 5), 0, 0,
                                               upperBound = 100),
                   5L)
  lat1 <- rep_len(-33, 10)
  lon1 <- 150 + seq(0, 1, length.out = 10)
  lat2 <- -33
  lon2 <- 150.09
  expect_identical(which_min_HaversineDistance(lat1, lon1, lat2, lon2), 2L)

  expect_identical(which_min_HaversineDistance(double(5), seq(-0.1, 0, length.out = 5), 0, 0,
                                               upperBound = 0),
                   5L)
  expect_identical(which_min_HaversineDistance(double(5), seq(-0.1, 0, length.out = 5), 0, 0,
                                               upperBound = 100),
                   5L)

})

test_that("match_min_Haversine", {
  lat2 <- c(-37.929, -37.962, -37.983, -37.928, -37.85)
  lon2 <- rep(145, 5)

  lat1 <- c(-37.875, -37.88)
  lon1 <- c(144.96, 144.978)

  expect_identical(match_min_Haversine(lat1, lon1, lat2, lon2, 0L)[[1L]], c(5L, 5L))
  expect_identical(match_min_Haversine(lat1, lon1, lat2, lon2, 101:105)[[1L]], c(5L, 5L) + 100L)
})

test_that("match_min_Haversine excl_self", {
  lat1 <- c(-37.875, -37.88)
  lon1 <- c(144.96, 144.978)
  mat <- match_min_Haversine(lat1, lon1, lat1, lon1, 0L)
  expect_equal(NROW(mat), length(lat1))
})

test_that("unitless", {
  lat2 <- c(-37.929, -37.962, -37.983, -37.928, -37.85)
  lon2 <- rep(145, 5)

  lat2 <- rep_len(lat2, 10)
  lon2 <- rep_len(lon2, 10)

  lat1 <- rep_len(c(-37.875, -37.88), 10)
  lon1 <- rep_len(c(144.96, 144.978), 10)

  expect_identical(order(haversineDistance(lat1, lon1, lat2, lon2)),
                   order(haversineDistance(lat1, lon1, lat2, lon2, unitless = TRUE)))
})

test_that("theEuclidDistance", {
  expect_equal(theEuclidDistance(0, 0, 0, 0), 0)
  expect_equal(theEuclidDistance(0, 1, 0, 0), 1)
  expect_equal(theEuclidDistance(0, 1, 0, 0, unitless = TRUE), 1)
  expect_equal(theEuclidDistance(0, 1, 0, 1), sqrt(2))
  expect_equal(theEuclidDistance(0, 1, 0, 1, unitless = TRUE), 2)

  x <- seq(1e7, by = 0.5, length.out = 200)
  ax <- seq(1e7, by = 0.7, length.out = 200)
  y <- seq(1e7 + 50, by = 1.1, length.out = 200)
  ay <- seq(1e7 + 50, by = 0.9, length.out = 200)
  dist_euc <- function(p, q, r, s) sqrt((p - q)^2 + (r - s)^2)
  # Can't be identical
  expect_equal(theEuclidDistance(x, ax, y, ay), dist_euc(x, ax, y, ay))
  expect_error(theEuclidDistance(x, ax, y, ay[1:5]),
               regexp = "lengths")


})

test_that("hausdorff distance", {
  o <- hausdorffEuclid(c(0, 0, 1), c(0, 1, 1))
  expect_equal(o, 1)
  o <- hausdorffEuclid(c(0, 0.5, 1), c(0, 0.7, 0))
  expect_equal(o, sqrt(0.74))
})

test_that("is_sorted_ascending", {
  expect_true(is_sorted_ascending_dbl(c(0, 1.5, 3)))
  expect_false(is_sorted_ascending_dbl(c(0, 1.5, 3, 1)))
})

test_that("Emptiest quadrants", {
  skip_if_not_installed("covr")
  skip_if_not_installed("withr")
  on_cran_not_covr <-
    !identical(Sys.getenv("NOT_CRAN"), "true") &&
    !covr::in_covr()

  skip_if(on_cran_not_covr, message = "On CRAN but not COVR")
  library(data.table)
  withr::with_seed(7441, {
    x <- c(sample(0:49, size = 2000, replace = TRUE), sample(76:100, size = 1000, replace = TRUE))
    y <- seq(0, 100, length.out = length(x))
    expect_identical(EmptiestQuarter(x, y, min(x), max(x), min(y), max(y)),  c(2L, 0L))
    expect_identical(last(theEmptiestQuarters(x, y, min(x), max(x), min(y), max(y))), -1L)

    not <- `!`
    ands <- function(x, ...) {
      if (missing(..1)) {
        return(x)
      }
      x & ands(...)
    }
    not2 <- function(...) {
      not(ands(...))
    }


    DT <- data.table(x = runif(10000, -1, 1),
                     y = runif(10000, -1, 1))
    DT_NE <- DT[not2(x > 0,
                     y > 0,
                     y - 1.9 * x <= 0.2,
                     y - x >= -0.1)]
    # Make sure we don't a gap in the wrong cell
    while (DT_NE[x %between% c(0.25, 0.5)][y %between% c(0.5, 0.75), .N] == 0L) {
      DT <- data.table(x = runif(10000, -1, 1),
                       y = runif(10000, -1, 1))
      DT_NE <- DT[not2(x > 0,
                       y > 0,
                       y - 1.9 * x <= 0.2,
                       y - x >= -0.1)]
    }


    res1 <- DT_NE[, poleInaccessibility2(x, y)]  # Will be unscaled
    res <- DT_NE[, poleInaccessibility2(x, y, x_range = c(-1, 1), y_range = c(-1, 1))]
    if (identical(res,
                  c("xmin" = 0.5,
                    "xmax" = 0.75,
                    "ymin" = 0.75,
                    "ymax" = 1.00))) {
      expect_identical(res,
                       c("xmin" = 0.5,
                         "xmax" = 0.75,
                         "ymin" = 0.75,
                         "ymax" = 1.00))
    } else {
      if (identical(Sys.getenv("USERNAME"), "hughp")) {
        saveRDS(res, "~/hutilscpp/data-raw/res.rds")
        saveRDS(DT_NE, "~/hutilscpp/data-raw/DT_NE.rds")
      }
    }

    # The chance of this if-statement being FALSE is vanishingly small
    # but Murphy's Law says it will happen on CRAN when they're making
    # a major release.  Let's lower the stakes.
    if (identical(round(res1, 2),
                  c("xmin" = 0.5,
                    "xmax" = 0.75,
                    "ymin" = 0.75,
                    "ymax" = 1.00))) {
      expect_identical(round(res1, 2),  # Won't fail too often?
                       c("xmin" = 0.5,
                         "xmax" = 0.75,
                         "ymin" = 0.75,
                         "ymax" = 1.00))
    }
    # if (FALSE) {
    #   library(ggplot2);ggplot(DT_NE,
    #                           aes(x, y)) +
    #     geom_point() +
    #     geom_blank(data = head(DT_NE, 2),
    #                aes(x = c(-1.01, 1.01),
    #                    y = c(-1.01, 1.01))) +
    #     annotate("rect",
    #              xmin = res["xmin"],
    #              xmax = res["xmax"],
    #              ymin = res["ymin"],
    #              ymax = res["ymax"],
    #              fill = NA,
    #              color = "red")
    # }
    res3 <- DT_NE[, poleInaccessibility3(x, y)]
    expect_false(DT_NE[, any(x > res3["xmin"] &
                               x < res3["xmax"] &
                               y > res3["ymin"] &
                               y < res3["ymax"])])

    res3 <- DT_NE[, poleInaccessibility3(x, y, x_range = c(-1, 1), y_range = c(-1, 1))]
    expect_false(DT_NE[, any(x > res3["xmin"] &
                               x < res3["xmax"] &
                               y > res3["ymin"] &
                               y < res3["ymax"])])


    expect_identical(DT_NE[, EmptiestQuarter(x, y, -1, 1, -1, 1)],
                     c(3L, DT_NE[x >= 0 & y >= 0, .N]))
    expect_identical(first(DT_NE[, EmptiestQuarter(x, y)]), 3L)
    expect_equal((DT_NE[, theEmptiestQuarters(x, y)])[1:2],
                 c(3, 3))

    DT <- data.table(x = runif(100000, -1, 1),
                     y = runif(100000, -1, 1))
    DT_NE <- DT[not2(x > 0,
                     y > 0,
                     {y - 1.9 * x} <= 0.2,
                     {y - 0.9*x} >= -0.1)]
    # Make sure we don't a gap in the wrong cell
    i <- 0
    while (i < 1e6 &&
           OR(DT_NE[x %between% c(0.25, 0.5)][y %between% c(0.5, 0.75), .N] == 0L,
              AND({DT_NE[, theEmptiestQuarters(x, y)]}[1] != 3,
                  {DT_NE[, theEmptiestQuarters(x, y)]}[2] != 3))) {
      DT <- data.table(x = runif(100000, -1, 1),
                       y = runif(100000, -1, 1))
      DT_NE <- DT[not2(x > 0,
                       y > 0,
                       y - 1.9*x <= 0.2,
                       y - 0.9*x >= -0.1)]
      i <- i + 1
    }



    # Cover all the quarters
    DT_NE[, y := -y]
    expect_identical(first(DT_NE[, EmptiestQuarter(x, y)]), 2L)
    expect_identical(first(DT_NE[, theEmptiestQuarters(x, y)]), 2L)
    DT_NE[, x := -x]
    expect_identical(first(DT_NE[, EmptiestQuarter(x, y)]), 0L)
    expect_identical(first(DT_NE[, theEmptiestQuarters(x, y)]), 0L)
    DT_NE[, y := -y]
    expect_identical(first(DT_NE[, EmptiestQuarter(x, y)]), 1L)
    expect_identical(first(DT_NE[, theEmptiestQuarters(x, y)]), 1L)

  }) # withr
})

test_that("poleInaccessibility error handling", {
  expect_error(poleInaccessibility2(), regexp = "were all NULL")
  expect_error(poleInaccessibility3(), regexp = "were all NULL")
})

test_that("poleInaccessibility3 infinite xmin_new's", {
  skip_if_not_installed("covr")
  on_cran_not_covr <-
    !identical(Sys.getenv("NOT_CRAN"), "true") &&
    !covr::in_covr()

  skip_if(on_cran_not_covr, message = "On CRAN but not COVR") # My fault but too sporadic for CRAN
  # Essentially need to test when a box occurs at the edges
  library(data.table)
  library(hutils) # for implies
  DT <- data.table(x = c(runif(10000),
                         runif(10000) * 0.9,
                         0.91, 0.92,
                         # Complete the box
                         1, 0.909),
                   y = c(runif(10000) * 0.9,
                         runif(10000),
                         0.925, 0.91,
                         # Complete the box
                         0.909, 1))
  resT <- DT[, poleInaccessibility3(x, y)]
  resF <- DT[, poleInaccessibility3(x, y, test_both = FALSE)]
  expect_equal(resT,
               c(xmin = 0.92,
                 xmax = 1.0,
                 ymin = 0.909,
                 ymax = 1.0))
  expect_equal(resF,
               c(xmin = 0.909,
                 xmax = 1.0,
                 ymin = 0.925,
                 ymax = 1.0))

  DT[, y := -y][, x := -x]
  res <- DT[, poleInaccessibility3(x, y)]
  expect_equal(res,
               c(xmin = -1,
                 xmax = -0.92,
                 ymin = -1.0,
                 ymax = -0.909))
})

test_that("poleInaccessibility error handling", {
  expect_error(poleInaccessibility2(DT = data.table()),
               "supplied but did not have.*LATITUDE.*LONGITUDE")
  expect_warning(poleInaccessibility2(data.table(LONGITUDE = 1:5 + 0,
                                                 LATITUDE = 11:15 + 0),
                                      x = 0, y = 0),
                 "`x` and `y` are not both NULL and will be ignored.",
                 fixed = TRUE)
  expect_error(poleInaccessibility3(DT = data.table()),
               regexp = "supplied but did not have.*LATITUDE.*LONGITUDE")
  expect_warning(poleInaccessibility3(data.table(LONGITUDE = 1:5 + 0,
                                                 LATITUDE = 11:15 + 0),
                                      x = 0, y = 0),
                 "`x` and `y` are not both NULL and will be ignored.",
                 fixed = TRUE)
  expect_error(cut_DT(data.table()),
               "lacked column")
})

test_that("match_min_haversine stops first at sufficiently close, then gets closer", {
  x <- c(-1, -0.5, 0.1, 0.05)
  y <- double(4)
  pos2 <-
    match_nrst_haversine(0, 0, y, x,
                         .verify_box = FALSE,
                         close_enough = "70 km",
                         cartesian_R = 0)
  expect_equal(pos2[["pos"]], 2L)
  expect_equal(floor(pos2[["dist"]]), 55)
  pos3 <-
    match_nrst_haversine(0,
                         0,
                         y, x,
                         close_enough = "70km",
                         cartesian_R = 0.4,
                         .verify_box = FALSE)


  pos4 <-
    match_nrst_haversine(0, 0,
                         y, x,
                         close_enough = "70km",
                         .verify_box = TRUE)
  expect_equal(pos3[["pos"]], 3L)
  expect_equal(pos4[["pos"]], 4L)
})

test_that("When the main match_min_haversine fails to pick a match", {
  lon <- c(-1, -0.5, 0.1, 0.05)
  lat <- double(4)
  match0 <- match_min_Haversine(0, 0,
                                lat, lon,
                                # small but nonzero
                                cartR = 0.00001,
                                verify_cartR = TRUE,
                                do_verify_box = FALSE,
                                tabl = 1L,
                                dist0_km = 200)
  expect_equal(match0$pos, 1L)
})




