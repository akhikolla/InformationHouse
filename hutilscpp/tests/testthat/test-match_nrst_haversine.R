context("test-match_nrst_haversine")

test_that("match_nrst_haversine error", {
  expect_error(match_nrst_haversine("a"))
  expect_error(match_nrst_haversine(1, "a"))
  expect_error(match_nrst_haversine(1, 1:2))
  expect_error(match_nrst_haversine(numeric(0), numeric(0)))
  expect_error(match_nrst_haversine(1, 1, "a"))
  expect_error(match_nrst_haversine(1, 1, 1, "a"))
  expect_error(match_nrst_haversine(1, 1, 1, 1:2))
  expect_error(match_nrst_haversine(1, 1, numeric(0), numeric(0)))

  lat2 <- seq(length.out = 5, -38, -37.8)
  lon2 <- rep(145, 5)

  lat1 <- c(-37.875, -37.91)
  lon1 <- c(144.96, 144.978)

  expect_error(match_nrst_haversine(lat1, lon1, lat2, lon2, excl_self = NULL),
               "excl_self")
  expect_error(match_nrst_haversine(lat1, lon1, lat2, lon2, as.data.table = NULL),
               "as.data.table")
  expect_error(match_nrst_haversine(lat1, lon1, lat2, lon2, close_enough = "a"),
               "`close_enough` was a character",
               fixed = TRUE)

  expect_warning(Table14_res <- match_nrst_haversine(lat1, lon1, lat2, lon2, Index = 1:4),
                 "not the same length as")
  expect_identical(Table14_res[[1L]], c(3L, 3L))


})

test_that("match_nrst_haversine works", {
  skip_if_not_installed("data.table")
  between <- data.table::between
  lat2 <- seq(from = -38, to = -37.8, length.out = 5)
  lon2 <- rep(145, 5)

  lat1 <- c(-37.855, -37.99)
  lon1 <- c(145.001, 144.98)

  oD <- match_nrst_haversine(lat1, lon1, lat2, lon2)
  oL <- match_nrst_haversine(lat1, lon1, lat2, lon2, as.data.table = FALSE)
  expect_equal(NROW(oD), 2L)
  expect_equal(NROW(oL), 2L)
  expect_identical(oL[[1L]], c(4L, 1L))
  expect_true(all(between(oL[[2L]], c(0.5, 2), c(0.6, 2.1))))


})

test_that("match_nrst_haversine Table", {
  lat2 <- seq(from = -38, to = -37.8, length.out = 6)
  lon2 <- rep(145, 6)

  lat1 <- c(-37.855, -37.99)
  lon1 <- c(145.001, 144.98)

  res <- match_nrst_haversine(lat1, lon1, lat2, lon2, Index = letters[1:6])
  expect_true(is.character(res[[1]]))
})

test_that("Suffixes", {
  res10 <- match_nrst_haversine(c(-37, -38), c(150, 149),
                                seq(-39, -36, length.out = 100),
                                seq(148, 151, length.out = 100))
  res10_RInf <- match_nrst_haversine(c(-37, -38), c(150, 149),
                                     seq(-39, -36, length.out = 100),
                                     seq(148, 151, length.out = 100))
  res10m <- match_nrst_haversine(c(-37, -38), c(150, 149),
                                 seq(-39, -36, length.out = 100),
                                 seq(148, 151, length.out = 100),
                                 .verify_box = FALSE,
                                 close_enough = "10m")
  res10m_in_km <-  match_nrst_haversine(c(-37, -38), c(150, 149),
                                        seq(-39, -36, length.out = 100),
                                        seq(148, 151, length.out = 100),
                                        .verify_box = FALSE,
                                        close_enough = "0.010km")
  expect_identical(res10, res10_RInf)
  expect_identical(res10, res10m)
  expect_identical(res10, res10m_in_km)
  res_equator_0m <- match_nrst_haversine(c(0, 0),
                                          c(149.95, 150.05),
                                          seq(-1, 1, length.out = 1000),
                                          seq(150, 150, length.out = 1000),
                                         cartesian_R = 0,
                                         .verify_box = FALSE,
                                          close_enough = 0)
  expect_identical(res_equator_0m[[1]], rep(500L, 2L))

  res_equator_10km <- match_nrst_haversine(c(0, 0),
                                           c(149.95, 150.05),
                                           seq(-1, 1, length.out = 1000),
                                           seq(150, 150, length.out = 1000),
                                           cartesian = 0,
                                           .verify_box = FALSE,
                                           close_enough = 10e3)
  expect_true(all(res_equator_10km[[2]] < 10))
  expect_true(all(res_equator_10km[[2]] > 5))
  res_equator_123km <- match_nrst_haversine(c(0, 0),
                                            c(149.95, 150.05),
                                            seq(-1, 1, length.out = 1000),
                                            rep(150, length.out = 1000),
                                            cartesian_R = 90,
                                            .verify_box = FALSE,
                                            close_enough = "53.0 km")
  expect_true(all(res_equator_123km[[2]] < 53))
  expect_true(all(res_equator_123km[[2]] > 52))

})
