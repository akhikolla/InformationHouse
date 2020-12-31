
test_that("wkt_has_missing() / wkt_is_finite() works", {
  expect_identical(wkt_has_missing(character(0)), logical(0))
  expect_false(wkt_has_missing("POINT EMPTY"))
  expect_false(wkt_has_missing("POINT (0 1)"))
  expect_true(wkt_has_missing("POINT (nan nan)"))
  expect_false(wkt_has_missing("POINT (inf inf)"))

  expect_true(wkt_is_finite("POINT EMPTY"))
  expect_true(wkt_is_finite("POINT (0 1)"))
  expect_false(wkt_is_finite("POINT (nan nan)"))
  expect_false(wkt_is_finite("POINT (inf inf)"))
})

test_that("wkb_has_missing() / wkb_is_finite() works", {
  expect_identical(wkb_has_missing(list()), logical(0))
  # WKB can't do empty point without missing coords
  expect_false(wkb_has_missing(as_wkb("LINESTRING EMPTY")))
  expect_false(wkb_has_missing(as_wkb("POINT (0 1)")))
  expect_true(wkb_has_missing(as_wkb("POINT (nan nan)")))
  expect_false(wkb_has_missing(as_wkb("POINT (inf inf)")))

  expect_identical(wkb_is_finite(list()), logical(0))
  expect_true(wkb_is_finite(as_wkb("LINESTRING EMPTY")))
  expect_true(wkb_is_finite(as_wkb("POINT (0 1)")))
  expect_false(wkb_is_finite(as_wkb("POINT (nan nan)")))
  expect_false(wkb_is_finite(as_wkb("POINT (inf inf)")))
})

test_that("wksxp_has_missing() / wksxp_is_finite() works", {
  expect_identical(wksxp_has_missing(list()), logical(0))
  expect_false(wksxp_has_missing(as_wksxp("POINT EMPTY")))
  expect_false(wksxp_has_missing(as_wksxp("POINT (0 1)")))
  expect_true(wksxp_has_missing(as_wksxp("POINT (nan nan)")))
  expect_false(wksxp_has_missing(as_wksxp("POINT (inf inf)")))

  expect_identical(wksxp_is_finite(list()), logical(0))
  expect_true(wksxp_is_finite(as_wksxp("POINT EMPTY")))
  expect_true(wksxp_is_finite(as_wksxp("POINT (0 1)")))
  expect_false(wksxp_is_finite(as_wksxp("POINT (nan nan)")))
  expect_false(wksxp_is_finite(as_wksxp("POINT (inf inf)")))
})
