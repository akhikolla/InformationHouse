
test_that("wkt_set_srid works", {
  expect_identical(wkt_set_srid(character(0), 1234), character(0))
  expect_identical(wkt_set_srid(NA_character_, 1234), NA_character_)
  expect_identical(
    wkt_set_srid("POINT (30 10)", 1234),
    "SRID=1234;POINT (30 10)"
  )

  expect_identical(
    wkt_set_srid("POINT (30 10)", c(1234, 5678)),
    c("SRID=1234;POINT (30 10)", "SRID=5678;POINT (30 10)")
  )

  expect_identical(
    wkt_set_srid(c("POINT (30 10)", "POINT (10 10)"), c(1234, 5678)),
    c("SRID=1234;POINT (30 10)", "SRID=5678;POINT (10 10)")
  )

  expect_identical(
    wkt_set_srid(c("POINT (30 10)", "POINT (10 10)"), c(1234)),
    c("SRID=1234;POINT (30 10)", "SRID=1234;POINT (10 10)")
  )
})

test_that("wkb_set_srid works", {
  expect_identical(wkb_set_srid(list(), 1234), list())
  expect_identical(wkb_set_srid(list(NULL), 1234), list(NULL))
  expect_identical(
    wkb_set_srid(as_wkb("POINT (30 10)"), 1234),
    unclass(as_wkb("SRID=1234;POINT (30 10)"))
  )
})

test_that("wksxp_set_srid works", {
  expect_identical(wksxp_set_srid(list(), 1234), list())
  expect_identical(wksxp_set_srid(list(NULL), 1234), list(NULL))
  expect_identical(
    wksxp_set_srid(as_wksxp("POINT (30 10)"), 1234),
    unclass(as_wksxp("SRID=1234;POINT (30 10)"))
  )
})

test_that("wkt_set_z works", {
  expect_identical(wkt_set_z(character(0), 1234), character(0))
  expect_identical(wkt_set_z(NA_character_, 1234), NA_character_)
  expect_identical(
    wkt_set_z("POINT (30 10)", 1234),
    "POINT Z (30 10 1234)"
  )

  expect_identical(
    wkt_set_z("POINT (30 10)", c(1234, 5678)),
    c("POINT Z (30 10 1234)", "POINT Z (30 10 5678)")
  )

  expect_identical(
    wkt_set_z(c("POINT (30 10)", "POINT (10 10)"), c(1234, 5678)),
    c("POINT Z (30 10 1234)", "POINT Z (10 10 5678)")
  )

  expect_identical(
    wkt_set_z(c("POINT (30 10)", "POINT (10 10)"), c(1234)),
    c("POINT Z (30 10 1234)", "POINT Z (10 10 1234)")
  )
})

test_that("wkb_set_z works", {
  expect_identical(wkb_set_z(list(), 1234), list())
  expect_identical(wkb_set_z(list(NULL), 1234), list(NULL))
  expect_identical(
    wkb_set_z(as_wkb("POINT (30 10)"), 1234),
    unclass(as_wkb("POINT Z (30 10 1234)"))
  )
})

test_that("wksxp_set_z works", {
  expect_identical(wksxp_set_z(list(), 1234), list())
  expect_identical(wksxp_set_z(list(NULL), 1234), list(NULL))
  expect_identical(
    wksxp_set_z(as_wksxp("POINT (30 10)"), 1234),
    unclass(as_wksxp("POINT Z (30 10 1234)"))
  )
})

test_that("wkt_transform works", {
  tx <- 3
  ty <- 7
  t_test_trans <- matrix(c(1, 0, 0, 0, 1, 0, tx, ty, 1), ncol = 3)

  t_test_rotate45 <- matrix(
    c(0.707106781186548, -0.707106781186547, 0, 0.707106781186547,
      0.707106781186548, 0, 0, 0, 1),
    ncol = 3
  )

  expect_identical(wkt_transform("POINT (0 0)", t_test_trans), "POINT (3 7)")
  expect_identical(wkt_transform("POINT (0 0)", t_test_rotate45), "POINT (0 0)")
  expect_identical(
    wkt_transform(sprintf("POINT (%s 0)", sqrt(2)), t_test_rotate45, precision = 10),
    "POINT (1 -1)"
  )

  expect_identical(
    wkt_transform(
      sprintf("POINT (%s 0)", sqrt(2)),
      t_test_trans %*% t_test_rotate45,
      precision = 10
    ),
    "POINT (4 6)"
  )
})

test_that("wkb_transform works", {
  tx <- 3
  ty <- 7
  t_test_trans <- matrix(c(1, 0, 0, 0, 1, 0, tx, ty, 1), ncol = 3)
  expect_identical(
    wkb_transform(as_wkb("POINT (0 0)"), t_test_trans),
    wkt_translate_wkb("POINT (3 7)")
  )
})

test_that("wksxp_transform works", {
  tx <- 3
  ty <- 7
  t_test_trans <- matrix(c(1, 0, 0, 0, 1, 0, tx, ty, 1), ncol = 3)
  expect_identical(
    wksxp_transform(as_wksxp("POINT (0 0)"), t_test_trans),
    wkt_translate_wksxp("POINT (3 7)")
  )
})
