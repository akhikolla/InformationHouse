
test_that("wk*_unnest() works", {
  expect_identical(wkt_unnest(NA_character_), structure(NA_character_, lengths = 1L))
  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)",
      keep_multi = FALSE, keep_empty = FALSE, max_depth = 2
    ),
    structure(c("POINT (30 10)", "POINT (10 10)", "LINESTRING (0 0, 1 1)"), lengths = 3L)
  )

  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)",
      keep_multi = FALSE, keep_empty = TRUE, max_depth = 2
    ),
    structure(
      c("POINT (30 10)", "POINT (10 10)", "LINESTRING (0 0, 1 1)", "GEOMETRYCOLLECTION EMPTY"),
      lengths = 4L
    )
  )

  expect_identical(
    wkt_unnest(
      "SRID=12;GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)",
      keep_multi = FALSE, max_depth = 2
    ),
    structure(
      c("SRID=12;POINT (30 10)", "SRID=12;POINT (10 10)", "SRID=12;LINESTRING (0 0, 1 1)"),
      lengths = 3L
    )
  )

  expect_identical(
    wkb_unnest(
      as_wkb("GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)"),
      keep_multi = TRUE, max_depth = 2, keep_empty = FALSE
    ),
    structure(
      wkt_translate_wkb(c("MULTIPOINT ((30 10), (10 10))", "LINESTRING (0 0, 1 1)")),
      lengths = 2L
    )
  )

  expect_identical(
    wkb_unnest(
      as_wkb("GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)"),
      keep_multi = TRUE, max_depth = 2, keep_empty = TRUE
    ),
    structure(
      wkt_translate_wkb(c("MULTIPOINT ((30 10), (10 10))", "LINESTRING (0 0, 1 1)", "GEOMETRYCOLLECTION EMPTY")),
      lengths =  3L
    )
  )

  expect_identical(
    wkb_unnest(
      as_wkb("GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)"),
      keep_multi = FALSE, max_depth = 2, keep_empty = FALSE
    ),
    structure(
      wkt_translate_wkb(c("POINT (30 10)", "POINT (10 10)", "LINESTRING (0 0, 1 1)")),
      lengths = 3L
    )
  )

  expect_identical(
    wksxp_unnest(
      as_wksxp("GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)"),
      keep_multi = TRUE, max_depth = 2, keep_empty = TRUE
    ),
    structure(
      wkt_translate_wksxp(c("MULTIPOINT ((30 10), (10 10))", "LINESTRING (0 0, 1 1)", "GEOMETRYCOLLECTION EMPTY")),
      lengths = 3L
    )
  )

  expect_identical(
    wksxp_unnest(
      as_wksxp("GEOMETRYCOLLECTION(MULTIPOINT (30 10, 10 10), LINESTRING (0 0, 1 1), GEOMETRYCOLLECTION EMPTY)"),
      keep_multi = FALSE, max_depth = 2, keep_empty = FALSE
    ),
    structure(
      wkt_translate_wksxp(c("POINT (30 10)", "POINT (10 10)", "LINESTRING (0 0, 1 1)")),
      lengths = 3L
    )
  )
})

test_that("wk*_unnest(max_depth) is respected", {
  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      max_depth = 0
    ),
    structure(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      lengths = 1L
    )
  )

  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      max_depth = 1
    ),
    structure("GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1)))", lengths = 1L)
  )

  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      max_depth = 2
    ),
    structure("GEOMETRYCOLLECTION (POINT (0 1))", lengths = 1L)
  )

  expect_identical(
    wkt_unnest(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      max_depth = 3
    ),
    structure("POINT (0 1)", lengths = 1L)
  )

  expect_identical(
    wkt_unnest(
      "SRID=21;GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0 1))))",
      max_depth = 3
    ),
    structure("SRID=21;POINT (0 1)", lengths = 1L)
  )
})
