
test_that("wkt_grob() works", {
  expect_is(wkt_grob(character(0)), "gTree")
  expect_is(wkt_grob("POINT EMPTY"), "gTree")

  grob_points <- wkt_grob(
    c("POINT (0.1 0.1)", "POINT (0.9 0.9)"),
    pch = c(1, 16), col = c("black", "red"),
    default.units = "npc"
  )
  expect_is(grob_points, "points")
  expect_equal(grob_points$pch, c(1, 16))
  expect_equal(grob_points$gp$col, c("black", "red"))

  # EMPTY items should be recycled along but not included in grob
  grob_points_with_empty <- wkt_grob(
    c("GEOMETRYCOLLECTION EMPTY", "POINT (0.1 0.1)", "POINT (0.9 0.9)"),
    pch = c(3, 1, 16), col = c("magenta", "black", "red"),
    default.units = "npc"
  )
  expect_is(grob_points_with_empty, "points")
  expect_equal(grob_points_with_empty$pch, c(1, 16))
  expect_equal(grob_points_with_empty$gp$col, c("black", "red"))

  # MULTIPOINT has to be handled slightly differently
  grob_multipoint <- wkt_grob(
    c("MULTIPOINT EMPTY", "GEOMETRYCOLLECTION EMPTY",
      "MULTIPOINT ((0.1 0.1), (0.2 0.1))", "POINT (0.5 0.6)"),
    col = c("magenta", "magenta", "black", "green"),
    default.units = "npc"
  )
  expect_is(grob_multipoint, "points")
  expect_identical(grob_multipoint$gp$col, c("black", "black", "green"))

  grob_lines <- wkt_grob(
    c("LINESTRING (0.1 0.1, 0.9 0.9)", "LINESTRING (0.1 0.9, 0.9 0.1)"),
    col = c("black", "red"),
    default.units = "npc"
  )
  expect_is(grob_lines, "polyline")
  expect_equal(grob_lines$gp$col, c("black", "red"))

  # multilines are handled slightly differently because polylineGrob has no
  # pathId argument
  grob_multiline <- wkt_grob(
    c(
      "MULTILINESTRING ((0.1 0.1, 0.9 0.9), (0.9 0.1, 0.1 0.9))",
      "MULTILINESTRING ((0.5 0.5, 0.5 1))"
    ),
    col = c("green", "blue"),
    default.units = "npc"
  )
  expect_is(grob_multiline, "polyline")
  expect_equal(grob_multiline$gp$col, c("green", "green", "blue"))

  grob_poly <- wkt_grob(
    c("POLYGON ((0.1 0.1, 0.9 0.1, 0.9 0.9, 0.1 0.9, 0.1 0.1))"),
    fill = "grey90",
    default.units = "npc"
  )
  expect_is(grob_poly, "pathgrob")
  expect_equal(grob_poly$gp$fill, "grey90")


  grob_mixed <- wkt_grob(
    c("POLYGON ((0.1 0.1, 0.9 0.1, 0.9 0.9, 0.1 0.9, 0.1 0.1))", "POINT (0.5 0.5)"),
    col = c("blue", "green"),
    default.units = "npc"
  )
  expect_is(grob_mixed, "gTree")
  expect_length(grob_mixed$children, 2)
  expect_is(grob_mixed$children[[1]], "pathgrob")
  expect_identical(grob_mixed$children[[1]]$gp$col, "blue")
  expect_is(grob_mixed$children[[2]], "points")
  expect_identical(grob_mixed$children[[2]]$gp$col, "green")

  grob_nested_point <- wkt_grob(
    c(
      "GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (GEOMETRYCOLLECTION (POINT (0.2 0.2))))",
      "LINESTRING EMPTY",
      "POINT (0.5 0.5)",
      "GEOMETRYCOLLECTION (POINT (0.7 0.7))"
    ),
    col = c("magenta", "black", "cyan", "red"),
    default.units = "npc"
  )
  expect_is(grob_nested_point, "points")
  expect_identical(grob_nested_point$gp$col, c("magenta", "cyan", "red"))

  grid::grid.newpage()
  grid::grid.draw(grob_points)
  grid::grid.draw(grob_points_with_empty)
  grid::grid.draw(grob_multipoint)
  grid::grid.draw(grob_lines)
  grid::grid.draw(grob_multiline)
  grid::grid.draw(grob_poly)
  grid::grid.draw(grob_mixed)
  grid::grid.draw(grob_nested_point)
})

test_that("wkb_grob() works", {
  grob_points <- wkb_grob(
    as_wkb(c("POINT (0.1 0.1)", "POINT (0.9 0.9)")),
    pch = c(1, 16), col = c("black", "red"),
    default.units = "npc"
  )
  expect_is(grob_points, "points")
  expect_equal(grob_points$pch, c(1, 16))
  expect_equal(grob_points$gp$col, c("black", "red"))
})

test_that("wksxp_grob() works", {
  grob_points <- wksxp_grob(
    as_wksxp(c("POINT (0.1 0.1)", "POINT (0.9 0.9)")),
    pch = c(1, 16), col = c("black", "red"),
    default.units = "npc"
  )
  expect_is(grob_points, "points")
  expect_equal(grob_points$pch, c(1, 16))
  expect_equal(grob_points$gp$col, c("black", "red"))
})
