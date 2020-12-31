test_that("encode_sf works", {

  # 3D point
  point3d <- sf::st_point(
    matrix(c(8.69821, 50.10228, 10), ncol = 3, byrow = TRUE), dim = "XYZ")

  # 2D linestring
  line2d <- sf::st_linestring(
    matrix(c(8.69821, 50.10228,
             8.69567, 50.10201,
             8.68752, 50.09878), ncol = 2, byrow = TRUE))

  # 3D polygon
  poly3d <- sf::st_polygon(list(
    matrix(c(8.69821, 50.10228, 10,
             8.69567, 50.10201, 20,
             8.69150, 50.10063, 30,
             8.69821, 50.10228, 10), ncol = 3, byrow = TRUE)), dim = "XYM")

  # Test encode_sf()
  expect_error(encode_sf(sf::st_multilinestring()),
               "Invalid geometry type 'MULTILINESTRING' of input, only 'POINT', 'LINESTRING' and 'POLYGON' is supported.")
  ## sfg (XY, XYZ and XYM)
  expect_type(encode_sf(point3d), "character")
  expect_type(encode_sf(line2d), "character")
  expect_type(encode_sf(poly3d), "character")
  ## sfc (XY and XYZ)
  point3d_sfc <- sf::st_as_sfc(list(point3d, point3d), crs = 4326)
  line2d_sfc <- sf::st_as_sfc(list(line2d, line2d), crs = 4326)
  poly3d_sfc <- sf::st_as_sfc(list(poly3d, poly3d), crs = 4326)
  expect_type(encode_sf(point3d_sfc), "character")
  expect_type(encode_sf(line2d_sfc), "character")
  expect_type(encode_sf(poly3d_sfc), "character")
  ## sf (XY and XYZ)
  expect_type(encode_sf(sf::st_as_sf(point3d_sfc)), "character")
  expect_type(encode_sf(sf::st_as_sf(line2d_sfc)), "character")
  expect_type(encode_sf(sf::st_as_sf(poly3d_sfc)), "character")

})
