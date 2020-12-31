test_that("decode_sf works", {

  # Encoded lines
  n <- 5
  encodedXYZ <- "B1Voz5xJ67i1Bgkh9B"
  encodedXY <- "BFoz5xJ67i1B1B7PlU9yB"
  encodedXYM <- as.factor("BlXoz5xJ67i1Bgkh9B1B7Pgkh9BzIhagkh9BqK-pB_ni6D")

  # Decode
  pointXYZ <- decode_sf(rep(encodedXYZ, n))
  lineXY <- decode_sf(rep(encodedXY, n))
  polyXYM <- decode_sf(rep(encodedXYM, n))

  # Test decode_sf()
  expect_s3_class(pointXYZ, c("sf", "data.frame"), exact = TRUE)
  expect_s3_class(lineXY, c("sf", "data.frame"), exact = TRUE)
  expect_s3_class(polyXYM, c("sf", "data.frame"), exact = TRUE)
  expect_s3_class(decode_sf(c("BlXoz5xJ67i1Bgkh9B", "BlXoz5xJ67i1Bgkh9B1B7Pgkh9BzIhagkh9B", "BlXoz5xJ67i1Bgkh9B1B7Pgkh9BzIhagkh9BqK-pB_ni6D")), c("sf", "data.frame"), exact = TRUE)
  expect_true(all(sf::st_geometry_type(pointXYZ) == "POINT"))
  expect_true(all(sf::st_geometry_type(lineXY) == "LINESTRING"))
  expect_true(all(sf::st_geometry_type(polyXYM) == "POLYGON"))
  expect_equal(nrow(pointXYZ), n)
  expect_equal(nrow(lineXY), n)
  expect_equal(nrow(polyXYM), n)

})
