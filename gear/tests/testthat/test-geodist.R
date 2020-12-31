if (requireNamespace("sp", quietly = TRUE)) {
  coords = matrix(rnorm(60), ncol = 2)
  coords1 = rbind(coords, coords[1:5,])
  coords2 = rbind(coords, coords[6:16,])

  # Euclidean 1 and 2 coords
  de1 = geodist(coords1)
  de2 = geodist(coords1, coords2)
  se1 = sp::spDists(coords1)
  se2 = sp::spDists(coords1, coords2)

  # Great circle 1 and 2 coords
  dg1 = geodist(coords1, longlat = TRUE)
  ds1 = sp::spDists(coords1, longlat = TRUE)
  dg2 = geodist(coords1, coords2, longlat = TRUE)
  ds2 = sp::spDists(coords1, coords2, longlat = TRUE)

  test_that("geodist accuracy w/ sp::spDists", {
    expect_equal(de1, se1)
    expect_equal(de2, se2)
    expect_equal(dg1, ds1)
    expect_equal(dg2, ds2)
  })

} else {
  warning("geodist tests not performed because sp package not installed")
}
