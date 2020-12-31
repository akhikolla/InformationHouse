dtf = data.frame(a = 1:2, b = 3:4)

# create geardf with matrix coords (note column names)
coords = matrix(rnorm(4), ncol = 2)
colnames(coords) = c("u", "v")
gdf1 = geardf(dtf, coords)

# create geardf with data.frame coords
coords = as.data.frame(coords)
gdf2 = geardf(dtf, coords)

# create geardf using coordnames
dtf2 = cbind(dtf, coords)
# vector form of coordnames
gdf3 = geardf(dtf2, coordnames = c("u", "v"))
# formula form of coordnames
gdf4 = geardf(dtf2, coordnames = ~ u + v)
# column index forum of coordnames
gdf5 = geardf(dtf2, coordnames = 3:4)

test_that("geardf behaves as expected", {
  # data must be a data.frame
  expect_error(geardf(matrix(1:4, ncol = 2), matrix(1:4, ncol = 2)))
  dtf = data.frame(a = 1:2, b = 1:2)
  # too many columns in coords
  expect_error(geardf(dtf, matrix(1:6, ncol = 3)))
  # not enough columns in coords
  expect_error(geardf(dtf, matrix(1:2, ncol = 1)))
  # no coordinate names in coords
  expect_error(geardf(dtf, matrix(1:6, nrow = 3)))
  # wrong number of rows for coords
  m = matrix(1:6, ncol = 2)
  colnames(m) = letters[1:2]
  expect_error(geardf(dtf, m))

  # wrong number of rows for coords (data.frame)
  m = data.frame(u = 1:3, v = 4:6)
  expect_error(geardf(dtf, m))

  # wrong number of rows for coords (data.frame)
  m = data.frame(u = 1:2, v = 4:5)
  # only coords or coordnames can be non-NULL
  expect_error(geardf(dtf, m, coordnames = c("u", "v")))

  # columns in dtf not available
  expect_error(geardf(dtf, coordnames = c("u", "z")))

  expect_equal(gdf1, gdf2)
  expect_equal(gdf1, gdf3)
  expect_equal(gdf1, gdf4)
  expect_equal(gdf1, gdf5)
})




