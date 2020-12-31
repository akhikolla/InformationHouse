# load saved data
fpath = system.file("testdata",  "evgram_data.rda", package = "gear")
load(fpath)

test_that("check accuracy of evgram function", {
  # for omnidirectional standard semivarigoram
  gear_v1 = evgram(cadmium ~ 1, meuse, nbins = 10, verbose = FALSE)

  expect_equal(gear_v1$semi$np, gstat_v1$np)
  expect_equal(gear_v1$semi$dist, gstat_v1$dist)
  # expect_true(max(abs(gear_v1$semi$semivariance - gstat_v1$gamma)) < 1e-10)
  expect_equal(gear_v1$semi$semivariance, gstat_v1$gamma)

  gear_v1b = evgram(cadmium ~ 1, meuse_df,
                    coordnames = ~ x + y, nbins = 10,
                    verbose = FALSE)
  # make sure estimated semivariogram matches
  expect_equal(gear_v1$semi, gear_v1b$semi)

  # for omnidirectional cressie semivarigoram
  gear_v2 = evgram(cadmium ~ 1, meuse, nbins = 10, type = "cressie", verbose = FALSE)
  expect_equal(gear_v2$semi$np, gstat_v2$np)
  expect_equal(gear_v2$semi$dist, gstat_v2$dist)
  expect_equal(gear_v2$semi$semivariance, gstat_v2$gamma)

  # for directional standard semivariogram
  # very similar, but some negligible discrepancies based on how distances are calculated
  # angles difference because gstat does clockwise directions, using straight north as 0
  gear_v3 = evgram(cadmium ~ 1, meuse, nbins = 10,
                   angle = 22.5, ndir = 4, npmin = 1,
                   verbose = FALSE, invert = TRUE)
  v3i = which(gear_v3$semivariogram$np == gstat_v3$np)
  expect_equal(gear_v3$semi$dist[v3i], gstat_v3$dist[v3i])
  expect_equal(gear_v3$semi$semivariance[v3i], gstat_v3$gamma[v3i])

  gear_v4 = evgram(cadmium ~ 1, meuse, nbins = 10,
                   angle = 70, ndir = 4, npmin = 1,
                   verbose = FALSE, invert = TRUE)
  v4i = which(gear_v4$semivariogram$np == gstat_v4$np)
  expect_equal(gear_v4$semi$dist[v4i], gstat_v4$dist[v4i])
  expect_equal(gear_v4$semi$semivariance[v4i], gstat_v4$gamma[v4i])

  gear_v5 = evgram(cadmium ~ 1, meuse, nbins = 10,
                   angle = 0, ndir = 4, npmin = 1,
                   invert = TRUE, verbose = FALSE)
  v5i = which(gear_v5$semivariogram$np == gstat_v5$np)
  expect_equal(gear_v5$semi$dist[v5i], gstat_v5$dist[v5i])
  expect_equal(gear_v5$semi$semivariance[v5i], gstat_v5$gamma[v5i])

  gear_v6 = evgram(cadmium ~ 1, meuse, nbins = 10,
                   angle = 115, ndir = 4, npmin = 1,
                   invert = TRUE, verbose = FALSE)
  v6i = which(gear_v6$semivariogram$np == gstat_v6$np)
  expect_equal(gear_v6$semi$dist[v6i], gstat_v6$dist[v6i])
  expect_equal(gear_v6$semi$semivariance[v6i], gstat_v6$gamma[v6i])

  gear_v7 = evgram(cadmium ~ 1, meuse, nbins = 10,
                   angle = 160, ndir = 4, npmin = 1,
                   invert = TRUE, verbose = FALSE)
  v7i = which(gear_v7$semivariogram$np == gstat_v7$np)
  expect_equal(gear_v7$semi$dist[v7i], gstat_v7$dist[v7i])
  expect_equal(gear_v7$semi$semivariance[v7i], gstat_v7$gamma[v7i])

  # test with trend
  gear_v8 = evgram(cadmium ~ x + y, meuse, nbins = 10, verbose = FALSE)
  expect_equal(gear_v8$semi$np, gstat_v8$np)
  expect_equal(gear_v8$semi$dist, gstat_v8$dist)
  expect_equal(gear_v8$semi$semivariance, gstat_v8$gamma)

  # test cloud with and without maximum distance
  gear_cloud = evgram(cadmium ~ 1, meuse, nbins = 10,
                      maxd = 5000, type = "cloud",
                      verbose = FALSE)
  expect_equal(sort(gstat_cloud$dist), sort(gear_cloud$semivariogram$distance))
  expect_equal(sort(gstat_cloud$gamma), sort(gear_cloud$semivariogram$semivariance))
  # expect_equal(geo_cloud$u, gear_cloud$semi$distance)
  # expect_equal(geo_cloud$v, gear_cloud$semi$semivariance)

  gear_cloud_maxd2000 = evgram(cadmium ~ 1, meuse,
                               nbins = 10, type = "cloud",
                               verbose = FALSE, maxd = 2000)
  expect_equal(sort(gstat_cloud_maxd2000$dist), sort(gear_cloud_maxd2000$semivariogram$distance))
  expect_equal(sort(gstat_cloud_maxd2000$gamma), sort(gear_cloud_maxd2000$semivariogram$semivariance))
  # expect_equal(geo_cloud_maxd2000$u, gear_cloud_maxd2000$semi$distance)
  # expect_equal(geo_cloud_maxd2000$v, gear_cloud_maxd2000$semi$semivariance)
})
