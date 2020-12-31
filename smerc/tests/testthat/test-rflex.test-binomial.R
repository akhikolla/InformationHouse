context("check accuracy of rflex.test binomial")
set.seed(15)
data(nydf)
data(nyw)
outb = rflex.test(coords = cbind(nydf$longitude, nydf$latitude),
                  cases = floor(nydf$cases),
                  pop = nydf$population,
                  w = nyw, k = 50,
                  nsim = 99, alpha = 1, longlat = FALSE,
                  alpha1 = 0.2,
                  type = "binomial")
# results taken from rflex_test_ny_binomial_50nn_cartesian
test_that("check accuracy for rflex.test binomial ", {
  expect_equal(sort(outb$clusters[[1]]$locids),
               c(1:2, 13, 15, 27, 35, 37:38, 43, 46:47, 49, 51:53))
  expect_equal(round(outb$clusters[[1]]$max_dist, 5), 0.22483)
  expect_equal(outb$clusters[[1]]$cases, 79)
  expect_equal(round(outb$clusters[[1]]$exp, 4), 36.1025)
  expect_equal(round(outb$clusters[[1]]$smr, 5), 2.18822)
  expect_equal(round(outb$clusters[[1]]$loglik, 4), 20.8156)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[1]]$pvalue, 0.01)

  expect_equal(sort(outb$clusters[[2]]$locids), c(89, 90))
  expect_equal(round(outb$clusters[[2]]$max_dist, 6), 0.125172)
  expect_equal(outb$clusters[[2]]$cases, 13)
  expect_equal(round(outb$clusters[[2]]$exp, 5), 3.46124)
  expect_equal(round(outb$clusters[[2]]$smr, 5), 3.75588)
  expect_equal(round(outb$clusters[[2]]$loglik, 5), 7.75475)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(round(outb$clusters[[2]]$pvalue, 1), 0.3)

  expect_equal(sort(outb$clusters[[12]]$locids),
               c(135, 146, 208, 210))
  expect_equal(round(outb$clusters[[12]]$max_dist, 6), 0.030934)
  expect_equal(outb$clusters[[12]]$cases, 12)
  expect_equal(round(outb$clusters[[12]]$exp, 5), 5.59582)
  expect_equal(round(outb$clusters[[12]]$smr, 5), 2.14446)
  expect_equal(round(outb$clusters[[12]]$loglik, 5), 2.79008)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[12]]$pvalue, 1)
})
