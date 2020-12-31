context("check accuracy of flex.test poisson")
set.seed(16)
data(nydf)
data(nyw)
coords = cbind(nydf$longitude, nydf$latitude)
cases = floor(nydf$cases)
pop = nydf$population
outp = flex.test(coords = coords, cases = cases, pop = pop,
                  w = nyw, k = 10, nsim = 99, alpha = 1)

invisible(capture.output(outp_ <- flex_test(coords = coords, cases = cases, pop = pop,
                 w = nyw, k = 10, nsim = 99, alpha = 1)))

# results taken from flex_test_ny_poisson_10nn_cartesian
test_that("check accuracy for flex.test poisson", {
  expect_equal(sort(outp$clusters[[1]]$locids),
               c(85, 86, 88, 89, 90, 92, 93))
  expect_equal(round(outp$clusters[[1]]$max_dist, 6), 0.245662)
  expect_equal(outp$clusters[[1]]$cases, 39)
  expect_equal(round(outp$clusters[[1]]$exp, 4), 16.3981)
  expect_equal(round(outp$clusters[[1]]$smr, 5), 2.37832)
  expect_equal(round(outp$clusters[[1]]$test_statistic, 4), 11.6713)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[1]]$pvalue, 0.01)

  expect_equal(sort(outp$clusters[[2]]$locids),
               c(1, 2, 13, 15, 47, 49, 51))
  expect_equal(round(outp$clusters[[2]]$max_dist, 7),  0.0507453)
  expect_equal(outp$clusters[[2]]$cases, 31)
  expect_equal(round(outp$clusters[[2]]$exp, 4), 13.4462)
  expect_equal(round(outp$clusters[[2]]$smr, 5), 2.30548)
  expect_equal(round(outp$clusters[[2]]$test_statistic, 5), 8.62939)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[2]]$pvalue, 0.05)

  expect_equal(sort(outp$clusters[[9]]$locids),
               c(102, 103, 106))
  expect_equal(round(outp$clusters[[9]]$max_dist, 6), 0.194769)
  expect_equal(outp$clusters[[9]]$cases, 11)
  expect_equal(round(outp$clusters[[9]]$exp, 5), 4.76234)
  expect_equal(round(outp$clusters[[9]]$smr, 5), 2.30979)
  expect_equal(round(outp$clusters[[9]]$test_statistic, 5), 3.00674)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp$clusters[[9]]$pvalue, 1)
})

# results taken from flex_test_ny_poisson_10nn_cartesian
test_that("check accuracy for flex_test poisson", {
  expect_equal(sort(outp_$clusters[[1]]$locids),
               c(85, 86, 88, 89, 90, 92, 93))
  expect_equal(round(outp_$clusters[[1]]$max_dist, 6), 0.245662)
  expect_equal(outp_$clusters[[1]]$cases, 39)
  expect_equal(round(outp_$clusters[[1]]$exp, 4), 16.3981)
  expect_equal(round(outp_$clusters[[1]]$smr, 5), 2.37832)
  expect_equal(round(outp_$clusters[[1]]$test_statistic, 4), 11.6713)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp_$clusters[[1]]$pvalue, 0.02)

  expect_equal(sort(outp_$clusters[[2]]$locids),
               c(1, 2, 13, 15, 47, 49, 51))
  expect_equal(round(outp_$clusters[[2]]$max_dist, 7),  0.0507453)
  expect_equal(outp_$clusters[[2]]$cases, 31)
  expect_equal(round(outp_$clusters[[2]]$exp, 4), 13.4462)
  expect_equal(round(outp_$clusters[[2]]$smr, 5), 2.30548)
  expect_equal(round(outp_$clusters[[2]]$test_statistic, 5), 8.62939)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp_$clusters[[2]]$pvalue, 0.09)

  expect_equal(sort(outp_$clusters[[9]]$locids),
               c(102, 103, 106))
  expect_equal(round(outp_$clusters[[9]]$max_dist, 6), 0.194769)
  expect_equal(outp_$clusters[[9]]$cases, 11)
  expect_equal(round(outp_$clusters[[9]]$exp, 5), 4.76234)
  expect_equal(round(outp_$clusters[[9]]$smr, 5), 2.30979)
  expect_equal(round(outp_$clusters[[9]]$test_statistic, 5), 3.00674)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outp_$clusters[[9]]$pvalue, 1)
})

