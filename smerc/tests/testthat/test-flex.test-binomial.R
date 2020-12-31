context("check accuracy of flex.test binomial")
set.seed(16)
data(nydf)
data(nyw)
coords = cbind(nydf$longitude, nydf$latitude)
cases = floor(nydf$cases)
pop = nydf$population
outb = flex.test(coords = coords, cases = cases, pop = pop,
                 w = nyw, k = 10, nsim = 99, alpha = 1,
                 type = "binomial")

# results taken from flex_test_ny_binomial_10nn_cartesian
test_that("check accuracy for flex.test binomial", {
  expect_equal(sort(outb$clusters[[1]]$locids),
               c(85, 86, 88, 89, 90, 92, 93))
  expect_equal(round(outb$clusters[[1]]$max_dist, 6), 0.245662)
  expect_equal(outb$clusters[[1]]$cases, 39)
  expect_equal(outb$clusters[[1]]$pop, 31420)
  expect_equal(round(outb$clusters[[1]]$smr, 5), 2.37832)
  expect_equal(round(outb$clusters[[1]]$test_statistic, 4), 11.6797)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[1]]$pvalue, 0.01)

  expect_equal(sort(outb$clusters[[2]]$locids),
               c(1, 2, 13, 15, 47, 49, 51))
  expect_equal(round(outb$clusters[[2]]$max_dist, 7),  0.0507453)
  expect_equal(outb$clusters[[2]]$cases, 31)
  expect_equal(outb$clusters[[2]]$pop, 25764)
  expect_equal(round(outb$clusters[[2]]$smr, 5), 2.30548)
  expect_equal(round(outb$clusters[[2]]$test_statistic, 5), 8.63552)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[2]]$pvalue, 0.05)

  expect_equal(sort(outb$clusters[[9]]$locids),
               c(102, 103, 106))
  expect_equal(round(outb$clusters[[9]]$max_dist, 6), 0.194769)
  expect_equal(outb$clusters[[9]]$cases, 11)
  expect_equal(outb$clusters[[9]]$pop, 9125)
  expect_equal(round(outb$clusters[[9]]$smr, 5), 2.30979)
  expect_equal(round(outb$clusters[[9]]$test_statistic, 5), 3.00889)
  # p-values are tough to test, make sure results don't change
  # in future versions since these were manually checked
  expect_equal(outb$clusters[[9]]$pvalue, 1)
})
