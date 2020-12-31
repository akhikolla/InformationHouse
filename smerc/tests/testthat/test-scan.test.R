context("check scan.test accuracy")
set.seed(2)
coords = runif(8)
coords = data.frame(x = runif(4), y = runif(4))
data(nydf)
## update test to be faster by restricting population
out = scan.test(coords = cbind(nydf$longitude, nydf$latitude),
                cases = floor(nydf$cases), pop = nydf$population,
                longlat = FALSE, nsim = 49, alpha = 1,
                ubpop = 0.1)
cl1 = c(49, 48, 50, 47, 1, 15, 51, 2, 13, 14, 3, 35, 16,
        52, 36, 37, 12, 11, 5, 17, 10, 4, 6, 38, 18, 9)
cl2 = c(89, 88, 87, 84, 90, 86, 92, 85)
cl7 = c(166, 159, 167)

test_that("check accuracy for scan.test with SatScan for NY data", {
  expect_equal(out$clusters[[1]]$locids, cl1)
  expect_equal(round(out$clusters[[1]]$r, 3), 0.070)
  expect_equal(out$clusters[[1]]$pop, 95249)
  expect_equal(out$clusters[[1]]$cases, 85)
  expect_equal(round(out$clusters[[1]]$exp, 2), 49.71)
  expect_equal(round(out$clusters[[1]]$smr, 2), 1.71)
  expect_equal(round(out$clusters[[1]]$rr, 2), 1.84)
  expect_equal(round(out$clusters[[1]]$loglik, 6), 11.577255)
  expect_equal(round(out$clusters[[1]]$test_statistic, 6), 11.577255)
  # p-value from satscan 0.001
  expect_equal(out$clusters[[1]]$pvalue, 0.02)

  expect_equal(out$clusters[[2]]$locids, cl2)
  expect_equal(round(out$clusters[[2]]$r, 2), 0.13)
  expect_equal(out$clusters[[2]]$pop, 34083)
  expect_equal(out$clusters[[2]]$cases, 38)
  expect_equal(round(out$clusters[[2]]$exp, 2), 17.79)
  expect_equal(round(out$clusters[[2]]$smr, 2), 2.14)
  expect_equal(round(out$clusters[[2]]$rr, 2), 2.22)
  expect_equal(round(out$clusters[[2]]$loglik, 6), 9.019716)
  expect_equal(round(out$clusters[[2]]$test_statistic, 6), 9.019716)
  # p-value from satscan 0.013
  expect_equal(out$clusters[[2]]$pvalue, 0.02)

  expect_equal(out$clusters[[7]]$locids, cl7)
  expect_equal(round(out$clusters[[7]]$r, 3), 0.012)
  expect_equal(out$clusters[[7]]$pop, 8839)
  expect_equal(out$clusters[[7]]$cases, 11)
  expect_equal(round(out$clusters[[7]]$exp, 2), 4.61)
  expect_equal(round(out$clusters[[7]]$smr, 2), 2.38)
  expect_equal(round(out$clusters[[7]]$rr, 2), 2.41)
  expect_equal(round(out$clusters[[7]]$loglik, 6), 3.209485)
  expect_equal(round(out$clusters[[7]]$test_statistic, 6), 3.209485)
  # p-value from satscan 0.976
  expect_equal(out$clusters[[7]]$pvalue, 0.98)
})
