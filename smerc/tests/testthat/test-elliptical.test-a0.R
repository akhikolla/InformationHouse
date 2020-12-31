context("check accuracy of elliptical.test a = 0")
set.seed(10)
data(nydf)
coords = nydf[, c("longitude", "latitude")]
pop = nydf$population
cases = floor(nydf$cases)
shape = c(1, 1.5, 2, 3, 4, 5)
nangle = c(1, 4, 6, 9, 12, 15)
ex = sum(cases) / sum(pop) * pop
ubpop = 0.1
cl = NULL
min.cases = 2
alpha = 1

nsim = 19
pbapply::pboptions(type = "none")
out0 = elliptic.test(coords, cases, pop, nsim = nsim,
                     alpha = 1, a = 0, ubpop = 0.1)
locids01 = c(46, 53, 45, 41, 44, 38, 52, 43, 54, 39, 50, 15, 40, 49, 16, 14, 1,
             48, 37, 13, 17, 11, 2, 12, 47)
locids05 = c(64, 65, 67, 62, 68)
test_that("check accuracy for elliptical.test a = 0", {
  expect_equal(locids01,
               out0$clusters[[1]]$locids)
  expect_equal(0.046,
               round(out0$clusters[[1]]$semiminor_axis, 3))
  expect_equal(0.18,
               round(out0$clusters[[1]]$semimajor_axis, 2))
  expect_equal(90, out0$clusters[[1]]$angle - 90)
  expect_equal(4, out0$clusters[[1]]$shape)
  expect_equal(102678, out0$clusters[[1]]$pop)
  expect_equal(96, out0$clusters[[1]]$cases)
  expect_equal(53.59, round(out0$clusters[[1]]$ex, 2))
  expect_equal(1.79, round(out0$clusters[[1]]$smr, 2))
  expect_equal(1.96, round(out0$clusters[[1]]$rr, 2))
  expect_equal(15.416469, round(out0$clusters[[1]]$test_statistic, 6))
  expect_equal(15.416469, round(out0$clusters[[1]]$loglikrat, 6))
  # satscan p-value 0.00046
  expect_equal(0.05, out0$clusters[[1]]$pvalue)

  expect_equal(locids05,
               out0$clusters[[5]]$locids)
  expect_equal(0.046,
               round(out0$clusters[[5]]$semiminor_axis, 3))
  expect_equal(0.18,
               round(out0$clusters[[5]]$semimajor_axis, 2))
  expect_equal(4, out0$clusters[[5]]$shape)
  expect_equal(27490, out0$clusters[[5]]$pop)
  expect_equal(28, out0$clusters[[5]]$cases)
  expect_equal(14.35, round(out0$clusters[[5]]$ex, 2))
  expect_equal(1.95, round(out0$clusters[[5]]$smr, 2))
  expect_equal(2.00, round(out0$clusters[[5]]$rr, 2))
  expect_equal(5.244377, round(out0$clusters[[5]]$test_statistic, 6))
  expect_equal(5.244377, round(out0$clusters[[5]]$loglikrat, 6))
  # satscan p-value 0.872
  expect_equal(0.9, out0$clusters[[5]]$pvalue)
})
