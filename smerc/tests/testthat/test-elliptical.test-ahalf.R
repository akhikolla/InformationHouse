context("check accuracy of elliptical.test a = 0.5")
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
out0.5 = elliptic.test(coords, cases, pop, nsim = nsim,
                       alpha = 1, a = 0.5, ubpop = 0.1)
locids0.51 = c(52, 50, 53, 38, 49, 15, 48, 39, 1, 37, 16, 44,
               14, 47, 40, 2, 13, 43, 51, 45, 17, 11, 12, 3, 46)
locids0.52 = c(87, 88, 86, 89, 92, 85, 90)
locids0.53 = c(115, 116, 114, 111, 117, 113, 123, 120, 110,
               122, 112, 118, 121, 124, 220, 133, 131, 119,
               130, 125, 132, 219, 126, 127, 135)
locids0.54 = c(170, 171, 166, 167)

test_that("check accuracy of elliptical.test a = 0.5", {
  expect_equal(locids0.51,
               out0.5$clusters[[1]]$locids)
  expect_equal(0.058,
               round(out0.5$clusters[[1]]$semiminor_axis, 3))
  expect_equal(0.087,
               round(out0.5$clusters[[1]]$semimajor_axis, 3))
  expect_equal(90, out0.5$clusters[[1]]$angle - 90)
  expect_equal(1.5, out0.5$clusters[[1]]$shape)
  expect_equal(99685, out0.5$clusters[[1]]$pop)
  expect_equal(93, out0.5$clusters[[1]]$cases)
  expect_equal(52.03, round(out0.5$clusters[[1]]$ex, 2))
  expect_equal(1.79, round(out0.5$clusters[[1]]$smr, 2))
  expect_equal(1.95, round(out0.5$clusters[[1]]$rr, 2))
  expect_equal(14.772705, round(out0.5$clusters[[1]]$loglikrat, 6))
  expect_equal(14.474236, round(out0.5$clusters[[1]]$test_statistic, 6))
  # true p-value 0.00019
  expect_equal(0.05, out0.5$clusters[[1]]$pvalue)

  expect_equal(locids0.52,
               out0.5$clusters[[2]]$locids)
  expect_equal(locids0.53,
               out0.5$clusters[[3]]$locids)
  expect_equal(locids0.54,
               out0.5$clusters[[4]]$locids)
})
