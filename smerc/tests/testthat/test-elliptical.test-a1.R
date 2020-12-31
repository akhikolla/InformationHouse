context("check accuracy of elliptical.test a = 1")
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
out1 = elliptic.test(coords, cases, pop, nsim = nsim,
                     alpha = 1, a = 1, ubpop = 0.1)
locids1.1 = c(52, 50, 53, 38, 49, 15, 48, 39, 1, 37, 16, 44,
              14, 47, 40, 2, 13, 43, 51, 45, 17, 11, 12, 3, 46)
locids1.7 = c(106, 103, 102, 77, 230)

test_that("check accuracy of elliptical.test a = 1", {
  expect_equal(locids1.1, out1$clusters[[1]]$locids)
  expect_equal(0.058,
               round(out1$clusters[[1]]$semiminor_axis, 3))
  expect_equal(0.087,
               round(out1$clusters[[1]]$semimajor_axis, 3))
  expect_equal(90 + 90, out1$clusters[[1]]$angle)
  expect_equal(1.5, out1$clusters[[1]]$shape)
  expect_equal(99685, out1$clusters[[1]]$pop)
  expect_equal(93, out1$clusters[[1]]$cases)
  expect_equal(52.03, round(out1$clusters[[1]]$ex, 2))
  expect_equal(1.79, round(out1$clusters[[1]]$smr, 2))
  expect_equal(1.95, round(out1$clusters[[1]]$rr, 2))
  expect_equal(14.181797, round(out1$clusters[[1]]$test_stat, 6))
  expect_equal(14.772705, round(out1$clusters[[1]]$loglikrat, 6))
  # true p-value  0.00016
  expect_equal(0.05, out1$clusters[[1]]$pvalue)

  expect_equal(locids1.7, out1$clusters[[7]]$locids)
  expect_equal(0.087,
               round(out1$clusters[[7]]$semiminor_axis, 3))
  expect_equal(0.35,
               round(out1$clusters[[7]]$semimajor_axis, 2))
  expect_equal(90 + 45, out1$clusters[[7]]$angle)
  expect_equal(4.00, out1$clusters[[7]]$shape)
  expect_equal(15576, out1$clusters[[7]]$pop)
  expect_equal(19, out1$clusters[[7]]$cases)
  expect_equal(8.13, round(out1$clusters[[7]]$ex, 2))
  expect_equal(2.34, round(out1$clusters[[7]]$smr, 2))
  expect_equal(2.38, round(out1$clusters[[7]]$rr, 2))
  expect_equal(3.436309, round(out1$clusters[[7]]$test_stat, 6))
  expect_equal(5.369233, round(out1$clusters[[7]]$loglikrat, 6))
  # true p-value  0.990
  expect_equal(1, out1$clusters[[7]]$pvalue)
})
