set.seed(2)
coords = runif(8)
context("sanity check scan.test")
test_that("sanity checks for scan.test arguments", {
  expect_that(scan.test(coords), throws_error())
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_error(scan.test(coords, cases = cases))
  cases = 1:4
  expect_error(scan.test(coords, cases = as.factor(cases)))
  pop = list(1:3)
  expect_error(scan.test(coords, cases = cases, pop = pop))
  pop = list(1:4)
  expect_error(scan.test(coords, cases = cases, pop = pop))
  pop = 1:4
  expect_error(scan.test(coords, cases = cases, pop = factor(pop)))
  ex = 1:3
  expect_error(scan.test(coords, cases = cases, pop = pop, ex = ex))
  ex = 1:4
  alpha = -1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = 1.1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = c(0.1, 0.3)
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = 0.1
  nsim = c(10, 20)
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim))
  nsim = 10
  ubpop = -0.1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop))
  ubpop = 1.1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop))
  ubpop = 1.1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop))
  ubpop = c(0.1, 0.2)
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop))
  ubpop = 0.5
  longlat = 1:2
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop, longlat = longlat))
  longlat = 1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop, longlat = longlat))
  longlat = FALSE
  parallel = 1:2
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop, longlat = longlat,
                         parallel = parallel))
  parallel = 1
  expect_error(scan.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         ubpop = ubpop, longlat = longlat,
                         parallel = parallel))
})
