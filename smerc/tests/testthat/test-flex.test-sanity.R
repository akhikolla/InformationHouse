set.seed(1)
coords = runif(8)
context("sanity check flex.test")
test_that("sanity checks for flex.test arguments", {
  expect_error(flex.test(coords))
  coords = data.frame(x = runif(4), y = runif(4))
  cases = 1:3
  expect_error(flex.test(coords, cases = cases))
  cases = 1:4
  expect_error(flex.test(coords, cases = as.factor(cases)))
  pop = list(1:3)
  expect_error(flex.test(coords, cases = cases, pop = pop))
  pop = list(1:4)
  expect_error(flex.test(coords, cases = cases, pop = pop))
  pop = 1:4
  expect_error(flex.test(coords, cases = cases, pop = factor(pop)))
  ex = 1:3
  expect_error(flex.test(coords, cases = cases, pop = pop, ex = ex))
  ex = 1:4
  alpha = -1
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = 1.1
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = c(0.1, 0.3)
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha))
  alpha = 0.1
  nsim = 0
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim))
  nsim = c(10, 20)
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim))
  nsim = 10
  k = 0.5
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         k = k))
  k = c(1, 2)
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         k = k))
  k = 2
  longlat = 1:2
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         k = k, longlat = longlat))
  longlat = 1
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         k = k, longlat = longlat))
  longlat = FALSE
  parallel = 1:2
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                         k = k, longlat = longlat,
                         parallel = parallel))
  parallel = 1
  expect_error(flex.test(coords, cases = cases, pop = pop,
                         ex = ex, alpha = alpha, nsim = nsim,
                        k = k, longlat = longlat,
                        parallel = parallel))
  parallel = TRUE
  w = 1:4
  expect_error(flex.test(coords, cases = cases, pop = pop, w = w))
  w = matrix(1:4)
  expect_error(flex.test(coords, cases = cases, pop = pop, w = w))
  w = diag(3)
  expect_error(flex.test(coords, cases = cases, pop = pop, w = w))
  w = matrix(factor(diag(4)), nrow = 4)
  expect_error(flex.test(coords, cases = cases, pop = pop, w = w))
  w = diag(4)
  expect_error(flex.test(coords, cases = cases, pop = pop, w = w, k = 10))
})
