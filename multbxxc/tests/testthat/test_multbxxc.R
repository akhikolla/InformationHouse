context("bop")
n=3
a=matrix(1, n, n)
b=matrix(1:n, n, 1)

a0=a+0
bop(a, 2, "+=", b)
test_that("bop +=", {
  expect_equal(a, a0+rep(b, n))
})
