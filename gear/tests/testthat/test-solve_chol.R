set.seed(10)
# create positive definite matrix a
a = crossprod(matrix(rnorm(25^2), nrow = 25))
# create vector x and matrix b
# x can be used to check the stability of the solution
x = matrix(rnorm(25))
b = a %*% x

# standard solve
x1 = solve(a, b)

# solve using cholesky decomposition
x2 = solve_chol(chol(a), b)

# solve using qr decomposition
x3 = solve.qr(qr(a), b)

# compare direct inversion
ai1 = solve(a)
ai2 = solve_chol(chol(a)) #using cholesky decomposition

test_that("solve_chol is correct", {
  expect_equal(x, x1)
  expect_equal(x, x2)
  expect_equal(x, x3)
  expect_equal(ai1, ai2)
})
