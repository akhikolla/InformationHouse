# generate data
n = 100
coords = matrix(runif(n*2), nrow = n, ncol = 2)

# create covariance matrix
d = as.matrix(dist(coords))
V = 3*exp(-d/2) + 0.1*diag(n)

d1 = decomp_cov(V, "chol")
d2 = decomp_cov(V, "eigen")
d3 = decomp_cov(V, "svd")

range(V - tcrossprod(d1))
range(V - tcrossprod(d2))
range(V - tcrossprod(d3))

test_that("decomp_cov takes valid arguments",{
  # "V should be a matrix or Matrix"
  expect_that(decomp_cov(as.data.frame(V), "blah"), throws_error())
  # "V must be a square numeric matrix"
  expect_that(decomp_cov(V[-1,], "blah"), throws_error())
  # "method must be 'chol', 'eigen', or 'svd'"
  expect_that(decomp_cov(V, "blah"), throws_error())
})

test_that("decomp_cov is accurate for Matrix class matrices",{
  expect_true(max(abs(range(V - tcrossprod(d1)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d2)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d3)))) < 1e-10)
})

V = Matrix::Matrix(V)
d1 = decomp_cov(V, "chol")
d2 = decomp_cov(V, "eigen")
d3 = decomp_cov(V, "svd")

test_that("decomp_cov is accurate for Matrix matrices",{
  expect_true(max(abs(range(V - tcrossprod(d1)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d2)))) < 1e-10)
  expect_true(max(abs(range(V - tcrossprod(d3)))) < 1e-10)
})
