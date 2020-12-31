context("Testing the class ConvexCombCopula...")
library(cort)

## ----fig.cap="constructing a copula"-----------------------
dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
copulas <- list(
  cbCopula(dataset[,2:3],m=10),
  cbCopula(dataset[,2:3],m=5)
)
alpha <- c(1,4)

cop <- ConvexCombCopula(copulas,alpha)

quiet(show(cop))

testthat::test_that("inputed argument for copulas is handeld in the right way.", {
  testthat::expect_warning(ConvexCombCopula(copulas = "abc"))
  testthat::expect_equal(suppressWarnings(ConvexCombCopula(copulas[[1]])),copulas[[1]])
  testthat::expect_equal(suppressWarnings(ConvexCombCopula(copulas[1])),copulas[[1]])
  testthat::expect_warning(ConvexCombCopula(copulas[[1]]))
  testthat::expect_warning(ConvexCombCopula(copulas[1]))
})

testthat::test_that("ConvexCombCopula class herited properly", {
  testthat::expect_is(cop,"ConvexCombCopula")
})

testthat::test_that("method dim for ConvexCombCopula class works properly", {
  testthat::expect_equal(dim(cop),dim(copulas[[1]]))
})

testthat::test_that("rCopula output is ok for COnvexCombCopula",{
  testthat::expect_equivalent(rCopula(0,cop),matrix(ncol=dim(cop),nrow=0))
  testthat::expect_is(rCopula(10,cop),"matrix")
  testthat::expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

testthat::test_that("pCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 4),nrow=2),cop) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 4),nrow=2),cop) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 4),nrow=2),cop),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop))
})
