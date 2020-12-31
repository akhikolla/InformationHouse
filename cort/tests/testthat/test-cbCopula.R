context("Testing the class cbCopula...")
library(cort)
data(LifeCycleSavings)
cop <- cbCopula(LifeCycleSavings,m=10)
u=matrix(rep(0,15),ncol=5)
v=matrix(seq(0,1,length.out=15),ncol=5)
w=matrix(rep(1,15),ncol=5)

quiet(show(cop))

testthat::test_that("zero-row or null data.frame are coerced to indepcopula", {
  testthat::expect_error(cbCopula(as.data.frame(NULL)))
  testthat::expect_error(cbCopula(as.data.frame(matrix(0,nrow=0,ncol=5))))
})

testthat::test_that("data must be provided", {
  testthat::expect_error(cbCopula())
})

testthat::test_that("Inheritance and methods are there", {
  testthat::expect_s4_class(cop,"empiricalCopula")
})

testthat::test_that("non-dividors m are not allowed", {
  testthat::expect_error(cbCopula(LifeCycleSavings,m=3))
  testthat::expect_error(cbCopula(LifeCycleSavings,m=7))
  testthat::expect_error(cbCopula(LifeCycleSavings,m=11))
  testthat::expect_error(cbCopula(LifeCycleSavings,m=49))
})

testthat::test_that("dimention of cbCopula is equal dimention of data", {
  testthat::expect_equal(dim(cbCopula(LifeCycleSavings)),ncol(LifeCycleSavings))
})

testthat::test_that("pCopula values are OK",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop))

  testthat::expect_equal(round(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop),3),c(0,0.024))
  testthat::expect_equal(pCopula(matrix(seq(1,1,length.out = 10),nrow=2),cop),c(1,1))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
})

testthat::test_that("vCopula did not change",{
  testthat::expect_equal(vCopula(u,v,cop)[1],0)
  testthat::expect_error(vCopula(v,u,cop))
  testthat::expect_equal(vCopula(u,w,cop),rep(1,3))
})

testthat::test_that("dim is ok",{
  testthat::expect_equal(dim(cop),ncol(LifeCycleSavings))
})

testthat::test_that("rCopula output is ok",{
  testthat::expect_is(rCopula(10,cop),"matrix")
  testthat::expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

testthat::test_that("dCopula is ok",{
  testthat::expect_equal(dCopula(u,cop),rep(0,3))
  testthat::expect_error(dCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop))
  testthat::expect_equal(dCopula(c(0.65,0.35,0.65,0.85,0.45),cop),2000)
})


testthat::test_that("rCopula output is ok for cbCopula",{
  testthat::expect_equivalent(rCopula(0,cop),matrix(ncol=4,nrow=0))
  testthat::expect_is(rCopula(10,cop),"matrix")
  testthat::expect_equal(ncol(rCopula(10,cop)),dim(cop))
})

testthat::test_that("pCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop))
})






















