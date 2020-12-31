context("Testing the class cbkmCopula...")
library(cort)

dataset <- apply(LifeCycleSavings,2,rank)/(nrow(LifeCycleSavings)+1)
true_copula <- cbCopula(dataset[,2:3],m=5)

cop <- cbkmCopula(x = dataset,
                  m = 2,
                  pseudo = TRUE,
                  margins_numbers = c(2,3),
                  known_cop = true_copula)
cop2 <- cbkmCopula(x = dataset,
                  m = 2,
                  pseudo = TRUE,
                  margins_numbers = c(2,3),
                  known_cop = true_copula)


quiet(show(cop))

u=matrix(rep(0,15),ncol=5)
v=matrix(seq(0,1,length.out=15),ncol=5)
w=matrix(rep(1,15),ncol=5)

# testthat::test_that("zero-row or null data.frame are coerced to indepcopula", {
#   testthat::expect_error(cbkmCopula(as.data.frame(NULL),))
#   testthat::expect_equal(cbkmCopula(as.data.frame(matrix(0,nrow=0,ncol=5))),indepCopula(5))
# })

testthat::test_that("Inheritance and methods are there", {
  testthat::expect_s4_class(cop,"empiricalCopula")
  testthat::expect_s4_class(cop2,"empiricalCopula")
})

testthat::test_that("non-dividors m are not allowed", {
  testthat::expect_error(cbkmCopula(dataset,m=3))
  testthat::expect_error(cbkmCopula(dataset,m=7))
  testthat::expect_error(cbkmCopula(dataset,m=11))
  testthat::expect_error(cbkmCopula(dataset,m=49))
})

testthat::test_that("dimention of cbCopula is equal dimention of data", {
  testthat::expect_equal(dim(cbCopula(dataset)),ncol(dataset))
})

testthat::test_that("pCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop2) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),cop2) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),cop2),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),cop2))
})

testthat::test_that("vCopula did not change",{
  testthat::expect_equal(vCopula(u,v,cop)[1],0)
  testthat::expect_error(vCopula(v,u,cop))
  testthat::expect_equal(vCopula(u,v,cop2)[1],0)
  testthat::expect_error(vCopula(v,u,cop2))
  #testthat::expect_equal(vCopula(u,w,cop),rep(1,3))
  #testthat::expect_equal(vCopula(u,w,cop2),rep(1,3))
})

testthat::test_that("dim is ok",{
  testthat::expect_equal(dim(cop),5)
  testthat::expect_equal(dim(cop2),5)
})

testthat::test_that("rCopula output is ok for cbkmCopula",{
  testthat::expect_equivalent(rCopula(0,cop),matrix(ncol=5,nrow=0))
  testthat::expect_is(rCopula(10,cop),"matrix")
  testthat::expect_equal(ncol(rCopula(10,cop)),dim(cop))
  testthat::expect_equivalent(rCopula(0,cop2),matrix(ncol=5,nrow=0))
  testthat::expect_is(rCopula(10,cop2),"matrix")
  testthat::expect_equal(ncol(rCopula(10,cop2)),dim(cop2))
})


testthat::test_that("dCopula is not defined for cbkm copula",{
  testthat::expect_error(dCopula(rCopula(10,cop),cop))
})

















