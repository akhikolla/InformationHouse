context("Testing the method vCopula")
quiet(library(cort))

# define some copulas :
data(LifeCycleSavings)
cop1 <- cbCopula(LifeCycleSavings[,2:4],m=10)
cop2 <- Cort(LifeCycleSavings[1:20,2:4],verbose_lvl = 0)
v=matrix(seq(0,1,length.out=12),ncol=3)
u=matrix(rep(0,12),ncol=3)

testthat::test_that("return 0 for equals inputs", {
  testthat::expect_equal(vCopula(rep(0.5,3),rep(0.5,3),cop1),0)
  testthat::expect_equal(vCopula(u,u,cop1),rep(0,4))
  testthat::expect_equal(vCopula(v,v,cop1),rep(0,4))
  testthat::expect_equal(vCopula(rep(0.5,3),rep(0.5,3),cop2),0)
  testthat::expect_equal(vCopula(u,u,cop2),rep(0,4))
  testthat::expect_equal(vCopula(v,v,cop2),rep(0,4))
})

test_that("return 1 for complete inputs", {
  testthat::expect_equal(vCopula(rep(0,3),rep(1,3),cop1),1)
  testthat::expect_equal(vCopula(rep(0,3),rep(1,3),cop2),1)

})

testthat::test_that("returns error for non-agreing dimentional inputs", {
  testthat::expect_error(vCopula(1,c(1,2),cop1))
  testthat::expect_error(vCopula(c(1,2,3),c(1,2),cop1))
  testthat::expect_error(vCopula(u,1,cop1))
  testthat::expect_error(vCopula(1,c(1,2),cop2))
  testthat::expect_error(vCopula(c(1,2,3),c(1,2),cop2))
  testthat::expect_error(vCopula(u,1,cop2))
})

