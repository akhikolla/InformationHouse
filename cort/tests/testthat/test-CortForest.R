context("Testing the class CortForest...")

quiet(library(cort))
quiet(suppressWarnings(library(purrr)))
model = CortForest(LifeCycleSavings[,1:2],verbose_lvl = 0,n_trees=3)

quiet(show(model))

model2 = quiet(CortForest(LifeCycleSavings[,1:3],verbose_lvl = 10,number_max_dim = 2,n_trees=2))
u=matrix(rep(0,9),ncol=3)
v=matrix(seq(0,1,length.out=9),ncol=3)
w=matrix(rep(1,9),ncol=3)

testthat::test_that("initialisation check of Cort are ok", {
  testthat::expect_error(model = CortForest(LifeCycleSavings,pseudo_data=TRUE,verbose_lvl = 0))
})

testthat::test_that("pCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 6),nrow=2),model2) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 6),nrow=2),model2) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(0,nrow=2,ncol=3),model2),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),model2))
})

testthat::test_that("dCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(dCopula(matrix(seq(0.3,1,length.out = 6),nrow=2),model2) >= c(0,0)))
  testthat::expect_error(dCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),model2))
})

testthat::test_that("vCopula did not change",{
  testthat::expect_equal(vCopula(u,v,model2)[1],0)
  testthat::expect_error(vCopula(v,u,model2))
})

testthat::test_that("dim is ok",{
  testthat::expect_equal(dim(model),2)
  testthat::expect_equal(dim(model2),3)
})

testthat::test_that("rCopula output is ok for Cort",{
  testthat::expect_equivalent(rCopula(0,model2),matrix(ncol=5,nrow=0))
  testthat::expect_is(rCopula(10,model2),"matrix")
  testthat::expect_equal(ncol(rCopula(10,model2)),dim(model2))
})
