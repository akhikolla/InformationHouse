context("Testing the class Cort...")
library(cort)
model = Cort(LifeCycleSavings[,1:3],verbose_lvl = 0)
model2 = quiet(Cort(LifeCycleSavings,verbose_lvl = 10,number_max_dim = 2))
u=matrix(rep(0,15),ncol=5)
v=matrix(seq(0,1,length.out=15),ncol=5)
w=matrix(rep(1,15),ncol=5)

quiet(show(model))
quiet(plot(model))
quiet(pairs(model))

rho = biv_rho(model)
tau = biv_tau(model)

loss = loss(model)
ctr_infl = constraint_infl(model)
qp = quad_prod_with_data(model)

smaller_model = project_on_dims(model,c(1,2))
kendall = kendall_func(model,seq(0,1,length.out=10),M=10)

testthat::test_that("initialisation check of Cort are ok", {
  testthat::expect_error(Cort(LifeCycleSavings,pseudo_data=TRUE,verbose_lvl = 0))
})

testthat::test_that("the quandratic norm is coherent with the quad prod", {
  testthat::expect_equal(quad_prod(model,model),quad_norm(model))
})


testthat::test_that("pCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),model2) >= c(0,0)))
  testthat::expect_true(all(pCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),model2) <= c(1,1)))
  testthat::expect_equal(pCopula(matrix(seq(0,0,length.out = 10),nrow=2),model2),c(0,0))
  testthat::expect_error(pCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),model2))
})

testthat::test_that("dCopula values are between 0 and 1 with OK bounds.",{
  testthat::expect_true(all(dCopula(matrix(seq(0.3,1,length.out = 10),nrow=2),model2) >= c(0,0)))
  testthat::expect_error(dCopula(matrix(seq(0.3,1,length.out = 8),nrow=2),model2))
})

testthat::test_that("vCopula did not change",{
  testthat::expect_equal(vCopula(u,v,model2)[1],0)
  testthat::expect_error(vCopula(v,u,model2))
})

testthat::test_that("dim is ok",{
  testthat::expect_equal(dim(model2),5)
  testthat::expect_equal(dim(model2),5)
})

testthat::test_that("rCopula output is ok for Cort",{
  testthat::expect_equivalent(rCopula(0,model2),matrix(ncol=5,nrow=0))
  testthat::expect_is(rCopula(10,model2),"matrix")
  testthat::expect_equal(ncol(rCopula(10,model2)),dim(model2))
})

testthat::test_that("projection is ok",{
  testthat::expect_equal(dim(smaller_model),2)
})

testthat::test_that("tau works",{
  testthat::expect_equal(dim(tau), c(dim(model),dim(model)))
  testthat::expect_equal(diag(tau), rep(1,dim(model)))
  expect_true(all(tau <= 1) & all(tau >=-1))
  testthat::expect_equal(biv_tau(smaller_model)[2],tau[2])
})

testthat::test_that("rho works",{
  testthat::expect_equal(dim(rho), c(dim(model),dim(model)))
  testthat::expect_equal(diag(rho), rep(1,dim(model)))
  expect_true(all(rho <= 1) & all(rho >=-1))
  testthat::expect_equal(biv_rho(smaller_model)[2],rho[2])
})






