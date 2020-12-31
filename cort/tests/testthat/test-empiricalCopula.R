context("Testing the class empiricalCopula...")
library(cort)

# we could check that validity is OK.
data("LifeCycleSavings")
cop <- cbCopula(LifeCycleSavings)

quiet(show(cop))

cop_wrong <- cop
cop_wrong@data <- matrix("plop")
class(cop_wrong) <- "empiricalCopula"
cop_wrong2 <- cop
class(cop_wrong2) <- "empiricalCopula"
cop_wrong2@data <- cop_wrong2@data*50

cop_wrong3 <- cop
cop_wrong3@data <- matrix(0,nrow=0,ncol=5)

testthat::test_that("validity checks for the empiricalCopula class", {
  testthat::expect_true(validObject(as(cop,"empiricalCopula")))
  testthat::expect_error(supressWarnings(validObject(cop_wrong)))
  testthat::expect_error(validObject(cop_wrong2))
  testthat::expect_error(validObject(cop_wrong4))
})
