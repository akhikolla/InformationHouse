N_TEST <- 5 # replications of tests with random input data

test_that("Set of detected change points increases with decreasing epsilon", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    ts <- list(testData(model="blocks"),
               testData(model="fms"),
               testData(model="mix"),
               testData(model="stairs10"),
               testData(model="teeth10"))
    for (x in ts) {
      cpts.prev <- numeric(0)
      cpts.prev.b <- numeric(0) # Include boundary extension
      G <- floor(runif(1, 5, 15))
      for (epsilon in (11:0)/10) { #TODO include 0
        cpts.curr <- mosum(x, G, boundary.extension=FALSE,
                                alpha=alpha, criterion="epsilon", epsilon=epsilon)$cpts
        cpts.curr.b <- mosum(x, G, boundary.extension=TRUE,
                                  alpha=alpha, criterion="epsilon", epsilon=epsilon)$cpts
        expect_true(all(cpts.prev %in% cpts.curr))
        expect_true(all(cpts.prev.b %in% cpts.curr.b))
        cpts.prev <- cpts.curr
        cpts.prev.b <- cpts.curr.b
      }
    }
  }
})

test_that("Set of detected change points increases with decreasing eta", {
  for (i in 1:N_TEST) {
    alpha <- runif(1, 0, 1)
    ts <- list(testData(model="blocks"),
               testData(model="fms"),
               testData(model="mix"),
               testData(model="stairs10"),
               testData(model="teeth10"))
    for (x in ts) {
      cpts.prev <- numeric(0)
      cpts.prev.b <- numeric(0) # Include boundary extension
      G <- floor(runif(1, 5, 15))
      for (eta in (20:0)/10) {
        cpts.curr <- mosum(x, G, boundary.extension=FALSE,
                                alpha=alpha, criterion="eta", eta=eta)$cpts
        cpts.curr.b <- mosum(x, G, boundary.extension=TRUE,
                                alpha=alpha, criterion="eta", eta=eta)$cpts
        expect_true(all(cpts.prev %in% cpts.curr))
        expect_true(all(cpts.prev.b %in% cpts.curr.b))
        cpts.prev <- cpts.curr
        cpts.prev.b <- cpts.curr.b
      }
    }
  }
})