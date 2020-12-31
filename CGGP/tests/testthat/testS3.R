context("testS3")

test_that("print and predict work", {
  
  # Check print
  SG = CGGPcreate(3,201, corr="cauchy")
  tp <- capture.output(print(SG))
  expect_is(tp, "character")
  expect_gt(length(tp), 5)
  
  # Check predict
  f <- function(x) {x[1]^1.2*sin(x[2])}
  SG <- CGGPfit(SG, apply(SG$design, 1, f))
  xp <- matrix(runif(10*3), ncol=3)
  pp <- predict(SG, xp)
  expect_is(pp, "list")
  expect_equal(names(pp), c("mean", "var"))
  expect_equal(length(pp$mean), 10)
  expect_equal(length(pp$var), 10)
  
  # Check plot
  pl <- plot(SG)
  expect_is(pl, "ggplot")
})