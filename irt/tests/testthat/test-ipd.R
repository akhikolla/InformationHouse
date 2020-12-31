
test_that("ipd", {
  n <- 20
  ip1 <- generate_ip(n = n)
  ipdf1 <- ipdf2 <- as.data.frame(ip1)
  ipdf2$a <- ipdf1$a + runif(n, -0.25, .25)
  ipdf2$b <- ipdf1$b + runif(n, -0.25, .25)
  ipdf2$a[1] <- ipdf1$a[1] + 1
  ipdf2$b[n] <- ipdf1$b[n] - 1
  ip2 <- itempool(ipdf2)
  result <- ipd(ip1 = ip1, ip2 = ip2)
  expect_true("Item-1" %in% result$a$unstable)
  expect_true(paste0("Item-", n) %in% result$b$unstable)
})
