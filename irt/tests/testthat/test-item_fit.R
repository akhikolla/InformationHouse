
test_that("item_fit - Itempool", {
  ip <- generate_ip(model = "2PL")
  theta <- rnorm(100)
  resp <- sim_resp(ip = ip, theta = theta)

  expect_is(item_fit(ip = ip, resp = resp, theta = theta), "matrix")
})
