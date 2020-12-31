
test_that("person_fit - Item", {
  item <- generate_item(model = "3PL")
  expect_is(person_fit(ip = item, resp = 0, theta = 0), "numeric")
})

test_that("person_fit - Itempool", {
  ip <- generate_ip(model = "2PL")
  theta <- rnorm(1)
  resp <- sim_resp(ip = ip, theta = theta)
  expect_is(person_fit(ip = ip, resp = resp, theta = theta), "numeric")
})


test_that("person_fit - Testlet", {
  t1 <- testlet(generate_ip(n = 3))
  theta <- rnorm(1)
  resp <- sim_resp(ip = t1, theta = theta)
  expect_is(person_fit(ip = t1, resp = resp, theta = theta), "numeric")

})
