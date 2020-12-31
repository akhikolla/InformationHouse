
test_that("dif", {
  ip <- generate_ip()
  group <- c(rep("M", sample(200:400, 1)), rep("F", sample(200:400, 1)))
  theta <- c(rnorm(sum(group == "M"), -1), rnorm(sum(group == "M"), 1))
  resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
  expect_is(dif(resp = resp, group = group, focal_name = "M",
                ip = NULL, type = "mh"), "data.frame")
})
