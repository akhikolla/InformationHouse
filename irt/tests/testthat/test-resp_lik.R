
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp_lik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# resp_lik (Item) @#################################
###############################################################################@

test_that("resp_lik - Item", {
  ##  Single theta
  ip <- item(a = rlnorm(1, 0, .3), b = rnorm(1), c = 0)
  theta <- rnorm(1)
  expect_equivalent(resp_lik(resp = 1, ip = ip, theta = theta),
                    prob(ip = ip, theta = theta))
  expect_equivalent(resp_lik(resp = 0, ip = ip, theta = theta),
                    1-prob(ip = ip, theta = theta))

  ##  GRM
  ip <- item(a = rlnorm(1, 0, .3), b = sort(rnorm(3)))
  theta <- rnorm(1)
  resp <- 1
  expect_equivalent(resp_lik(ip = ip, resp = resp, theta = theta),
                    prob(ip = ip, theta = theta)[resp+1])
  resp <- 2
  expect_equivalent(resp_lik(ip = ip, resp = resp, theta = theta),
                    prob(ip = ip, theta = theta)[resp+1])

  ##  GPCM
  ip <- item(a = rlnorm(1, 0, .3), b = sort(rnorm(3)), model = "GPCM")
  theta <- rnorm(1)
  resp <- 1
  expect_equivalent(resp_lik(ip = ip, resp = resp, theta = theta),
                    prob(ip = ip, theta = theta)[resp+1])
  resp <- 2
  expect_equivalent(resp_lik(ip = ip, resp = resp, theta = theta),
                    prob(ip = ip, theta = theta)[resp+1])

})

