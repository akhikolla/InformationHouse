
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp_loglik %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# resp_loglik (Item) @##############################
###############################################################################@

test_that("resp_loglik - Item", {
  ##  Single theta
  ip <- item(a = 1, b = 0.02173983, c = 0, D = 0.7636914)
  expect_equal(
    resp_loglik(resp = 1, ip = ip, theta = 0.4130181),
    -0.5548593,
    check.attributes = FALSE, tolerance = 1e-5)
  # 3PL
  ip <- item(a = 1.366882, b = -1.10991, c =  0.1462888,
                           D = 0.8379783)
  expect_equal(resp_loglik(resp = 1, ip = ip, theta = -0.2072357), -0.2535339,
               check.attributes = FALSE, tolerance = 1e-5)

  # -------------------------------------------------------------------------- #
  theta = -0.228
  a = 1.29
  b = -1.029
  D = 1.702
  ip = item(a = a, b = b, D = D, model = '2PL', id = paste0('i', 1:1))
  expected <- -1.91760642
  object <- resp_loglik(ip = ip, resp = 0, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  theta = -0.62
  a = 1.165
  b = 1.022
  D = 1.702
  ip = item(a = a, b = b, D = D, model = '2PL', id = paste0('i', 1:1))
  expected <- -3.29363208
  object <- resp_loglik(ip = ip, resp = 1, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # Multiple theta
  ip <- item(a = 1.428541, b = 0.6719602, c =  0.05310894,
                           d = 0.8918666, D = 1.570897)
  expect_equal(resp_loglik(resp = c(0, 0, 0, 0, 1, 1, 0), ip = ip,
                            theta = -3:3)[4], -0.2295754,
               check.attributes = FALSE, tolerance = 1e-5)
  expect_equal(resp_loglik(resp = tibble::tibble(item1 = c(0, 0, 0, 0, 1, 1, 0)),
                           ip = ip, theta = -3:3)[4], -0.2295754,
               check.attributes = FALSE, tolerance = 1e-5)

  ## Response Likelihood
  # It is simply exp(resp_loglik)
  # Single theta
  ip <- item(a = 1, b = 0.02173983, c = 0, D = 1.7)
  theta <- 0.4130181
  p <- prob(ip, theta)
  resp <- 1
  expect_equal(exp(resp_loglik(resp = resp, ip = ip, theta = theta)), p,
               check.attributes = FALSE, tolerance = 1e-5)
  # 3PL
  ip <- item(a = 1.366882, b = -1.10991, c =  0.1462888, D = 1)
  theta <- -.53
  resp <- 0
  p <- prob(ip, theta)
  expect_equal(exp(resp_loglik(resp = resp, ip = ip, theta = theta)), 1-p,
               check.attributes = FALSE, tolerance = 1e-5)

  # Multiple theta
  ip <- item(a = 1.4, b = 0.6, c =  0.05, d = 0.96, D = 1)
  resp = c(0, 0, 0, 0, 1, 1, 0)
  theta = -3:3
  p <- prob(ip, theta)
  expect_equal(exp(resp_loglik(resp = resp, ip = ip, theta = theta)),
               p^resp * (1-p)^(1-resp),
               check.attributes = FALSE, tolerance = 1e-5)
  # item pool
  ip <- itempool(data.frame(
    a = c(1.36, 1.19, 0.96, 1.41, 1.1, 1.44, 2.42),
    b = c(-0.22, 1.56, -1.34, -0.37, -0.16, 0.09, -0.75),
    c = c(0.13, 0.25, 0.25, 0.20, 0.19, 0.05, 0.05),
    D = 1.702))
  resp = c(0, 0, 1, 0, 0, 0, 1)
  theta = -0.55
  p <- prob(ip, theta)
  expect_equal(
    exp(resp_loglik(resp = resp, ip = ip, theta = theta)),
    prod(p^resp * (1-p)^(1-resp)), check.attributes = FALSE, tolerance = 1e-5)

  ### Missing Responses ###
  # -------------------------------------------------------------------------- #
  # Missing Responses should return NA
  ip <- item(a = 1.4, b = 0.6, c =  0.05, d = 0.96, D = 1)
  resp = c(NA, 0, 0, 0, 1, NA, 0)
  theta = -3:3
  expect_equal(which(is.na(resp_loglik(resp = resp, ip = ip, theta = theta))),
               c(1, 6))

  # -------------------------------------------------------------------------- #
  # If resp is all missing data, then resp_loglik should be NA.
  ip <- item(a = 1.4, b = 0.6, c =  0.05, d = 0.96, D = 1)
  theta = -3:3
  resp = rep(NA, length(theta))
  ll <- resp_loglik(ip = ip, resp = resp, theta = theta)
  expect_equal(length(ll), length(theta))
  expect_true(all(is.na(ll)))

  # -------------------------------------------------------------------------- #
  # Missing responses should not change the log likelihood or a response string
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .7, 1.5), b = runif(n, -2, 2))
  ip_list <- ip$item_list
  resp <- sim_resp(ip = ip, theta = rnorm(1, 0, .4))
  # Randomly drop two items and create an item pool
  selected_items <- sample(1:n, 2)
  ip_small <- itempool(ip_list[-selected_items])
  resp_small <- resp[, -selected_items]
  resp_na <- resp
  resp_na[, selected_items] <- NA
  theta <- rnorm(1)
  expect_equal(resp_loglik(ip = ip, resp = resp_na, theta = theta),
               resp_loglik(ip = ip_small, resp = resp_small, theta = theta))

  # -------------------------------------------------------------------------- #
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .7, 1.5), b = runif(n, -2, 2))
  ip_list <- ip$item_list
  resp <- tibble::as_tibble(sim_resp(ip = ip, theta = rnorm(1, 0, .4)))
  # Randomly drop two items and create an item pool
  selected_items <- sample(1:n, 2)
  ip_small <- itempool(ip_list[-selected_items])
  resp_small <- resp[, -selected_items]
  resp_na <- resp
  resp_na[, selected_items] <- NA
  theta <- rnorm(1)
  expect_equal(resp_loglik(ip = ip, resp = resp_na, theta = theta),
               resp_loglik(ip = ip_small, resp = resp_small, theta = theta))


  ### Graded Response Model ###
  # Item
  # -------------------------------------------------------------------------- #
  theta = -0.472
  D = 1.702
  a = 1.08
  b = c(-1.231, 0.805, 1.294)
  ip = item(a = a, b = b, D = D, model = 'GRM')
  expected <- -0.336680359157912
  object <- resp_loglik(ip = ip, resp = 1, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # Item Pool
  # -------------------------------------------------------------------------- #
  theta = 1.376
  D = 1.702
  a = 1.27
  b = c(-2.663, 1.025, 1.257)
  item1 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.111
  b = c(0.772, 0.797, 1.259)
  item2 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.128
  b = c(-1.106, -0.626, -0.225)
  item3 = item(a = a, b = b, D = D, model = 'GRM')
  ip = itempool(c(item1, item2, item3), id = paste0('i', 1:3))
  expected <- -3.60888518
  object <- resp_loglik(ip = ip, resp = c(2, 0, 3), theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)


  ### Generalized Partial Credit Model ###
  # Item
  # -------------------------------------------------------------------------- #
  theta = 1.279
  D = 1.702
  a = 1.119
  b = c(-0.39, -0.494, -0.315)
  ip = item(a = a, b = b, D = D, model = 'GPCM')
  expected <- -3.08437963
  object <- resp_loglik(ip = ip, resp = 2, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # Item Pool
  # -------------------------------------------------------------------------- #
  theta = -0.894
  D = 1.702
  a = 0.93
  b = c(-1.053, -1.833, -1.299)
  item1 = item(a = a, b = b, D = D, model = 'GPCM')
  a = 1.053
  b = c(0.438, 1.254, -0.433)
  item2 = item(a = a, b = b, D = D, model = 'GPCM')
  a = 0.97
  b = c(0.259, -0.12, 0.564)
  item3 = item(a = a, b = b, D = D, model = 'GPCM')
  ip = itempool(c(item1, item2, item3), id = paste0('i', 1:3))
  expected <- -7.05057445
  object <- resp_loglik(ip = ip, resp = c(2, 0, 3), theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)


  ### Mixture of Models ###
  # -------------------------------------------------------------------------- #
  # Here the item pool consist of two 2PL, three GRM and two GPCM items.
  # This test checks whether the log likelihood is calculated correctly.
  theta = 1.272
  a = c(0.983, 1.03)
  b = c(-0.581, 1.643)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b), D = D, model = '2PL',
                   id = paste0('irt', 1:2))
  expected <- -3.56417023
  resp <- c(0, 0)
  # GRM
  a = 1.027
  b = c(-2.439, -0.21, 0.693)
  item1 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.024
  b = c(0.241, 0.257, 0.311)
  item2 = item(a = a, b = b, D = D, model = 'GRM')
  a = 0.896
  b = c(-0.225, 0.458, 0.764)
  item3 = item(a = a, b = b, D = D, model = 'GRM')
  ip = c(ip, itempool(c(item1, item2, item3), id = paste0('grm', 1:3)))
  resp = c(resp, 3, 2, 1)
  expected <- sum(expected, -6.75156645)
  a = 0.927
  b = c(2.086, 0.107)
  item1 = item(a = a, b = b, D = D, model = 'GPCM')
  a = 1.201
  b = c(-0.349, 1.162)
  item2 = item(a = a, b = b, D = D, model = 'GPCM')
  ip = c(ip, itempool(c(item1, item2), id = paste0('gpcm', 1:2)))
  expected <- sum(expected, -1.70721466)
  resp = c(resp, 0, 2)
  object <- resp_loglik(ip = ip, resp = resp, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
})

###############################################################################@
############################# resp_loglik (Testlet) %###########################
###############################################################################@

test_that("resp_loglik - Testlet", {
  # Ability Estimation with the testlets
  t1 <- testlet(itempool(b = rnorm(2), id = c("t1-i1", "t1-i2")),
                   id = "t1")
  n_theta <- sample(4:9, 1)
  resp <- sim_resp(ip = t1, theta = rnorm(n_theta))
  theta <- rnorm(n_theta)
  expect_equal(
    resp_loglik(resp = resp, ip = t1, theta = theta),
    resp_loglik(resp = resp, ip = t1@item_list, theta = theta))

  # -------------------------------------------------------------------------- #
  t1 <- testlet(itempool(b = rnorm(1), id = c("t1-i1")), id = "t1")
  n_theta <- sample(4:9, 1)
  resp <- sim_resp(ip = t1, theta = rnorm(n_theta))
  theta <- rnorm(n_theta)
  expect_equal(
    resp_loglik(resp = resp, ip = t1, theta = theta),
    resp_loglik(resp = resp, ip = t1@item_list, theta = theta))

  # -------------------------------------------------------------------------- #
  t1 <- testlet(itempool(b = rnorm(1), id = c("t1-i1")), id = "t1")
  n_theta <- 1
  resp <- as.vector(sim_resp(ip = t1, theta = rnorm(n_theta)))
  theta <- rnorm(n_theta)
  resp_loglik(resp = resp, ip = t1, theta = theta)
  expect_equal(
    resp_loglik(resp = resp, ip = t1, theta = theta),
    resp_loglik(resp = resp, ip = t1@item_list, theta = theta))

})

###############################################################################@
############################# resp_loglik (Itempool) %#########################
###############################################################################@

test_that("resp_loglik - Itempool", {

  ip <- itempool(data.frame(
    a = c(1.36947439378127, 1.19103434309363, 0.966255453648046,
          1.41118730790913, 1.03645211085677, 1.44584470195696,
          1.42817918118089),
    b = c(-0.221226601859114, 1.56170696559919, -1.34820736712213,
          -0.374769023260681, -0.164801504563551, 0.0918581153824936,
          -0.758465659387833),
    c = c(0.132225778349675, 0.252476210868917, 0.251043332414702,
          0.20630893916823, 0.190127589460462, 0.051204833923839,
          0.0552326969336718),
    D = 2.129246))
  expect_equal(
    resp_loglik(resp = c(0, 0, 1, 0, 0, 0, 1), ip = ip, theta = -0.5545929),
    -2.72455, check.attributes = FALSE, tolerance = 1e-5)

  # -------------------------------------------------------------------------- #
  # resp (numeric), Multiple theta, 1 item (Itempool)
  resp <- c(1,0)
  theta <- c(-1, 1, 2)
  expect_error(resp_loglik(resp = resp, ip = ip, theta = theta),
               regexp = "Invalid arguments.")

  # -------------------------------------------------------------------------- #
  # The length of the ip and the length of the resp should correspond
  ip <- itempool(b = rnorm(10))
  resp <- sample(0:1, 7, replace = TRUE)
  expect_error(resp_loglik(resp =resp, ip = ip, theta = rnorm(1)),
               regexp = "Invalid arguments.")
  resp <- tibble::as_tibble(matrix(resp, nrow = 1,
                           dimnames = list(c(), paste0("i", 1:7))))
  expect_error(resp_loglik(resp = resp, ip = ip, theta = rnorm(1)),
               regexp = "Invalid arguments.")

  # -------------------------------------------------------------------------- #
  # If resp is all missing data, then resp_loglik should be NA.
  expect_true(is.na(resp_loglik(ip = ip, resp = rep(NA, length(ip)),
                                theta = rnorm(1))))
  # If there is one item and multiple thetas and all resp are missing
  # it should return a vector of NAs
  ip <- ip[[1]]
  theta <- c(-1, 1, 2)
  ll <- resp_loglik(ip = ip, resp = rep(NA, length(theta)), theta = theta)
  expect_equal(length(ll), 3)
  expect_true(all(is.na(ll)))

})


###############################################################################@
############################# resp_loglik (REST) @##############################
###############################################################################@
test_that("resp_loglik - REST", {
  # Try character, 1 theta
  expect_error(resp_loglik(resp = 1, ip = "1", theta = 2),
               regexp = "Cannot convert object to an")
  expect_error(resp_loglik(resp = 3, ip = 1, theta = 2),
               regexp = "Insufficient parameters. Please give more information")
})



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp_loglik (First Derivative) %%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# resp_loglik (First Derivative) (Item) @###########
###############################################################################@

test_that("resp_loglik (First Derivative) - Item", {
  ##  Single theta
  ip <- item(a = 1.288916, b = 0.2149486, D = 2.393605)
  expect_equal(resp_loglik(resp = 0, ip = ip, theta = -0.7466429,
                           derivative = 1),
               -0.1510338, check.attributes = FALSE, tolerance = 1e-5)

  ##  Multiple theta
  # Multiple theta
  ip <- item(a = 0.9580645, b = -0.4588905, c = 0.08207303,
                           d = 0.8559631, D = 1.011083)
  expect_equal(resp_loglik(resp = c(0, 0, 0, 0, 0, 1, 1), ip = ip, theta = -3:3,
                           derivative = 1)[4], -0.3997864,
               check.attributes = FALSE, tolerance = 1e-5)

  ### Missing Responses ###
  # -------------------------------------------------------------------------- #
  # Missing Responses should return NA
  ip <- item(a = 1.4, b = 0.6, c =  0.05, d = 0.96, D = 1)
  resp = c(NA, 0, 0, 0, 1, NA, 0)
  theta = -3:3
  expect_equal(which(is.na(resp_loglik(resp = resp, ip = ip, theta = theta,
                                       derivative = 1))), c(1, 6))

  # -------------------------------------------------------------------------- #
  # Missing responses should not change the log likelihood or a response string
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .7, 1.5), b = runif(n, -2, 2))
  ip_list <- ip$item_list
  resp <- sim_resp(ip = ip, theta = rnorm(1, 0, .4))
  # Randomly drop two items and create an item pool
  selected_items <- sample(1:n, 2)
  ip_small <- itempool(ip_list[-selected_items])
  resp_small <- resp[, -selected_items]
  resp_na <- resp
  resp_na[, selected_items] <- NA
  theta <- rnorm(1)
  expect_equal(resp_loglik(ip = ip, resp = resp_na, theta = theta,
                           derivative = 1),
               resp_loglik(ip = ip_small, resp = resp_small, theta = theta,
                           derivative = 1))

  # -------------------------------------------------------------------------- #
  # All of the elements of response matrix should be numeric
  n <- sample(10:20, 1) # number of examinees
  ip <- item(b = rnorm(1))
  resp <- data.frame(i1 = sim_resp(ip = ip, theta = rnorm(n)))
  resp[1, 1] <- 'a'
  expect_warning(rllfd <- resp_loglik(ip = ip, resp = resp, theta = rnorm(n),
                                      derivative = 1))
  expect_true(is.na(rllfd[1]))

})


###############################################################################@
############################# resp_loglik (First Derivative) (Itempool) @######
###############################################################################@

test_that("resp_loglik (First Derivative) - Itempool", {
  # -------------------------------------------------------------------------- #
  # All of the elements of response matrix should be numeric
  n <- sample(10:20, 1) # number of examinees
  ip <- itempool(b = rnorm(sample(5:10, 1)))
  theta <- rnorm(n)
  resp1 <- resp2 <- data.frame(sim_resp(ip = ip, theta = rnorm(n)))
  resp2[1, 1] <- 'a'
  fd1 <- resp_loglik(ip = ip, resp = resp1, theta = theta, derivative = 1)
  expect_warning(fd2 <- resp_loglik(ip = ip, resp = resp2, theta = theta,
                                    derivative = 1))
  expect_equivalent(fd1[-1], fd2[-1])
  expect_equal(fd2[1], resp_loglik(ip = ip[-1], resp = resp1[, -1],
                                   theta = theta, derivative = 1)[1])
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp_loglik (Second Derivative) %%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# resp_loglik (Second Derivative) (Item) @##########
###############################################################################@


test_that("resp_loglik (Second Derivative) - Item", {
  ##  Single theta
  # 3PL
  ip <- item(a = 0.7356062, b = -0.5102405, c = 0.2897516,
                           D = 2.057638)
  expect_equal(resp_loglik(resp = 1, ip = ip, theta = -0.5288874,
                           derivative = 2),
               -0.1673824, check.attributes = FALSE, tolerance = 1e-5)

  # -------------------------------------------------------------------------- #
  # All of the elements of response matrix should be numeric, otherwise NAs
  # will be introduced for non-numeric response values.
  n <- sample(10:20, 1) # number of examinees
  ip <- item(b = rnorm(1))
  resp <- data.frame(i1 = sim_resp(ip = ip, theta = rnorm(n)))
  resp[1, 1] <- 'a'
  expect_warning(rllsd <- resp_loglik(ip = ip, resp = resp, theta = rnorm(n),
                                      derivative = 2))
  expect_true(is.na(rllsd[1]))

  ### Missing Responses ###
  # -------------------------------------------------------------------------- #
  # Missing Responses should return NA
  ip <- item(a = 1.4, b = 0.6, c =  0.05, d = 0.96, D = 1)
  resp = c(NA, 0, 0, 0, 1, NA, 0)
  theta = -3:3
  expect_equal(which(is.na(resp_loglik(resp = resp, ip = ip, theta = theta,
                                       derivative = 2))), c(1, 6))
})


###############################################################################@
############################# resp_loglik (Second Derivative) (Itempool) @#####
###############################################################################@

test_that("resp_loglik (Second Derivative) - Itempool", {
  # -------------------------------------------------------------------------- #
  # All of the elements of response matrix should be numeric
  n <- sample(10:20, 1) # number of examinees
  ip <- itempool(b = rnorm(sample(5:10, 1)))
  theta <- rnorm(n)
  resp1 <- resp2 <- data.frame(sim_resp(ip = ip, theta = rnorm(n)))
  resp2[1, 1] <- 'a'
  sd1 <- resp_loglik(ip = ip, resp = resp1, theta = theta, derivative = 2)
  expect_warning(sd2 <- resp_loglik(ip = ip, resp = resp2, theta = theta,
                                    derivative = 2))
  expect_equivalent(sd1[-1], sd2[-1])
  expect_equal(sd2[1],
    resp_loglik(ip = ip[-1], resp = resp1[, -1], theta = theta,
                derivative = 2)[1])

  # -------------------------------------------------------------------------- #
  n_examinee <- sample(20:30, 1)
  D = 1.702
  ip = itempool(data.frame(a = c(0.983, 1.03), b = c(-0.581, 1.643)),
                   D = D, model = '2PL', id = paste0('irt', 1:2))
  item1 = item(a = 1.027, b = c(-2.439, -0.21, 0.693), D = D, model = 'GPCM')
  item2 = item(a = 1.024, b = c(0.241, 0.257, 0.311), D = D, model = 'GPCM')
  item3 = item(a = 0.896, b = c(-0.225, 0.458, 0.764), D = D, model = 'GPCM')
  ip = c(ip, itempool(c(item1, item2, item3), id = paste0('grm', 1:3)))
  item1 = item(a = 0.927, b = c(2.086, 0.107), D = D, model = 'GPCM')
  item2 = item(a = 1.201, b = c(-0.349, 1.162), D = D, model = 'GPCM')
  ip = c(ip, itempool(c(item1, item2), id = paste0('gpcm', 1:2)))
  theta <- rnorm(n_examinee)
  resp <- sim_resp(ip = ip, theta = theta)
  # Check regular response matrix
  output <- resp_loglik(ip = ip, resp = resp, theta = theta, derivative = 2)
  expect_equal(length(output), n_examinee)
  # Check with NAs
  resp[sample(1:(prod(dim(resp))), nrow(resp) * 2)] <- NA
  output <- resp_loglik(ip = ip, resp = tibble::as_tibble(resp), theta = theta,
                        derivative = 2)
  expect_equal(length(output), n_examinee)
  # Check tibble
  output <- resp_loglik(ip = ip, resp = tibble::as_tibble(resp), theta = theta,
                        derivative = 2)
  expect_equal(length(output), n_examinee)



  # -------------------------------------------------------------------------- #
  # Missing responses should not change the log likelihood or a response string
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .7, 1.5), b = runif(n, -2, 2))
  ip_list <- ip$item_list
  resp <- sim_resp(ip = ip, theta = rnorm(1, 0, .4))
  # Randomly drop two items and create an item pool
  selected_items <- sample(1:n, 2)
  ip_small <- itempool(ip_list[-selected_items])
  resp_small <- resp[, -selected_items]
  resp_na <- resp
  resp_na[, selected_items] <- NA
  theta <- rnorm(1)
  expect_equal(resp_loglik(ip = ip, resp = resp_na, theta = theta,
                           derivative = 2),
               resp_loglik(ip = ip_small, resp = resp_small, theta = theta,
                           derivative = 2))

})


