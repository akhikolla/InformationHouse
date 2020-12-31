
require(tibble)
# library(testthat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% est_ability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# The following test set is for overall functioning of est_ability.
test_that("est_ability", {
  # Check whether the function can also accept a 'tibble' as resp.
  ip <- itempool(b = rnorm(2))
  resp <- as_tibble(sim_resp(ip = ip, theta = rnorm(10)))
  expect_type(est_ability(ip = ip, resp = resp, method = "ml")$est, "double")

  # -------------------------------------------------------------------------- #
  # Ability estimate output will have the names of the response matrix.
  ip <- itempool(b = rnorm(2))
  theta_names <- paste0("Examinee-", 1:3)
  theta <- rnorm(3)
  names(theta) <- theta_names
  resp <- sim_resp(ip = ip, theta = theta)
  est <- est_ability(ip = ip, resp = resp, method = "ml")
  expect_equal(names(est$est), theta_names)
  expect_equal(names(est$se), theta_names)

  # -------------------------------------------------------------------------- #
  # Check whether both outputs, se and est, are vectors.
  n_examinee <- sample(10:20, 1)
  ip <- itempool(a = rlnorm(10, 0, .3), b = rnorm(10))
  resp <- sim_resp(ip = ip, theta = rnorm(n_examinee))
  est <- est_ability(ip = ip, resp = resp, method = "ml")
  expect_is(est$est, "numeric")
  expect_is(est$se, "numeric")
  expect_true(is.null(dim(est$est)))
  expect_true(is.null(dim(est$se)))
  expect_equal(length(est$est), n_examinee)
  expect_equal(length(est$se), n_examinee)


  # -------------------------------------------------------------------------- #
  # Standard errors should be different for the same response string but
  # different NA values.
  n_items <- 10 # Number of items
  n_theta <- 6
  ip <- itempool(a = rlnorm(n_items, 0, .25), b = rnorm(n_items))
  resp <- matrix(rep(sim_resp(ip, theta = rnorm(1, 0, .25)), n_theta),
                 nrow = n_theta, byrow = TRUE)
  for (i in 2:n_theta)  # Impose some missingness to data.
    resp[i, 1:(i-1)] <- NA
  est <- est_ability(ip = ip, resp = resp, method = "ml")
  # Only the first and second row should be the same
  expect_equal(sum(duplicated(est$se)), 0)
  # The following fails when all correct responses.
  # expect_equal(sum(duplicated(est$est)), 0)
  # # Standard errors should increase because smaller number of items
  # se <- est$se
  # for (i in 3:N) expect_lt(se[i-1], se[i])

  # -------------------------------------------------------------------------- #
  # If response string is all missing, i.e. all NA, then ability estimate and
  # standard error should be NA as well.
  est <- est_ability(ip = item(b = rnorm(1)), resp = NA, method = 'owen')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = item(b = rnorm(1)), resp = NA, method = 'ml')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = itempool(b = rnorm(10)), resp = rep(NA, 10),
                           method = 'ml')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))

  # -------------------------------------------------------------------------- #
  # If ip is not an item or Itempool object, it will be converted to an
  # Itempool object:
  ip <- list(item(b = rnorm(1)), item(b = rnorm(1)))
  resp <- sim_resp(ip = ip, theta = rnorm(10))
  expect_type(est_ability(ip = ip, resp = resp, method = "ml")$est, "double")

  # -------------------------------------------------------------------------- #
  # Function can deal with missing responses without an issue. If there are
  # missing responses they should be ignored.
  # Generate an item pool
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
  # Test for ML
  est_small <- est_ability(ip = ip_small, resp = resp_small, method = 'ml')
  est_na <- est_ability(ip = ip, resp = resp_na, method = 'ml')
  expect_equivalent(est_small$est, est_na$est)
  expect_equivalent(est_small$se, est_na$se)
  # Test for EAP
  est_small <- est_ability(ip = ip_small, resp = resp_small, method = 'eap')
  est_na <- est_ability(ip = ip, resp = resp_na, method = 'eap')
  expect_equivalent(est_small$est, est_na$est, tolerance = 1e-10)
  expect_equivalent(est_small$se, est_na$se, tolerance = 1e-10)
  # Test for Owen
  est_small <- est_ability(ip = ip_small, resp = resp_small, method = 'owen')
  est_na <- est_ability(ip = ip, resp = resp_na, method = 'owen')
  expect_equivalent(est_small$est, est_na$est, tolerance = 1e-10)
  expect_equivalent(est_small$se, est_na$se, tolerance = 1e-10)

  ## Missing response to one item should return NA ##
  est <- est_ability(ip = generate_item(), resp = NA, method = 'ml')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = generate_item(), resp = NA, method = 'eap')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = generate_item(), resp = NA, method = 'owen')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))


  ## All Missing responses should return NA ##
  n <- sample(10:20, 1)
  ip <- generate_ip(n = n, model = "2PL")
  resp <- rep(NA, n)
  est <- est_ability(ip = ip, resp = resp, method = 'ml')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))
  est <- est_ability(ip = ip, resp = resp, method = 'owen')
  expect_true(is.na(est$est))
  expect_true(is.na(est$se))

  ### Mixture of Models ###
  # -------------------------------------------------------------------------- #
  # Function can deal with the ability estimation of mixture of item models
  # Here the item pool consist of two 2PL, three GRM and two GPCM items.
  # This test checks whether the log likelihood is calculated correctly.
  D = 1.702
  ip = itempool(data.frame(a = c(0.983, 1.03), b = c(-0.581, 1.643)),
                   D = D, model = '2PL', id = paste0('irt', 1:2))
  # GRM
  item1 = item(a = 1.027, b = c(-2.439, -0.21, 0.693), D = D, model = 'GRM')
  item2 = item(a = 1.024, b = c(0.241, 0.257, 0.311), D = D, model = 'GRM')
  item3 = item(a = 0.896, b = c(-0.225, 0.458, 0.764), D = D, model = 'GRM')
  ip = c(ip, itempool(c(item1, item2, item3), id = paste0('grm', 1:3)))
  item1 = item(a = 0.927, b = c(2.086, 0.107), D = D, model = 'GPCM')
  item2 = item(a = 1.201, b = c(-0.349, 1.162), D = D, model = 'GPCM')
  ip = c(ip, itempool(c(item1, item2), id = paste0('gpcm', 1:2)))
  resp = c(0, 0, 3, 2, 1, 0, 2)
  est <- est_ability(ip = ip, resp = resp, method = 'ml')
  expect_false(is.na(est$est))
  expect_false(is.na(est$se))
  # Sum score calculated correctly
  expect_equal(est_ability(ip = ip, resp = resp, method = 'sum_score')$est,
               sum(resp, na.rm = TRUE))

  # -------------------------------------------------------------------------- #
  # Owen's ability estimation cannot be used models other than dichotomous items
  expect_error(est_ability(ip = ip, resp = resp, method = 'owen'),
               regexp = paste0("Owen's Bayesian ability estimation method can ",
                               "only be used for dichotomous IRT models"))

  # Sum score calculated correctly
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .7, 1.5), b = runif(n, -2, 2))
  resp <- sim_resp(ip = ip, theta = rnorm(8, 0, .4))
  resp[sample(1:prod(dim(resp)), 9)] <- NA
  expect_equal(est_ability(ip = ip, resp = resp, method = 'sum_score')$est,
               rowSums(resp, na.rm = TRUE))


  # -------------------------------------------------------------------------- #
  # If a response string is all NA's, the 'est' and 'se' should be NA as well.
  n_items <- 5
  n_examinee <- 4
  ip <- itempool(a = runif(n_items, .7, 1.5), b = runif(n_items, -2, 2))
  resp <- sim_resp(ip = ip, theta = rnorm(n_examinee, 0, .4))
  resp[1, ] <- 0
  resp[2, ] <- NA
  resp[3, sample(1:n_items, 2)] <- NA

  # Test for all estimation methods
  for (method in c("sum_score", "ml", "eap", "owen")) {
    output <- est_ability(ip = ip, resp = resp, method = method)
    expect_true(is.na(output$est[2]))
    expect_true(is.na(output$se[2]))
  }


  # -------------------------------------------------------------------------- #
  # Ability Estimation with mixed model items
  i1 <- item(b = rnorm(1), id = "i1")
  i2 <- item(a = rlnorm(1, 0, .3), b = rnorm(1), c = .2, id = "i2")
  i3 <- item(a = rlnorm(1, 0, .3), b = sort(runif(3)), id = "i3")
  ip <- c(i1, i2, i3)
  n_examinee <- 4
  resp <- sim_resp(ip = ip, theta = rnorm(n_examinee))
  observed <- est_ability(ip = ip, resp = resp, method = "ml")
  observed <- est_ability(ip = ip, resp = resp, method = "eap")

  # -------------------------------------------------------------------------- #
  # Ability Estimation with the testlets
  t1 <- testlet(itempool(b = rnorm(2), id = c("t1-i1", "t1-i2")),
                   id = "t1")
  t2 <- testlet(itempool(a = rlnorm(3, 0, .3), b = rnorm(3),
                                id = c("t2-i1", "t2-i2", "t2-i3")), id = "t2")
  i1 <- item(b = rnorm(1), id = "i1")
  i2 <- item(a = rlnorm(1, 0, .3), b = rnorm(1), c = .2, id = "i2")
  i3 <- item(a = rlnorm(1, 0, .3), b = sort(runif(3)), id = "i3")
  ip <- c(t1, t2, i1, i2, i3)
  ip1 <- c(t1@item_list, t2@item_list, i1, i2, i3)
  n_examinee <- 4
  resp <- sim_resp(ip = ip, theta = rnorm(n_examinee))
  # ML
  observed <- est_ability(ip = ip, resp = resp, method = "ml")
  expected <- est_ability(ip = ip1, resp = resp, method = "ml")
  expect_equal(expected, observed)
  # # EAP
  # observed <- est_ability(ip = ip, resp = resp, method = "eap")
  # expected <- est_ability(ip = ip1, resp = resp, method = "eap")
  # expect_equal(expected, observed)

  # -------------------------------------------------------------------------- #
  # The estimate numbers are rounded to the tolerance level.
  n_items <- sample(5:10, 1)
  ip <- itempool(b = rnorm(n_items))
  resp <- rep(1, times = n_items)
  est <- est_ability(ip = ip, resp = resp, method = "ml", tol = 10e-13)$est
  for (tol in 10^(-seq(1, 7, 1))) {
    est1 <- est_ability(ip = ip, resp = resp, method = "ml", tol = tol)$est
    expect_equivalent(est, est1, tol = tol)
    # cat(tol, est, est1, "\n", sep = "\t")
  }
  resp <- sim_resp(ip = ip, theta = rnorm(1))
  est <- est_ability(ip = ip, resp = resp, method = "ml", tol = 10e-13)$est
  for (tol in 10^(-seq(1, 7, 1))) {
    est1 <- est_ability(ip = ip, resp = resp, method = "ml", tol = tol)$est
    expect_equivalent(est, est1, tol = tol)
    # cat(tol, est, est1, "\n", sep = "\t")
  }

  # -------------------------------------------------------------------------- #
  # Ability Estimation with the testlets

})

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% est_ability_owen %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

############################# est_ability_owen (Item) ##########################
test_that("est_ability_owen (Item)", {
  ##  Single theta
  est <- est_ability(ip = item(b = -1.480974, D = 1.611941), resp = 1,
                     method = 'owen', prior_pars = c(1, sqrt(.5)))
  expect_equal(est$est, 1.045586, check.attributes = FALSE, tolerance = 1e-5)
  expect_equal(est$se, sqrt(0.4796707), check.attributes = FALSE, tolerance = 1e-5)

  est <- est_ability(ip = item(b = -0.6397636, D = 1.166008), resp = 1,
                     method = 'owen', prior_pars = c(1, sqrt(.5)))
  expect_equal(est$est, 1.080177, check.attributes = FALSE, tolerance = 1e-5)
  expect_equal(est$se, sqrt(0.458222), check.attributes = FALSE, tolerance = 1e-5)


  # -------------------------------------------------------------------------- #
  # Another single item example
  m0 <- -1.055991
  v0 <- 1.424225
  item <- item(a = 1.4792, b = 0.8699, c = 0.1767892, D = 1)
  est <- est_ability(ip = item, resp = 0, method = 'owen',
                     prior_pars = c(m0, sqrt(v0)))
  expect_equal(est$est, -1.224033, tol = 1e-5)
  expect_equal(est$se, sqrt(1.150978), tol = 1e-5)

  # -------------------------------------------------------------------------- #
  ##  Multiple theta raises error
  ip <- item(b = rnorm(1))
  expect_error(
    est_ability(ip = ip, resp = sim_resp(ip = ip, theta = rnorm(10)),
                method = 'owen'),
    regexp = paste0("The length of the response pattern should be ",
                    "equal to the number of items."))
})

############################# est_ability_owen (Itempool) ######################
test_that("est_ability_owen (Itempool)", {
  # Single Item
  est <- est_ability(ip = itempool(a = 0.6129655, b = -0.660839, c = 0.01146235,
                                      D = 2.960664), resp = 1,
                     method = 'owen', prior_pars = c(1, sqrt(.5)))
  expect_equal(est$est, 1.059813, check.attributes = FALSE, tolerance = 1e-5)
  expect_equal(est$se, sqrt(0.4943381), check.attributes = FALSE, tolerance = 1e-5)

  # -------------------------------------------------------------------------- #
  # One item and multiple theta raises error when resp is not a matrix.
  ip <- itempool(b = rnorm(1))
  resp <- as.vector(sim_resp(ip = ip, theta = rnorm(10)))
  expect_error(est_ability(ip = ip, resp = resp, method = 'owen'),
               regexp = paste0("The length of the response pattern should be ",
                               "equal to the number of items."))
  # If resp is a matrix, no error will be raised
  expect_equal(length(est_ability(ip = ip, resp = matrix(resp, ncol = 1),
                                  method = 'owen')$est), 10)

  # -------------------------------------------------------------------------- #
  # Multiple Items
  ip <- itempool(
    a = c(1.34730418561958, 0.838341692695394, 0.670091423671693,
          0.910351292463019, 0.841347561450675, 1.19060949212871,
          1.00658598938026),
    b = c(1.90739340746829, -0.108147490438633, -0.797692163077999,
          -0.0608525829373425, -1.4585962922829, 0.864645431892862,
          -0.26818135533447),
    c = c(0.152733101486228, 0.168994629196823, 0.160417800606228,
          0.0402568220626563, 0.0588782917242497, 0.174779037642293,
          0.274252452305518),
    D = 1.715379, id = paste0("Item-",1:7), content = rep("Algebra", 7))
  resp <- c(0, 1, 1, 1, 1, 1, 1)
  est <- est_ability(ip = ip, resp = resp, method = 'owen',
                     prior_pars = c(1, sqrt(.5)))
  expect_equal(est$est, 1.298185, check.attributes = FALSE, tolerance = 1e-5)
  expect_equal(est$se, sqrt(0.3571816), check.attributes = FALSE, tolerance = 1e-5)

  # -------------------------------------------------------------------------- #
  # One can either use an informative prior obtained by previous responses
  # 1 to n-1 or all item responses from 1 to n.
  m0 <- rnorm(1)
  v0 <- runif(1, .5, 2)
  ip <- generate_ip(n = 10)
  resp <- sim_resp(ip = ip, theta = rnorm(1))[1,]
  temp_est <- est_ability(ip = ip, resp = resp, method = "owen",
                          prior_pars = c(m0, sqrt(v0)))
  # add a new item and response
  item <- generate_item()
  resp_new <- sample(0:1, 1)
  # Two types of estimae: all items vs informative prior. Both should be the
  # same.
  est1 <- est_ability(ip = item, resp = resp_new, method = "owen",
                      prior_pars = c(temp_est$est, sqrt(temp_est$se^2)))
  est2 <- est_ability(ip = c(ip, item), resp = c(resp, resp_new),
                      method = "owen", prior_pars = c(m0, sqrt(v0)))
  expect_equal(est1$est, est2$est, tol = 1e-5)
  expect_equal(est1$se, est2$se, tol = 1e-5)

})


###############################################################################@
############################# est_ability_ml ###################################
###############################################################################@
test_that("est_ability_ml", {
  ip <- itempool(
    a = c(0.884197125327773, 0.641740406746976, 0.637108941096812,
          1.20244530553464, 1.34323987562675, 1.474029551493, 1.42986021412071,
          0.564554519834928, 0.579683377640322, 0.667258570319973,
          1.34760631772224, 1.075640355004, 1.45692114031408, 0.607514199451543,
          0.666346082813106, 1.8072933035437, 1.30535634083208, 1.93822752987035,
          1.73588501708582, 1.60583402158227, 1.09267979150172, 0.540868601528928,
          1.54821544664446, 1.16912602749653, 1.09414781618398, 1.40611228963826,
          0.920170411816798, 1.88622517616022, 0.689107569050975, 0.861420561093837,
          0.525530818500556, 0.844860291108489, 1.73652094230056, 1.23069458303507,
          1.73296662454959, 1.64101757179014, 1.76887892151717, 0.639569879742339),
    b = c(-0.383924056368274, -1.53373225397389, -1.97827922812951, 0.51093189825887,
          -1.481247183071, 0.282873812035544, -1.50036749678062, -1.14698898614752,
          0.604867489398726, 1.39114876765744, -0.323809346783167, -0.762786729612589,
          0.449696800365297, 1.99408201166263, 0.263998428457945, 2.68167705276291,
          -0.780200871008166, 0.424690332217461, 0.629369734172185, 1.52097072996385,
          -1.23893482185556, -1.72437664723673, -1.6663494787618, 0.1500752471024,
          0.573747122347429, -0.694487697707527, 0.101807324722603, 0.923059152892394,
          -1.08240830943348, 0.228139311222657, 0.634357465817595, -0.324786430696751,
          1.20426233805235, -0.422862774355197, 0.588777621523254, -0.546680192400087,
          0.970457541815255, -0.199916207031426),
    c = c(0.159727674443275, 0.181966272555292, 0.132189845200628, 0.0441575736505911,
          0.207619160204194, 0.145414250926115, 0.224728510645218, 0.273771667736582,
          0.192767260316759, 0.137916084146127, 0.280767961777747, 0.0396517782937735,
          0.0979830316035077, 0.294935363531113, 0.094270147732459, 0.197115818294697,
          0.16417963730637, 0.102108501945622, 0.225032159895636, 0.125950632593594,
          0.210073017189279, 0.123149601859041, 0.0160024092765525, 0.214372928347439,
          0.118859681440517, 0.088402352691628, 0.089805198716931, 0.238103704573587,
          0.0694027771940455, 0.100334205455147, 0.259730131994002, 0.268449394055642,
          0.159935772232711, 0.275540588493459, 0.0413019084138796, 0.113490730267949,
          0.142325833230279, 0.135479651740752),
    D = 1
  )
  resp <- c(1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1,
            1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1)
  est <- est_ability(ip = ip, resp = resp, method = 'ml')
  expect_equal(est$est, 0.135097, tolerance = 1e-4, check.attributes = FALSE)
  expect_equal(est$se, 0.3734483, tolerance = 1e-4, check.attributes = FALSE)

  # -------------------------------------------------------------------------- #
  # Multiple response patterns
  theta <- rnorm(3)
  resp <- sim_resp(ip = ip, theta = theta)
  est <- est_ability(ip = ip, resp = resp, method = 'ml')
  expect_equal(length(est$est), 3)
  expect_equal(length(est$se), 3)

  # -------------------------------------------------------------------------- #
  # One item and multiple theta raises error when resp is not a matrix.
  ip <- itempool(b = rnorm(1))
  resp <- as.vector(sim_resp(ip = ip, theta = rnorm(10)))
  expect_error(est_ability(ip = ip, resp = resp, method = 'ml'),
               regexp = paste0("The length of the response pattern should be ",
                               "equal to the number of items."))
  # If resp is a matrix, no error will be raised
  expect_equal(length(est_ability(ip = ip, resp = matrix(resp, ncol = 1),
                                  method = 'ml')$est), 10)

  # -------------------------------------------------------------------------- #
  # Example from Baker Ch.5, p.86-88.
  # http://echo.edres.org:8080/irt/baker/chapter5.pdf
  ip <- itempool(a = c(1, 1.2, .8), b = c(-1, 0, 1), D = 1)
  resp <- c(1, 0, 1)
  est <- est_ability(ip = ip, resp = resp, method = 'ml')
  expect_equal(est$est, .3249, tolerance = 1e-4, check.attributes = FALSE)
  expect_equal(est$se, 1.23, tolerance = 1e-2, check.attributes = FALSE)


  # -------------------------------------------------------------------------- #
  # For all correct and incorrect responses, the difference between the
  # estimate and the theta boundary should be smaller than the tolerance
  n_item <- sample(10:20, 1)
  ip <- itempool(a = rlnorm(n_item, 0, .3), b = rnorm(n_item), D = 1)
  resp1 <- rep(1, n_item)
  resp0 <- rep(0, n_item)
  tol <- runif(1, 1e-14, 1e-2)
  min_theta <- -4
  max_theta <- 4
  est1 <- est_ability(ip = ip, resp = resp1, method = 'ml',
                      theta_range = c(min_theta, max_theta), tol = tol)
  expect_true(max_theta - est1$est < tol)
  est0 <- est_ability(ip = ip, resp = resp0, method = 'ml',
                      theta_range = c(min_theta, max_theta), tol = tol)
  expect_true(min_theta - est0$est < tol)

  # -------------------------------------------------------------------------- #
  # A response pattern with a decimal number is acceptable.
  ip <- itempool(b = c(-0.442, -0.1172, 0.12614, -0.1912, -0.0452, -0.0771,
                           0.12607, -0.1798, 0.03754, -0.283, -0.2239),
                     D = 1)
  resp <- rep(1, 11)
  for (i in 1:length(resp)) {
    resp2 <- resp
    resp2[i] <- 0.5
    est <- est_ability(ip = ip, resp = resp2, method = 'ml')
    # cat("Responses: ",  paste0(resp2, collapse = ", "), "\n")
    # cat(paste0("Theta Estimate = ",  est$est, ",  SE = ", est$se, "\n\n" ))
    expect_equal(est$est, 2.94, tolerance = 1e-2, check.attributes = FALSE)
  }


  # -------------------------------------------------------------------------- #
  # A test case
  ip <- itempool(a = c(1.6807, 1.4203, 1.3582, 0.8569, 0.7989, 1.1182, 0.8105),
                 b = c(0.2829, 1.0154, 0.8887, 1.4957, 1.4821, 0.3831, 0.6588),
                 c = c(0.0970, 0.0971, 0.0748, 0.2116, 0.2184, 0.0163, 0.0174),
                 D = 1)
  resp <- c(1, 1, 1, 1, 0, 1, 1)
  # plot_resp_loglik(ip = ip, resp = resp)
  est <- est_ability(ip = ip, resp = resp, method = "ml")
  expect_equal(est$est, 2.72882)
  expect_equal(est$se, 1.240068)
})


###############################################################################@
############################# abilityEstEAP ####################################
###############################################################################@
test_that("abilityEstEAP", {
  # -------------------------------------------------------------------------- #
  a = c(0.604, 1.945, 1.174, 0.577, 0.722, 1.265, 0.837, 1.752, 1.409, 1.73,
        0.781, 0.527, 1.767, 0.911, 1.648, 1.703, 1.027, 0.947)
  b = c(-1.88, -1.274, -0.993, -0.956, -0.869, -0.621, -0.612, -0.42, -0.393,
        -0.357, -0.232, 0.033, 0.197, 0.331, 0.437, 0.486, 0.804, 1.885)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b), D = D, model = '2PL',
                   id = paste0('i', 1:18))
  resp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  expected <- 1.566382
  object <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_equivalent(object$est, expected, tolerance = 1e-3)

  # -------------------------------------------------------------------------- #
  a = c(0.684, 0.618, 0.824, 1.698, 1.334, 1.179, 1.697, 1.898, 0.772, 0.702,
        0.897, 0.6, 1.881, 1.655, 1.102, 1.514, 0.777, 1.605)
  b = c(-1.841, -1.799, -1.767, -1.629, -0.991, -0.781, -0.688, -0.574, -0.367,
        -0.328, -0.212, 0.044, 0.13, 0.195, 1.213, 1.355, 1.702, 1.893)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b), D = D, model = '2PL',
                   id = paste0('i', 1:18))
  resp = c(1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0)
  expected <- 0.383088
  object <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_equivalent(object$est, expected, tolerance = 1e-3)

  # # -------------------------------------------------------------------------- #
  # ip_dtf <- data.frame(
  #   a = c(0.82628, 1.18005, 0.93197, 1.09279, 0.62252, 0.55811, 0.54016,
  #         0.76499, 1.10399, 0.6466,  0.58406, 0.93681, 0.85042, 0.58623,
  #         0.47125, 0.57258),
  #   b = c(0.3248, 1.58142, 0.9219, 1.84474, 1.25307, 0.9381, 0.40202, 1.15291,
  #         2.75536, 1.45367, 1.67435, 1.9451, 1.20863, 1.66219, 1.87079,
  #         1.35354),
  #   c = c(0, 0.15767, 0.2311, 0.27658, 0, 0.13462, 0.12783, 0.27903, 0.32545,
  #         0.20548, 0.24894, 0.12961, 0.13346, 0.21813, 0, 0))
  # ip <- itempool(ip_dtf, D = 1.7)
  # resp <- rep(1, nrow(ip_dtf))
  # expected <- 2.5227
  # observed <- est_ability(ip = ip, resp = resp, method = 'eap',
  #                         theta_range = c(-3, 3),
  #                         number_of_quads = 10, tol = 0.0001)
  # expect_equivalent(observed$est, expected, tolerance = 1e-3)

  # -------------------------------------------------------------------------- #
  a = c(1.753, 1.295, 1.712, 1.408, 0.9, 0.568, 1.145, 1.407, 1.818, 0.842,
        1.404, 1.086, 0.999, 0.926, 1.689, 0.929, 1.597)
  b = c(-1.968, -1.663, -1.566, -1.454, -0.993, -0.883, -0.73, -0.573, -0.172,
        0.352, 0.947, 1.085, 1.192, 1.245, 1.689, 1.708, 1.713)
  c = c(0.334, 0.172, 0.131, 0.224, 0.297, 0.113, 0.315, 0.084, 0.312, 0.271,
        0.076, 0.082, 0.2, 0.223, 0.293, 0.077, 0.124)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c), D = D, model = '3PL',
                   id = paste0('i', 1:17))
  resp = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1)
  expected <- 1.9033609
  object <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_equivalent(object$est, expected, tolerance = 1e-3)
  expect_equivalent(object$se, 0.4766623, tolerance = 1e-3)

  # -------------------------------------------------------------------------- #
  # Check very low ability is limited
  ip <- itempool(b = rep(-5, 20), D = 1.702)
  resp <- rep(0, 20)
  expected <- -4.223
  est <- est_ability(ip, resp, method = "eap", theta_range = c(expected, 5))
  expect_equivalent(est$est, expected, tolerance = 1e-3)

  # -------------------------------------------------------------------------- #
  # This case seems problematic. Need double checking.
  prior_mean = 0
  prior_sd = 2
  a = c(1.48, 1.15, 1.03, 0.57, 0.69, 0.87, 0.57, 1.27)
  b = c(1.02, 1.17, -0.97, -1.95, -0.66, 0.59, -1.47, 2.09)
  c = c(0.17, 0.17, 0.21, 0.1, 0.24, 0.09, 0.27, 0.29)
  resp = c(0, 0, 1, 0, 1, 0, 0, 1)
  D = 1
  ip = itempool(a = a, b = b, c = c, D = D)
  est = est_ability(ip = ip, resp = resp, method = 'eap', number_of_quads = 91,
                    prior_pars = c(prior_mean, prior_sd), theta_range = c(-5, 5))
  expect_equivalent(est$est, -1.505777, tolerance = 1e-3)


  # -------------------------------------------------------------------------- #
  # For a repeating response string, all of the ability estimates should be
  # the same.
  n <- sample(10:20, 1)
  ip <- itempool(b = rnorm(8))
  resp <- matrix(rep(sim_resp(ip = ip, theta = rnorm(1)), n * 2), nrow = n * 2,
                 byrow = TRUE)
  est <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_true(length(unique(est$est)) == 1)
  expect_true(length(unique(est$se)) == 1)


  # -------------------------------------------------------------------------- #
  # Multiple response patterns
  n <- 18
  theta <- rnorm(n)
  ip <- itempool(b = rnorm(8))
  resp <- sim_resp(ip = ip, theta = theta)
  est <- est_ability(ip = ip, resp = resp, method = 'eap')
  expect_equal(length(est$est), n)
  for (i in 1:n) {
    temp_est <- est_ability(ip = ip, resp = resp[i,], method = 'eap')
    expect_equivalent(temp_est$est, est$est[i])
    expect_equivalent(temp_est$se, est$se[i])
  }
  irt:::est_ability_eap_cpp(resp = resp, ip = ip, theta_range = c(-4, 4),
                        no_of_quadrature = 41,
                        prior_dist = 'norm', prior_par = c(0,1))

  # -------------------------------------------------------------------------- #
  # One item and multiple theta raises error when resp is not a matrix.
  ip <- itempool(b = rnorm(1))
  resp <- as.vector(sim_resp(ip = ip, theta = rnorm(10)))
  expect_error(est_ability(ip = ip, resp = resp, method = 'eap'),
               regexp = paste0("The length of the response pattern should be ",
                               "equal to the number of items."))
  # If resp is a matrix, no error will be raised
  expect_equal(length(est_ability(ip = ip, resp = matrix(resp, ncol = 1),
                                  method = 'eap')$est), 10)

  # -------------------------------------------------------------------------- #
  # A response string with missing responses can be handled correctly.
  n_items <- 10
  ip <- itempool(a = rlnorm(n_items), b = rnorm(n_items))
  resp <- sim_resp(ip = ip, theta = rnorm(1, 0, .2))[1, ]
  missing_indices <- rep(F, n_items)
  missing_indices[sample(1:n_items, 3)] <- TRUE
  resp[missing_indices] <- NA
  ip_new <- ip[!missing_indices]
  resp_new <- resp[!missing_indices]
  est1 <- est_ability(ip = ip, resp = resp, method = "eap")
  est2 <- est_ability(ip = ip_new, resp = resp_new, method = "eap")
  expect_equivalent(est1$est, est2$est)
  expect_equivalent(est1$se, est2$se)

})

