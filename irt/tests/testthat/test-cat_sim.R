# library(microbenchmark, testthat)

# library(testthat)
# test_file(path = file.path(getwd(), 'tests', 'testthat', 'test-cat_sim.R'))

############################################################################@###
##################____ create_cat_design ____###################################
############################################################################@###
test_that("Test 'create_cat_design' function.", {
  ip <- generate_ip(n = 55)
  # Check the default:
  cd <- create_cat_design()
  expect_is(cd, "cat_design")
  # FIX ME:
  # expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  cd <- create_cat_design(ip = ip)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # Title can be entered
  title_text <- "xyz123abc"
  cd <- create_cat_design(ip = generate_ip(), title = title_text,
                          next_item_rule = 'random',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = 3))
  expect_is(cd, "cat_design")
  expect_equal(cd$title, title_text)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ####------------------ Item Pool Tests -----------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ### Item Pool Tests ###
  # Item pool size cannot be smaller than maximum item number.
  expect_error(create_cat_design(
    ip = itempool(data.frame(b = rnorm(5))),
    termination_rule = c('min_item', 'min_se', 'max_item'),
    termination_par = list(min_item = 10, min_se = .33, max_item = 20)))

  # -------------------------------------------------------------------------- #
  # The size of the true Itempool should be the same as ip.
  true_ip <- itempool(data.frame(a = runif(10, .5, 1.5), b = rnorm(10)))
  expect_error(create_cat_design(
    ip = itempool(data.frame(b = rnorm(5))),
    true_ip = true_ip,
    termination_rule = c('min_item', 'min_se', 'max_item'),
    termination_par = list(min_item = 10, min_se = .33, max_item = 20)))

  # -------------------------------------------------------------------------- #
  # true_ip works as it should
  n_ip <- 55
  true_ip <- itempool(data.frame(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip)))
  cd <- create_cat_design(ip = ip, true_ip = true_ip)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # id's of true_ip and ip should be the same.
  true_ip <- itempool(data.frame(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip)),
                      id = c(paste0(1:(n_ip - 1)), "last_item"))
  expect_error(create_cat_design(ip = ip, true_ip = true_ip))

  # -------------------------------------------------------------------------- #
  # A very short length test
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # Check the small item pool error:
  expect_error(create_cat_design(ip = ip),
               "Item pool size should not be smaller")

  # -------------------------------------------------------------------------- #
  # 'ip' should be a valid item pool object.
  expect_error(create_cat_design(
    ip = list(),
    termination_rule = c('min_item', 'min_se', 'max_item'),
    termination_par = list(min_item = 10, min_se = .33, max_item = 20)))

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Termination Rule ----------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # Unmatched list sizes
  n_ip <- 55
  ip <- itempool(data.frame(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip)))
  expect_error(create_cat_design(
    termination_rule = rep(list(c('min_item', 'min_se', 'max_item')), 10),
    termination_par = rep(list(list(min_item = 10, min_se = .33,
                                    max_item = 20)), 20)))
  # No error expected
  cd <- create_cat_design(
    ip = ip, termination_rule = c('min_item', 'min_se', 'max_item'),
    termination_par = list(min_item = 10, min_se = .33, max_item = 20))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # Maximum item number can be larger than the length(ip) if there are testlets
  temp_list <- list(ids = paste0("testlet-", 1:3), n = c(2, 3, 4))
  ip <- c(generate_ip(n = 5), itempool(sapply(1:length(temp_list$id), function(i)
    generate_testlet(id = temp_list$id[i], n = temp_list$n[i]))))
  cd <- create_cat_design(ip = ip, termination_rule = c('max_item'),
                          termination_par = list(max_item = 10))
  expect_is(cd, "cat_design")
  # FIX ME
  # expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Next Item Rule ------------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ### Next Item Rule ###
  n_ip <- 55
  ip <- generate_ip(n = n_ip)
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          next_item_par = NULL)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # Next item rule cannot be NULL
  expect_error(create_cat_design(next_item_rule = NULL))

  # -------------------------------------------------------------------------- #
  # "fixed" items selection rule
  cd <- create_cat_design(
    ip = itempool(b = rnorm(20)),
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 20, max_item = 20),
    next_item_rule = 'fixed',
    next_item_par = lapply(paste0("Item-", 1:20), function(x) list(item_id = x)))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # "fixed" selection rule needs valid item pool
  expect_error(create_cat_design(
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 20, max_item = 20),
    next_item_rule = 'fixed',
    next_item_par = lapply(paste0(1:20), function(x) list(item_id = x))),
    "There should be a valid item pool")

  # -------------------------------------------------------------------------- #
  # 'next_item_par' for "fixed" selection rule item parameter can be NULL.
  # In that case, items in the item pool will be administered in that order.
  temp_list <- list(ids = paste0("testlet-", 1:8), n = 1:8)
  ip <- itempool(sample(c(generate_ip(n = 10, output = "list"),
                          sapply(1:length(temp_list$id), function(i)
                            generate_testlet(id = temp_list$id[i],
                                             n = temp_list$n[i])))))

  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 10, max_item = 10),
    next_item_rule = 'fixed',
    next_item_par = NULL)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 8, max_item = 8),
    next_item_rule = 'fixed',
    next_item_par = list(item_id = "testlet-8"))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 8, max_item = 8),
    next_item_rule = 'fixed',
    next_item_par = list(item_id = c("Item-1", "testlet-8")))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  cd <- create_cat_design(
    ip = itempool(b = rnorm(20)),
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 20, max_item = 20),
    next_item_rule = 'fixed',
    next_item_par = NULL)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # sum_score
  cd <- create_cat_design(ip = generate_ip(n = 30),
                          next_item_rule = "fixed",
                          ability_est_rule = "sum_score",
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 30))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Ability Estimation Rule ---------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  ### Ability Estimation Rule ###
  # expect_is(create_cat_design(ip = ip, ability_est_rule = 'sum_score',
  #                             ability_est_par = NULL),
  #           "cat_design")
  cd <- create_cat_design(
    ip = generate_ip(n = 50),
    ability_est_rule = 'eap',
    ability_est_par = list(prior_dist = 'unif',
                           prior_par = list(min = -2, max = 2),
                           min_theta = -4, max_theta = 4,
                           no_of_quadrature = 30))
  expect_is(cd, "cat_design")
  # FIX ME:  WARNING the following line crashes R.
  # expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # If no ability_est_par given, defaults will be used for "eap"
  cd <- create_cat_design(ip = generate_ip(n = 50), ability_est_rule = 'eap')
  expect_true(all(sapply(cd$step,
                         function(x) x$ability_est_par$prior_dist) == "norm"))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$min_theta) == -4))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$max_theta) == 4))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$prior_par)[1] == 0))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$prior_par)[2] == 1))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$no_of_quadrature) == 50))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # If no ability_est_par given, defaults will be used for "owen"
  cd <- create_cat_design(ip = generate_ip(n = 50), ability_est_rule = 'owen')

  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$prior_mean) == 0))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$prior_var) == 1))
  # The function defaults should be overridden.
  expect_true(all(sapply(
    cd$step, function(x) is.null(x$ability_est_par$no_of_quadrature))))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # If no ability_est_par given, defaults will be used for "ml"
  cd <- create_cat_design(ip = generate_ip(n = 50), ability_est_rule = 'ml')
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$min_theta) == -4))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$max_theta) == 4))
  expect_true(
    all(sapply(cd$step, function(x) x$ability_est_par$criterion) == 0.001))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # 'ml' with all relevant parameters
  cd <- create_cat_design(ip = generate_ip(n = 50), ability_est_rule = 'ml',
                          ability_est_par = list(min_theta = -3, max_theta = 3,
                                                 criterion = 0.01))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # 'ml' with incorrect parameters should raise an error.
  expect_error(create_cat_design(
    ip = generate_ip(n = 50), ability_est_rule = 'ml',
    ability_est_par = list(m_theta = -3, x_theta = 3, criterion = 0.01)))

  # -------------------------------------------------------------------------- #
  # FIX ME: The lack of item pool should give a better error before running
  # 'cat_sim()' function
  cd <- create_cat_design(ip = generate_ip(n = 50), ability_est_rule = 'ml')
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Final Ability Estimation Rule ---------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # -------------------------------------------------------------------------- #
  # If no final_ability_est_par given, defaults will be used for "eap"
  cd <- create_cat_design(ip = generate_ip(n = 50),
                          final_ability_est_rule = 'eap')

  expect_true(cd$final_ability_est_par$prior_dist == "norm")
  expect_true(cd$final_ability_est_par$min_theta == -4)
  expect_true(cd$final_ability_est_par$max_theta == 4)
  expect_true(cd$final_ability_est_par$prior_par[1] == 0)
  expect_true(cd$final_ability_est_par$prior_par[2] == 1)
  expect_true(cd$final_ability_est_par$no_of_quadrature == 50)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # If no final_ability_est_par given, defaults will be used for "owen"
  cd <- create_cat_design(ip = generate_ip(n = 50),
                          final_ability_est_rule = 'owen')
  expect_true(cd$final_ability_est_par$prior_mean == 0)
  expect_true(cd$final_ability_est_par$prior_var == 1)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # If no final_ability_est_par given, defaults will be used for "ml"
  cd <- create_cat_design(ip = generate_ip(n = 50),
                          final_ability_est_rule = 'ml')
  expect_true(cd$final_ability_est_par$min_theta == -4)
  expect_true(cd$final_ability_est_par$max_theta == 4)
  expect_true(cd$final_ability_est_par$criterion == 0.001)
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Exposure Control ----------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ### Exposure Control ###
  cd <- create_cat_design(ip = generate_ip(n = 50),
                          exposure_control_rule = 'randomesque',
                          exposure_control_par = list(num_items = 1))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # 5-4-3-2-1 exposure control
  cd <- create_cat_design(
    ip = generate_ip(n = 50), exposure_control_rule = 'randomesque',
    exposure_control_par = lapply(c(5:1, rep(1, 15)),
                                  function(x) list(num_items = x)))
  expect_is(cd, "cat_design")
  expect_is(cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")

  # -------------------------------------------------------------------------- #
  # Randomesque without a parameter should throw an error
  expect_error(create_cat_design(exposure_control_rule = 'randomesque',
                               exposure_control_par = NULL))

  # -------------------------------------------------------------------------- #
  # A warning will be issued if the number of sympson-hetter-k values that are
  # equal to 1 are smaller than the maximum test length.
  n_ip <- 10
  ip <- generate_ip(n = n_ip,
                    misc = lapply(c(0, runif(n_ip, 0.4, 0.6)),
                                  function(k) list("sympson_hetter_k" = k)))
  expect_warning(create_cat_design(ip = ip,
                          termination_rule = c('max_item'),
                          termination_par = list(max_item = n_ip - 2),
                          next_item_rule = 'mepv',
                          exposure_control_rule = "sympson-hetter",
                          ),
                 "When using 'sympson-hetter' exposure control rule")






  # skip("Check the following test in the future.")
  # # A problem is with sympson hetter, the parameter is not acceptable
  # # because of tibble. if it is wrapped with "as.numeric()" it works, see ip3.
  # # In create_cat_design put a check for this.
  # n <- 10
  # ip_dtf <- dplyr::tibble(a = rlnorm(n, 0, 0.3), b = rnorm(n),
  #                         K = runif(n, .5, 1))
  # ip1 <- list()
  # for (i in seq_len(n)) {
  #   ip1 <- c(ip1, item(a = ip_dtf[i, "a"], b = ip_dtf[i, "b"],
  #                      misc = list(sympson_hetter_k = ip_dtf[i, "K"])))
  # }
  # ip2 <- list()
  # for (i in seq_len(n)) {
  #   ip2 <- c(ip2, item(a = as.numeric(ip_dtf[i, "a"]),
  #                      b = as.numeric(ip_dtf[i, "b"]),
  #                      misc = list(sympson_hetter_k = ip_dtf[i, "K"])))
  # }
  # ip3 <- list()
  # for (i in seq_len(n)) {
  #   ip3 <- c(ip3, item(
  #     a = as.numeric(ip_dtf[i, "a"]), b = as.numeric(ip_dtf[i, "b"]),
  #     misc = list(sympson_hetter_k = as.numeric(ip_dtf[i, "K"]))))
  # }
  # cd <- create_cat_design(
  #   ip = itempool(ip2),
  #   first_item_rule = "fixed_theta",
  #   first_item_par = list(theta = 0),
  #   next_item_rule = "mepv",
  #   termination_rule = c("max_item"),
  #   termination_par = list(max_item = 5),
  #   exposure_control_rule = "sympson-hetter"
  # )
  # true_ability <- rnorm(1)
  # expected <- cat_sim(true_ability = true_ability, cd = cd)
  #
  # cd <- create_cat_design(
  #   ip = ip3,
  #   first_item_rule = "fixed_theta",
  #   first_item_par = list(theta = 0),
  #   next_item_rule = "mepv",
  #   termination_rule = c("max_item"),
  #   termination_par = list(max_item = 5),
  #   exposure_control_rule = "sympson-hetter"
  # )
  # expected <- cat_sim(rnorm(1), cd = cd)
  #
  # ### Content Balancing ###
  # # expect_is(create_cat_design(
  # #   content_bal_rule = 'max_discrepancy',
  # #   content_bal_par = list(target_dist = c(
  # #     Geometry = .3, `Rational Numbers` = .2, Algebra = .5))),
  # #   "cat_design")

  }) # End of test_that



############################################################################@###
#################____ cat_sim ____##############################################
############################################################################@###
test_that("Test cat_sim Function", {

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- cat_sim Function Tests ----------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Test cat_sim for single true_ability
  n <- 100 # number of items
  ip <- itempool(a = runif(n, .5, 1.5), b = sort(rnorm(n)), c = runif(n, 0,.3))
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 5))
  true_ability <- rnorm(1)
  expected <- cat_sim(true_ability = true_ability, cd = cd)
  expect_is(expected, "cat_output")
  expect_equivalent(expected$true_ability, true_ability)
  expect_equal(length(expected$est_history), 5)

  # -------------------------------------------------------------------------- #
  # Test for multiple theta
  n <- 100 # number of items
  ip <- itempool(a = runif(n, .5, 1.5), b = sort(rnorm(n)),
                     c = runif(n, 0,.3))
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 5))
  true_ability <- rnorm(sample(5:20, 1))
  expected <- cat_sim(true_ability = true_ability, cd = cd)
  expect_is(expected, "list")
  expect_true(all(sapply(expected, inherits, "cat_output")))
  expect_equal(length(expected), length(true_ability))


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Ability Estimation Rule ---------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Simulation with ml (both intermediate and final ability est = "ml")
  n <- sample(20:40, 1) # number of items
  ip <- generate_ip(model = "3PL", n = n)
  min_theta <- runif(1, -6, -2)
  max_theta <- runif(1, 2, 6)
  max_item <- sample(5:12, 1)
  cd <- create_cat_design(
    ip = ip, next_item_rule = 'mfi',
    ability_est_rule = "ml",
    ability_est_par = list(min_theta = min_theta, max_theta = max_theta,
                           criterion = 0.001),
    final_ability_est_rule = "ml",
    final_ability_est_par = list(min_theta = min_theta, max_theta = max_theta,
                                 criterion = 0.001),
    termination_rule = 'max_item',
    termination_par = list(max_item = max_item))

  co <- cat_sim(true_ability = rnorm(1), cd = cd)
  # Since it is ML first item should be equal to either max or min theta.
  expect_true(co$est_history[[1]]$est_after %in% c(min_theta, max_theta))
  administered_ip <- itempool(sapply(co$est_history, function(x) x$item))
  resp <- sapply(co$est_history, function(x) x$resp)

  # The ability estimate just before the administration of the last item is
  # indeed calculated correctly:
  expect_equal(
    est_ability(administered_ip[-max_item], resp[-max_item], method = "ml")$est,
    co$est_history[[max_item]]$est_before, tol = 0.001)
  # Final ability estimate calculated correctly:
  expect_equal(
    est_ability(administered_ip, resp, method = "ml")$est,
    co$est_history[[max_item]]$est_after, tol = 0.001)

  # -------------------------------------------------------------------------- #
  # Simulation with testlets - eap
  n_testlet <- 10
  ip_list <- list()
  for (i in seq_len(n_testlet)) {
    temp <- 2
    ip_list <- c(ip_list, testlet(itempool(b = rnorm(temp),
                                       id = paste0("t", i, "-i", 1:temp)),
                          id = paste0("t", i)))
  }
  # Add individual items
  ip <- c(itempool(ip_list), itempool(b = rnorm(n_testlet)))
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 20))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 20)

  # -------------------------------------------------------------------------- #
  # Simulation with testlets - owen
  n_testlet <- 10
  ip_list <- list()
  for (i in seq_len(n_testlet)) {
    temp <- 2
    ip_list <- c(ip_list, testlet(itempool(
      b = rnorm(temp), id = paste0("t", i, "-i", 1:temp)),
      id = paste0("t", i)))
  }
  # Add individual items
  ip <- c(itempool(ip_list), itempool(b = rnorm(n_testlet)))
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0, prior_var = 4),
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 20))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 20)

  # -------------------------------------------------------------------------- #
  # sum_score
  ip <- generate_ip(n = 30)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "fixed",
                          ability_est_rule = "sum_score",
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 30))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 30)
  # An intermediate ability estimate is correctly calculated
  k <- sample(2:28, 1)
  expect_equal(sum(sapply(test$est_history, `[[`, "resp")[1:k]),
               test$est_history[[k]]$est_after)
  # Final ability estimate is correctly calculated
  expect_equal(sum(sapply(test$est_history, `[[`, "resp")),
               test$est_history[[30]]$est_after)


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Next Item Rule ------------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Test for a linear test where items order is predefined.
  n <- 5
  ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i",1:n))
  cd = create_cat_design(
    ip = ip,
    next_item_rule = 'fixed',
    next_item_par = lapply(c("i3", "i2", "i4", "i5", "i1"),
                           function(x) list(item_id = x)),
    # next_item_par = list(list(item_id = 'i3'), list(item_id = 'i2'),
    #                      list(item_id = 'i4'), list(item_id = 'i5'),
    #                      list(item_id = 'i1')),
    ability_est_rule = "eap",
    termination_rule = 'max_item', termination_par = list(max_item = n))
  test = cat_sim(true_ability = rnorm(1), cd = cd)
  expect_equal(test$est_history[[1]]$item$id, 'i3')
  expect_equal(test$est_history[[2]]$item$id, 'i2')
  expect_equal(test$est_history[[3]]$item$id, 'i4')
  expect_equal(test$est_history[[4]]$item$id, 'i5')
  expect_equal(test$est_history[[5]]$item$id, 'i1')
  # summary(test)

  # cat_sim returns a list when true_ability is a vector of numbers
  test = cat_sim(true_ability = rnorm(3), cd = cd)
  expect_is(test, 'list')
  expect_equal(length(test), 3)


  # -------------------------------------------------------------------------- #
  # Fixed where ip consist of standalone items and testlets and
  # next_item_par = NULL
  temp_list <- list(ids = paste0("testlet-", 1:8), n = 1:8)
  ip <- itempool(c(generate_ip(n = 5, output = "list"),
                   sapply(1:length(temp_list$id), function(i)
                     generate_testlet(id = temp_list$id[i],
                                      n = temp_list$n[i]))))
  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 10, max_item = 10),
    next_item_rule = 'fixed',
    next_item_par = NULL)
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')


  # -------------------------------------------------------------------------- #
  # Fixed where ip consist of standalone items and testlets.
  temp_list <- list(ids = paste0("testlet-", 1:8), n = 1:8)
  ip <- itempool(sample(c(generate_ip(n = 10, output = "list"),
                          sapply(1:length(temp_list$id), function(i)
                            generate_testlet(id = temp_list$id[i],
                                             n = temp_list$n[i])))))
  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 8, max_item = 8),
    next_item_rule = 'fixed',
    next_item_par = list(item_id = "testlet-8"))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  # All selected testlets are from 'testlet-8'
  expect_true(all(sapply(lapply(test$est_history, `[[`, "testlet"), identical,
                         ip["testlet-8"][[1]])))

  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 8, max_item = 8),
    next_item_rule = 'fixed',
    next_item_par = list(item_id = c("Item-1", "testlet-8")))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  eh <- test$est_history
  expect_equal(eh[[1]]$item, ip["Item-1"][[1]])
  expect_equal(eh[[2]]$testlet, ip["testlet-8"][[1]])
  expect_equal(eh[[2]]$item, ip["testlet-8"][[1]]$item_list[[1]])

  # -------------------------------------------------------------------------- #
  # Fixed where ip consist of standalone items and next_item_par is NULL.
  # Items are delivered as the items appear in the item pool
  ip <- generate_ip(n = 20)
  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 20, max_item = 20),
    next_item_rule = 'fixed',
    next_item_par = NULL)
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  n <- sample(1:20, 1)
  expect_equal(test$est_history[[n]]$item, ip[[n]])

  # -------------------------------------------------------------------------- #
  # Fixed where ip consist of standalone items and next_item_par is NULL.
  # Items are delivered as the items appear in the item pool
  ip <- generate_ip(n = 20)
  cd <- create_cat_design(
    ip = ip,
    termination_rule = c('min_item', 'max_item'),
    termination_par = list(min_item = 10, max_item = 15),
    next_item_rule = 'fixed',
    next_item_par = NULL)
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  n <- sample(1:15, 1)
  expect_equal(test$est_history[[n]]$item, ip[[n]])


  # -------------------------------------------------------------------------- #
  # Simulation with Rasch model - eap-mepv
  ip <- itempool(b = rnorm(10), model = "Rasch")
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 5))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_false(is.na(test$est_history[[2]]$est_after))

  # -------------------------------------------------------------------------- #
  # Simulation with Rasch model - owen-mepv
  ip <- itempool(b = rnorm(10), model = "Rasch")
  cd <- create_cat_design(ip = ip, next_item_rule = 'mepv',
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0, prior_var = 4),
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 5))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_false(is.na(test$est_history[[2]]$est_after))


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Final Ability Estimation Rule ---------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Simulation with different final ability estimate
  n <- 20 # number of items
  ip <- itempool(a = runif(n, .5, 1.5), b = sort(rnorm(n)),
                     c = runif(n, 0,.3))
  cd <- create_cat_design(
    ip = ip, next_item_rule = 'mfi',
    ability_est_rule = "owen",
    ability_est_par = list(prior_mean = 0, prior_var = 4),
    final_ability_est_rule = "eap",
    final_ability_est_par = list(
      prior_dist = "norm", prior_par = c(0, 1), min_theta = -4, max_theta = 4,
      no_of_quadrature = 50),
    termination_rule = 'max_item',
    termination_par = list(max_item = 5))
  true_ability <- rnorm(sample(5:20, 1))
  expected <- cat_sim(true_ability = true_ability, cd = cd)
  expect_is(expected, "list")
  expect_true(all(sapply(expected, inherits, "cat_output")))
  expect_equal(length(expected), length(true_ability))


  # -------------------------------------------------------------------------- #
  # Simulation with ml and eap (intermediate "ml", final ability est = "eap")
  n <- sample(20:40, 1) # number of items
  ip <- generate_ip(model = "3PL", n = n)
  min_theta <- runif(1, -6, -2)
  max_theta <- runif(1, 2, 6)
  max_item <- sample(5:12, 1)
  cd <- create_cat_design(
    ip = ip, next_item_rule = 'mfi',
    ability_est_rule = "ml",
    ability_est_par = list(min_theta = min_theta, max_theta = max_theta,
                           criterion = 0.001),
    final_ability_est_rule = "eap",
    final_ability_est_par = list(
      prior_dist = "norm", prior_par = c(0, 1), min_theta = min_theta,
      max_theta = max_theta, no_of_quadrature = 50),
    termination_rule = 'max_item',
    termination_par = list(max_item = max_item))

  co <- cat_sim(true_ability = rnorm(1), cd = cd)
  # Since it is ML first item should be equal to either max or min theta.
  expect_true(co$est_history[[1]]$est_after %in% c(min_theta, max_theta))
  administered_ip <- itempool(sapply(co$est_history, function(x) x$item))
  resp <- sapply(co$est_history, function(x) x$resp)

  # The ability estimate just before the administration of the last item is
  # indeed calculated correctly:
  est_r <- est_ability(administered_ip[-max_item], resp[-max_item],
                       method = "ml")$est
  est_cpp <- co$est_history[[max_item]]$est_before
  if (abs(est_r - est_cpp) > 0.001) {
    cat("\n\nML estiamte for the following item within CPP is not equal to ",
        "the ML estimate in R's est_ability() function. \n\n")
    print(administered_ip[-max_item])
    cat("\nResponses: ", resp[-max_item], "\n\n")
  }
  expect_equal(est_r, est_cpp, tol = 0.001)
  # Final ability estimate calculated correctly:
  expect_equal(
    est_ability(administered_ip, resp, method = "eap",
                theta_range = c(min_theta, max_theta))$est,
    co$est_history[[max_item]]$est_after, tol = 0.001)


  # -------------------------------------------------------------------------- #
  # sum_score
  ip <- generate_ip(n = 30)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "fixed",
                          ability_est_rule = "eap",
                          final_ability_est_rule = "sum_score",
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 30))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 30)
  # Final ability estimate is correctly calculated
  expect_equal(sum(sapply(test$est_history, `[[`, "resp")),
               test$est_history[[30]]$est_after)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Exposure Control ----------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Simulation with testlets with owen + sympson-hetter
  n_testlet <- 25
  ip_list <- list()
  for (i in seq_len(n_testlet)) {
    temp <- 2
    ip_list <- c(ip_list, testlet(itempool(
      b = rnorm(temp), id = paste0("t", i, "-i", 1:temp)),
      id = paste0("t", i), misc = list(
        sympson_hetter_k = ifelse(i %% 2 == 1, 1, runif(1, .5, 1)))))
  }
  # Add individual items
  ip <- c(itempool(ip_list), itempool(
    b = rnorm(n_testlet),
    misc = lapply(c(rep(1, ceiling(n_testlet/2)),
                    round(runif(n_testlet - ceiling(n_testlet/2), 0.5, 1), 2)),
                  function(x) list(sympson_hetter_k = x))
    ))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mepv',
                          exposure_control_rule = "sympson-hetter",
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0, prior_var = 4),
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 20))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 20)

  # -------------------------------------------------------------------------- #
  # Simulation with testlets with owen + sympson-hetter + fixed_theta
  n_theta <- 3 # Number of examinees
  true_theta <- rnorm(n_theta, mean = 0, sd = 1)
  cd <- create_cat_design(
    ip = ip,
    first_item_rule = "fixed_theta",
    first_item_par = list(theta = 0),
    next_item_rule = "mepv",
    ability_est_rule = "owen",
    ability_est_par = list(prior_mean = 0, prior_var = 4),
    termination_rule = c("max_item"),
    termination_par = list(max_item = 10),
    exposure_control_rule = "sympson-hetter"
  )
  expected <- cat_sim(true_ability = true_theta, cd = cd)
  expect_equal(length(expected), n_theta)
  expect_true(all(sapply(expected, is, 'cat_output')))

  # -------------------------------------------------------------------------- #
  # Sympson-Hetter will not select an item with 0 values
  n_ip <- 4
  ip <- generate_ip(n = n_ip,
                    misc = lapply(c(0, rep(1, n_ip)),
                                  function(k) list("sympson_hetter_k" = k)))
  ip[[1]]@parameters$a <- 3 # Make the item very eligible
  cd <- create_cat_design(ip = ip,
                          termination_rule = c('max_item'),
                          termination_par = list(max_item = 3),
                          next_item_rule = 'mepv',
                          exposure_control_rule = "sympson-hetter",
                          )
  expect_is(cd, "cat_design")
  expect_is(co <- cat_sim(true_ability = rnorm(1), cd = cd), "cat_output")
  expect_false("Item-1" %in% sapply(co$est_history, function(k) k$item$id))
  expect_true("Item-2" %in% sapply(co$est_history, function(k) k$item$id))

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Termination Rule ----------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  n <- 100 # number of items
  ip <- itempool(data.frame(a = runif(n, .5, 1.5), b = sort(rnorm(n)),
                               c = runif(n, 0,.3)), id = paste0("Item-",1:n),
                    content = sample(c("Algebra", "Arithmetic", "Geometry"),
                                     n, replace = TRUE))

  cd <- create_cat_design(ip = ip, next_item_rule = 'mfi',
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 5))
  test <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(test, 'cat_output')
  expect_equal(length(test$est_history), 5)

  # -------------------------------------------------------------------------- #
  ### SPRT ###
  ip <- generate_ip(n = sample(30:50, 1))
  sprt_par <- list(theta_0 = -.3, theta_1 = .3, alpha = .1, beta = .1)
  A <- (1 - sprt_par$beta) / sprt_par$alpha
  B <- sprt_par$beta/(1-sprt_par$alpha)
  cd <- create_cat_design(ip = ip, termination_rule = c("max_item", "sprt"),
                          termination_par = list("sprt" = sprt_par,
                                                 "max_item" = sample(17:25, 1)))
  co <- cat_sim(true_ability = runif(1, 1, 3), cd = cd)
  administered_ip <- get_cat_administered_items(co)
  resp <- get_cat_response_data(co)

  # Make sure that the ratio of likelihood (and log-likelihood) are outside the
  # interval [B, A] at the end of the CAT test.
  lik0 <- resp_lik(ip = administered_ip, resp = resp, theta = sprt_par$theta_0)
  lik1 <- resp_lik(ip = administered_ip, resp = resp, theta = sprt_par$theta_1)
  lik_ratio <- lik0/lik1
  expect_true(lik_ratio < B | lik_ratio > A)
  # Make sure that the ratio of likelihood (and log-likelihood) are WITHIN the
  # interval [B, A] just before the administration of the last item.
  test_length <- length(co$est_history)
  lik0 <- resp_lik(ip = administered_ip[-test_length],
                   resp = resp[-test_length], theta = sprt_par$theta_0)
  lik1 <- resp_lik(ip = administered_ip[-test_length],
                   resp = resp[-test_length], theta = sprt_par$theta_1)
  lik_ratio <- lik0/lik1
  expect_false(lik_ratio < B | lik_ratio > A)


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # -------------------- Misc Scenarios ------------------------------------####
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #





})


############################################################################@###
################### select_next_item #######################################@###
############################################################################@###
# test_that("Test select_next_item Function", {
#   skip("Skipping test of 'select_next_item' Function, R implementation")
#   n <- 25 # number of items
#   ip <- itempool(data.frame(a = runif(n, 1, 1),
#                                b = seq(from = -3, to = 3, length.out = n),
#                                c = runif(n, 0, 0) ),
#                     id = paste0("Item-",1:n),
#                     content = sample(c("Algebra", "Arithmetic", "Geometry"),
#                                      n, replace = TRUE))
#   # CAT Design
#   cd = create_cat_design(ip = ip, next_item_rule = 'mfi',
#                        termination_rule = 'max_item',
#                        termination_par = list(max_item = 5))
#   est_history = list('0' = list(
#     est = 0, se = NA, resp = NULL, item = NULL))
#   item = select_next_item(cd = cd, est_history = est_history)
#   expect_is(item, 'item')
#   expect_equal(item$id,  ip[[13]]$id)
#   est_history$`1` = list(est = NA, se = NA, item = item, resp = 1)
#   est_history$`2` = list(est = NA, se = NA, item = ip[[2]], resp = 0)
#   est_history$`3` = list(est = NA, se = NA, item = ip[[3]], resp = 0)
#   # test <- cat_sim(trueTheta = rnorm(1), ip = ip, maxTestLength = 20)
#
#
#   # Random Test
#   n_ip = 5
#   ip = itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
#   cd = create_cat_design(ip = ip, next_item_rule = 'random',
#                          termination_rule = 'max_item',
#                          termination_par = list('max_item' = n_ip))
#   est_history = list("0" = list(est = 0.2, se = 0.3,  resp = 0, item = ip[[1]]),
#                      "1" = list(est = .4, se = .2,  resp = 1, item = NULL))
#   selected_item = select_next_item_cpp(cd = cd, est_history = est_history)
#   expect_true(selected_item$id %in% ip[-1]$id)
#
#
#   # Fixed test
#   n_ip = 5
#   ip = itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
#   cd = create_cat_design(
#     ip = ip,
#     next_item_rule = 'fixed',
#     next_item_par = list(list(item_id = "Item-2"), list(item_id = "Item-5"),
#                          list(item_id = "Item-1"), list(item_id = "Item-3"),
#                          list(item_id = "Item-4")),
#     termination_rule = 'max_item',
#     termination_par = list('max_item' = n_ip),
#     exposure_control_rule = 'randomesque',
#     exposure_control_par = list(num_items = 1))
#
#   est_history = list("0" = list(est = 0.2, se = 0.3,  resp = 0, item = ip[[1]]),
#                      "1" = list(est = .4, se = .2,  resp = 1, item = NULL))
#   expect_equivalent(select_next_item_cpp(cd = cd, est_history = est_history),
#                     ip[[2]])
# })

############################################################################@###
################### generate_cat_response ##################################@###
############################################################################@###
# test_that("Test generate_cat_response Function", {
#   skip("Skipping test of 'generate_cat_response' Function, R implementation")
#   n <- 5 # number of items
#   ip <- itempool(data.frame(b = c(20, 20, 20, 20, 20)), id = paste0("i", 1:n))
#   # CAT Design
#   cd = create_cat_design(ip = ip, next_item_rule = 'random',
#                        termination_rule = 'max_item',
#                        termination_par = list(max_item = 3))
#   est_history = list('0' = list(
#     est = 0, se = NA, resp = NULL, item = NULL))
#   item = select_next_item(cd = cd, est_history = est_history)
#   est_history <- c(est_history, list(list(item = item)))
#   # Since true_ip is not used this should be definitely incorrect
#   expect_equal(generate_cat_response(cd = cd, est_history = est_history,
#                                      true_ability = 0), 0)
#
#   true_ip = itempool(data.frame(b = c(-20, -20, -20, -20, -20)),
#                       id = paste0("i", 1:n))
#   cd = create_cat_design(ip = ip, true_ip = true_ip, next_item_rule = 'random',
#                        termination_rule = 'max_item',
#                        termination_par = list(max_item = 3))
#   expect_equal(generate_cat_response(cd = cd, est_history = est_history,
#                                      true_ability = 0), 1)
# })

############################################################################@###
################### ability_est_cat  #######################################@###
############################################################################@###
# test_that("Test ability_est_cat Function", {
#   skip("Skipping test of 'ability_est_cat' Function, R implementation")
#   n <- 25 # number of items
#   ip_dtf = data.frame(a = 1, b = rnorm(n),
#                       #seq(from = -3, to = 3, length.out = n),
#                       c = 0, d = 1, D = 1.7)
#   ip <- itempool(ip_dtf,
#                     id = paste0("Item-",1:n),
#                     content = sample(c("Algebra", "Arithmetic", "Geometry"),
#                                      n, replace = TRUE))
#   # CAT Design
#   cd = create_cat_design(
#     ip = ip, ability_est_rule = 'eap',
#     ability_est_par = list(prior_dist = 'norm',
#                            prior_par = list(mean = 0, sd = 1),
#                            min_theta = -5, max_theta = 5),
#     termination_rule = 'max_item', termination_par = list(max_item = n))
#
#   resp = sample(c(0,1), n, replace = TRUE)# c(1, 0, 1)
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:n) {
#     est_history[paste(i)] = list(list(est = NA, se = NA, item = ip[[i]],
#                                       resp = resp[i]))
#     x1 = est_ability(resp = resp[1:i], ip = ip[1:i], method = 'eap',
#                      theta_range = c(-5, 5), prior_pars = c(0, 1))
#     x2 = ability_est_cat(cd = cd, est_history = est_history)
#     expect_equal(x1$est, x2$est)
#     expect_equal(x1$se, x2$se)
#   }
#
#   # ------------------------------------------------------------------------ #
#   skip("Skipping test of 'ability estimation' Function, R implementation")
#   # Ability estimation using sum score.
#   n <- 5
#   ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i",1:n))
#   cd = create_cat_design(
#     ip = ip, ability_est_rule = 'sum_score',
#     termination_rule = 'max_item', termination_par = list(max_item = n))
#   est_history = list('0' = list(
#     est = 0, se = NA, resp = NULL, item = NULL))
#   est_history$`1` = list(est = NA, se = NA, item = ip[[1]], resp = 1)
#   expect_equal(ability_est_cat(cd = cd, est_history = est_history)$est, 1)
#   est_history$`2` = list(est = NA, se = NA, item = ip[[2]], resp = 0)
#   expect_equal(ability_est_cat(cd = cd, est_history = est_history)$est, 1)
#   est_history$`3` = list(est = NA, se = NA, item = ip[[3]], resp = 1)
#   expect_equal(ability_est_cat(cd = cd, est_history = est_history)$est, 2)
#   est_history$`4` = list(est = NA, se = NA, item = ip[[4]], resp = 0)
#   expect_equal(ability_est_cat(cd = cd, est_history = est_history)$est, 2)
# })


############################################################################@###
################### terminate_cat ##########################################@###
############################################################################@###
# test_that("Test terminate_cat Function", {
#   skip("Skipping test of 'terminate_cat' Function, R implementation")
#   ### Test should terminate when maximum length reached. ###
#   n <- 15 # number of items
#   ip <- itempool(data.frame(b = rnorm(n)))
#   cd = create_cat_design(
#     ip = ip, termination_rule = 'min_se',
#     termination_par = list(min_se = 0.0001)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:n)
#     est_history[paste(i)] = list(list(est = NA, se = NA, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_true(terminate_cat(cd = cd, est_history = est_history))
#
#   ### Test should not end before minimum items administered, even other
#   # criteria met, like min_se below. ###
#   min_item = 5
#   cd = create_cat_design(
#     ip = ip, termination_rule = c('min_item', 'min_se', 'max_item'),
#     termination_par = list(min_item = min_item, min_se = .33,
#                            max_item = n)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:(min_item - 2))
#     est_history[paste(i)] = list(list(est = NA, se = 0.01, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_false(terminate_cat(cd = cd, est_history = est_history))
#
#   ### The order of criteria is important: ###
#   min_item = 5
#   cd = create_cat_design(
#     ip = ip, termination_rule = c('min_se', 'min_item', 'max_item'),
#     termination_par = list(min_se = .33, min_item = min_item,
#                            max_item = n)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:(min_item - 2))
#     est_history[paste(i)] = list(list(est = NA, se = 0.01, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_true(terminate_cat(cd = cd, est_history = est_history))
#
#   ### Only one termination rule and it is not satisfied, then test should
#   ### continue
#   min_item = 5
#   cd = create_cat_design(
#     ip = ip, termination_rule = c('min_item'),
#     termination_par = list(min_item = min_item)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:(min_item - 2))
#     est_history[paste(i)] = list(list(est = NA, se = 0.01, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_false(terminate_cat(cd = cd, est_history = est_history))
#
#   max_item = 5
#   cd = create_cat_design(
#     ip = ip, termination_rule = c('max_item'),
#     termination_par = list(max_item = max_item)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:(min_item - 2))
#     est_history[paste(i)] = list(list(est = NA, se = 0.01, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_false(terminate_cat(cd = cd, est_history = est_history))
#
#   min_se = .0001
#   cd = create_cat_design(
#     ip = ip, termination_rule = c('min_se'),
#     termination_par = list(min_se = min_se)) # A very low value
#   est_history = list('0' = list(est = 0, se = NA, resp = NULL, item = NULL))
#   for (i in 1:(min_item - 2))
#     est_history[paste(i)] = list(list(est = NA, se = 0.01, item = ip[[i]],
#                                       resp = sample(0:1, 1)))
#   expect_false(terminate_cat(cd = cd, est_history = est_history))
# })

