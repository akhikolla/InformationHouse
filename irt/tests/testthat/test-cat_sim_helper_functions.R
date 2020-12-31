

############################################################################@###
################### print.cat_design ###########################################
############################################################################@###

test_that("Test 'print.cat_design' function.", {
  # The following function tests cat_design print to the console.
  # cd: cat_design object.
  test_print_output <- function(cd) {
    p <- capture_output(print(cd))
    expect_output(print(p), sprintf("Item Pool Size: %d", length(cd$ip)))
    expect_output(print(p), sprintf("Maximum Test Length: %d", cd$max_test_length))
    expect_output(print(p), sprintf("First Item Rule: '%s'", cd$first_item_rule))
    expect_output(print(p), sprintf("First Item Parameters:"))
    for (i in seq_len(length(cd$first_item_par)))
      expect_output(print(p), paste0(names(cd$first_item_par)[i], ": ",
                                     cd$first_item_par[[i]]))
    expect_output(print(p), sprintf("Next Item Rule: '%s'",
                                    cd$step[[1]]$next_item_rule))
    expect_output(print(p), sprintf("Ability Estimation Rule: '%s'",
                                    cd$step[[1]]$ability_est_rule))
    expect_output(print(p), "Test Termination Rules and Parameters")
    for (i in seq_len(length(cd$termination_rule)))
      expect_output(print(p), paste0(
        names(cd$termination_par)[i], ": ",
        cd$termination_par[[names(cd$termination_par)[i]]] ))
  }
  n_items = 30
  ip <- itempool(data.frame(a = runif(n_items, .5, 1.5), b = rnorm(n_items)))
  cd <- create_cat_design(ip = ip, next_item_rule = "random",
                          termination_rule = "min_item",
                          termination_par = list(min_item = n_items))
  test_print_output(cd)
  cd <- create_cat_design(ip = ip, next_item_rule = "random")
  test_print_output(cd)

  cd <- create_cat_design(
    ip = ip,
    termination_rule = c("min_item", "min_se"),
    termination_par = list(min_item = 10, min_se = .33))
  test_print_output(cd)

  cd <- create_cat_design(
    ip = itempool(b = rnorm(20)),
    termination_rule = c("min_item", "max_item"),
    termination_par = list(min_item = 20, max_item = 20),
    next_item_rule = "fixed",
    next_item_par = lapply(paste0("Item-", 1:20), function(x) list(item_id = x)))
  test_print_output(cd)

  }) # End of test_that


############################################################################@###
################### summary.cat_output #########################################
############################################################################@###
test_that("Test summary.cat_output Function", {
  n <- 3 # number of items
  ip <- itempool(data.frame(b = rnorm(n)))

  cd = create_cat_design(ip = ip, next_item_rule = "mfi",
                       termination_rule = "max_item",
                       termination_par = list(max_item = n))
  # CAT summary works with single CAT output.
  t1 = cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(summary(t1), "data.frame")

  # CAT summary works with multiple CAT output.
  t2 = cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(summary(list(t1, t2)), "data.frame")

  # Error raised when inadmissable column name entered.
  expect_error(summary(t1, cols = c("true_theta", "somexyx")))

  # -------------------------------------------------------------------------- #
  # Test multilple theta
  nsim <- 5
  true_theta <- rnorm(nsim)
  t1 <- cat_sim(true_ability = true_theta, cd = cd)
  e <- summary(t1)
  expect_is(e, "data.frame")
  expect_equal(nrow(e), nsim)
  expect_equal(colnames(e),
               c("true_ability", "est_ability", "se", "test_length"))
  expect_equal(e$true_ability, true_theta)
  expect_true(all(e$test_length == n))

  desired_cols <- c("true_ability", "est_ability", "se", "test_length", "bias")
  e <- summary(t1, cols = desired_cols)
  expect_equal(colnames(e), desired_cols)
  expect_equal(e$true_ability, true_theta)
  expect_true(all(e$test_length == n))
  expect_equal(e$bias, e$est_ability - e$true_ability)
})



############################################################################@###
################### print.cat_output ###########################################
############################################################################@###

test_that("Test 'print.cat_output' function.", {
  n <- 3 # number of items
  ip <- itempool(data.frame(b = rnorm(n)))

  cd <- create_cat_design(ip = ip, next_item_rule = "mfi",
                          termination_rule = "max_item",
                          termination_par = list(max_item = n))
  # CAT summary works with single CAT output.
  co <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_output(expected <- print(co))
  expect_is(expected, "data.frame")
  expect_is(expected$item_id, "character")
  # When silent, it should not produce output
  expect_silent(print(co, silent = TRUE))
})


############################################################################@###
################### get_cat_response_data ######################################
############################################################################@###
test_that("Test get_cat_response_data Function", {
  # Function works with single CAT output.
  n <- 10 # number of items
  ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i-", 1:n))

  cd <- create_cat_design(ip = ip, next_item_rule = "random",
                         termination_rule = "max_item",
                         termination_par = list(max_item = n))
  co <- cat_sim(true_ability = rnorm(1), cd = cd)
  resp <- get_cat_response_data(cat_sim_output = co, cd = cd)
  expect_is(resp, "integer")
  for (i in seq_len(n)) {
    expect_equal(co$est_history[[i]]$resp, as.numeric(resp[i]))
    expect_equal(co$est_history[[i]]$item$id, names(resp[i]))
  }
  # ---------------------------------- #
  # Function works when "cd" not provided
  resp <- get_cat_response_data(cat_sim_output = co)
  expect_is(resp, "integer")
  for (i in seq_len(n)) {
    expect_equal(co$est_history[[i]]$resp, as.numeric(resp[i]))
    expect_equal(co$est_history[[i]]$item$id, names(resp[i]))
  }
  # ---------------------------------- #
  # Attaching CAT summary
  resp <- get_cat_response_data(cat_sim_output = co, attach_summary = TRUE)
  expect_is(resp, "data.frame")
  expect_equal(nrow(resp), 1)
  expect_true(all(c("true_ability", "est_ability", "se", "test_length") %in%
                    colnames(resp)))
  expect_true(all(ip$id %in% colnames(resp)))

  # -------------------------------------------------------------------------- #
  n <- 5 # number of items
  ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i-", 1:n))

  cd <- create_cat_design(ip = ip, next_item_rule = "random",
                          termination_rule = "max_item",
                          termination_par = list(max_item = 3))
  # Function works with multiple CAT output.
  co <- cat_sim(true_ability = rnorm(2), cd = cd)
  t1 <- co[[1]]
  t2 <- co[[2]]
  resp <- get_cat_response_data(cat_sim_output = co, cd = cd)
  expect_is(resp, "data.frame")
  for (i in 1:3) {
    expect_equal(t1$est_history[[i]]$resp, resp[1, t1$est_history[[i]]$item$id])
    expect_equal(t2$est_history[[i]]$resp, resp[2, t2$est_history[[i]]$item$id])
  }
  # ---------------------------------- #
  # Function works when "cd" not provided
  resp <- get_cat_response_data(cat_sim_output = co)
  expect_is(resp, "data.frame")
  for (i in 1:3) {
    expect_equal(t1$est_history[[i]]$resp, resp[1, t1$est_history[[i]]$item$id])
    expect_equal(t2$est_history[[i]]$resp, resp[2, t2$est_history[[i]]$item$id])
  }
  # ---------------------------------- #
  # Attaching CAT summary
  resp <- get_cat_response_data(cat_sim_output = co, attach_summary = TRUE)
  expect_true(all(c("true_ability", "est_ability", "se", "test_length") %in%
                    colnames(resp)))

  # ---------------------------------- #
  # Attaching CAT summary for a list
  resp <- get_cat_response_data(cat_sim_output = list(t1, t2),
                                attach_summary = TRUE)
  expect_is(resp, "data.frame")
  expect_true(all(c("true_ability", "est_ability", "se", "test_length") %in%
                    colnames(resp)))
  expect_true(sum(ip$id %in% colnames(resp)) %in% 3:5)
  # -------------------------------------------------------------------------- #

})

############################################################################@###
################### get_cat_administered_items #################################
############################################################################@###
test_that("Test get_cat_administered_items Function", {
  # Function works with single CAT output.
  n <- 10 # number of items
  ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i-", 1:n))

  cd <- create_cat_design(ip = ip, next_item_rule = "random",
                         termination_rule = "max_item",
                         termination_par = list(max_item = n))
  co <- cat_sim(true_ability = rnorm(1), cd = cd)
  expect_is(get_cat_administered_items(co), "Itempool")

  # -------------------------------------------------------------------------- #
  # Function works with multiple CAT outputs.
  ip <- itempool(data.frame(b = rnorm(n)), id = paste0("i-", 1:10))

  cd <- create_cat_design(ip = ip, next_item_rule = "random",
                          termination_rule = "max_item",
                          termination_par = list(max_item = 3))
  n_theta <- sample(2:6, 1)
  # Function works with multiple CAT output.
  co <- cat_sim(true_ability = rnorm(n_theta), cd = cd)
  administered_ip <- get_cat_administered_items(co)
  expect_equal(length(administered_ip), n_theta)
  expect_true(all(sapply(administered_ip, is, "Itempool")))
})

############################################################################@###
################### calculate_exposure_rates ###################################
############################################################################@###
test_that("Test calculate_exposure_rates Function", {
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")),
                   id = "t1")
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")),
                   id = "t2")
  t3 <- testlet(itempool(b = 0:1, id = c("t3-i1", "t3-i2")),
                   id = "t3")
  ip <- c(t1, t2, t3, itempool(b = rnorm(18)))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          termination_rule = "max_item",
                          termination_par = list("max_item" = 10))
  co <- cat_sim(true_ability = rnorm(3), cd = cd)
  observed <- calculate_exposure_rates(co, cd)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))

  # Run the same function with item_ids
  observed <- calculate_exposure_rates(co, item_ids = cd$ip$id)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))

  # cd <- create_cat_design(
  #   ip = ip,
  #   next_item_rule = "fixed",
  #   next_item_par = lapply(c("t3", "Item-1", "Item-3", "Item-11"),
  #                          function(x) list(item_id = x)),
  #   termination_rule = "max_item", termination_par = list("max_item" = 5))


})


############################################################################@###
################### calculate_overlap_rates ####################################
############################################################################@###
test_that("Test calculate_overlap_rates Function", {
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")))
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")))
  ip <- c(t1, t2, itempool(b = rnorm(18)))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          termination_rule = "max_item",
                          termination_par = list(max_item = 10))
  co <- cat_sim(true_ability = rnorm(5), cd = cd)
  observed <- calculate_overlap_rates(co, cd)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))

  # Run the same function with item_ids
  observed <- calculate_overlap_rates(co, item_ids = cd$ip$id)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))
})



