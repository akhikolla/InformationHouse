# library(rbenchmark, testthat)


#' This function creates an cat_output object in order to ease the writing
#' of est_history that will be used throughout this test module. Ideally
#' one does not need to
create_est_history <- function(num_of_steps, ip = NULL,
                               empty_last_step = TRUE, cat_design = NULL,
                               true_ability = rnorm(1)) {
  cd <- cat_design
  if (is.null(cat_design) && !is.null(ip)) {
    cd <- create_cat_design(
      ip = ip, termination_rule = 'max_item',
      termination_par = list('max_item' = min(length(ip), num_of_steps)))
  }
  co <- cat_sim(true_ability = true_ability, cd = cd)
  # Only get the portion of the estimate history that is declared by number
  # of steps
  co$est_history <- co$est_history[1:num_of_steps]
  if (empty_last_step) {
    i <- length(co$est_history)
    co$est_history[[i]]$resp <- NA
    co$est_history[[i]]$est_after <- NA
    co$est_history[[i]]$se_after <- NA
    co$est_history[[i]]["testlet"] <- list(NULL)
    co$est_history[[i]]["item"] <- list(NULL)
  }
  return(co)
}

print_est_history <- function(est_history, silent = FALSE) {
  eh <- list(est_history = est_history)
  class(eh) <- "cat_output"
  return(irt:::.print.cat_output(eh, silent = silent))
}

############################################################################@###
################### get_remaining_items ########################################
############################################################################@###
test_that("Test get_remaining_items Function", {
  n_ip <- 12
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = n_ip))
  selected_items <- sample(1:n_ip, 2)
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[selected_items[1]]], est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[selected_items[2]]], est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA, item = NULL,
         est_after = 0.3, se_after = .2)
    )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  selected_ids <- ip[selected_items]$id
  expect_false(any(selected_ids %in% remaining_ids))

  # -------------------------------------------------------------------------- #
  ### Testlet ###
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2")), id = "t1")
  t2 <- testlet(itempool(b = 2:4, id = c("t2-i1", "t2-i2", "t2-i3")),
                id = "t2")
  i1 <- item(b = -1, id = "i1")
  i2 <- item(b = 0, id = "i2")
  i3 <- item(b = 1, id = "i3")
  ip <- c(t1, t2, i1, i2, i3)
  expect_equal(length(ip), 5)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          next_item_par = list(var_calc_method = "eap"),
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = 4))
  # Select a new item after an administration of a testlet
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[3]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  expect_true("t1" %in% remaining_ids)
  expect_false("t2" %in% remaining_ids)
  expect_true("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_true("i3" %in% remaining_ids)


  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  remaining_ids <- remaining_ip$id
  expect_true("t1" %in% remaining_ids)
  expect_true("t2" %in% remaining_ids)
  expect_true("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_true("i3" %in% remaining_ids)


  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = i1, testlet = NULL, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[1]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[2]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  expect_false("t1" %in% remaining_ids)
  expect_true("t2" %in% remaining_ids)
  expect_false("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_true("i3" %in% remaining_ids)


  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = i1, testlet = NULL, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[1]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[2]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  expect_false("t1" %in% remaining_ids)
  expect_true("t2" %in% remaining_ids)
  expect_false("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_true("i3" %in% remaining_ids)


  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = i1, testlet = NULL, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[1]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t1$item_list[[2]], testlet = t1, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = i3, testlet = NULL, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  remaining_ip <- irt:::get_remaining_items(cd, est_history, list())
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  expect_false("t1" %in% remaining_ids)
  expect_true("t2" %in% remaining_ids)
  expect_false("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_false("i3" %in% remaining_ids)

  # -------------------------------------------------------------------------- #
  # Check whether the items in the set_aside_item_list are indeed removed from
  # remaining item's list.
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = t2$item_list[[3]], testlet = t2, est_after = -1, se_after = .6),
    list(est_before = -1, se_before = 0.6, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  additional_args <- list(set_aside_item_list = list(i1, t1))
  remaining_ip <- irt:::get_remaining_items(cd, est_history, additional_args)
  expect_is(remaining_ip, "Itempool")
  expect_true(validObject(remaining_ip))
  remaining_ids <- remaining_ip$id
  expect_false("t1" %in% remaining_ids)
  expect_false("t2" %in% remaining_ids)
  expect_false("i1" %in% remaining_ids)
  expect_true("i2" %in% remaining_ids)
  expect_true("i3" %in% remaining_ids)

})

############################################################################@###
################### get_administered_items_cpp #################################
############################################################################@###
test_that("Test get_administered_items_cpp Function", {
  n_ip <- 5
  ip <- itempool(b = -2:2)

  # -------------------------------------------------------------------------- #
  # No item administered yet.
  est_history <- list(
    list(est_before = 1.2, se_before = 0.3, resp = NA, item = NULL,
         est_after = NA, se_after = NA)
  )

  # This should be an empty Itempool object
  administered_ip <- irt:::get_administered_items_cpp(est_history = est_history)
  expect_is(administered_ip, "Itempool")
  # expect_true(validObject(administered_ip))
  expect_equal(length(administered_ip@item_list), 0)

  # -------------------------------------------------------------------------- #
  # Only one item administered
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[3]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = NA, item = NULL,
         est_after = NA, se_after = NA)
  )
  administered_ip <- irt:::get_administered_items_cpp(est_history = est_history)
  expect_is(administered_ip, "Itempool")
  expect_true(validObject(administered_ip))
  expect_equal(length(administered_ip), 1)
  expect_true(ip[[3]]$id %in% administered_ip$id)
  expect_false(ip[[1]]$id %in% administered_ip$id)
  expect_false(ip[[2]]$id %in% administered_ip$id)
  expect_false(ip[[4]]$id %in% administered_ip$id)
  expect_false(ip[[5]]$id %in% administered_ip$id)

  # -------------------------------------------------------------------------- #
  # Two items administered, but no third element added to est_history
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[3]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.2, resp = 1, item = ip[[2]],
         est_after = 2, se_after = -.1)
  )
  administered_ip <- irt:::get_administered_items_cpp(est_history = est_history)
  expect_is(administered_ip, "Itempool")
  expect_true(validObject(administered_ip))
  expect_equal(length(administered_ip), 2)
  expect_true(ip[[3]]$id %in% administered_ip$id)
  expect_true(ip[[2]]$id %in% administered_ip$id)
  expect_false(ip[[1]]$id %in% administered_ip$id)
  expect_false(ip[[4]]$id %in% administered_ip$id)
  expect_false(ip[[5]]$id %in% administered_ip$id)

  # -------------------------------------------------------------------------- #
  # Three items administered
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[3]],
         est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 1, item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[5]],
         est_after = 0.3, se_after = .2),
    list(est_before = -0.5, se_before = 0.3, resp = NA, item = NULL,
         est_after = 0.3, se_after = .2)
  )
  administered_ip <- irt:::get_administered_items_cpp(est_history = est_history)
  expect_is(administered_ip, "Itempool")
  expect_true(validObject(administered_ip))
  expect_equal(length(administered_ip), 3)
  expect_true(ip[[3]]$id %in% administered_ip$id)
  expect_true(ip[[1]]$id %in% administered_ip$id)
  expect_false(ip[[2]]$id %in% administered_ip$id)
  expect_false(ip[[4]]$id %in% administered_ip$id)
  expect_true(ip[[5]]$id %in% administered_ip$id)
})


############################################################################@###
################### loglik_est_history #########################################
############################################################################@###
test_that("Test loglik_est_history Function", {
  ip <- itempool(b = rnorm(10))
  theta <- rnorm(1)
  resp <- sim_resp(ip = ip, theta = rnorm(1))[1, ]

  # -------------------------------------------------------------------------- #
  # Only one item administered
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = resp[1], item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = NA, item = NULL,
         est_after = NA, se_after = NA)
  )
  expect_equal(irt:::loglik_est_history(est_history = est_history,
                                        theta = theta, calculate_loglik = TRUE),
               resp_loglik(ip = ip[[1]], resp = resp[1], theta = theta))
  expect_equal(irt:::loglik_est_history(
    est_history = est_history, theta = theta, calculate_loglik = FALSE),
               resp_lik(ip = ip[[1]], resp = resp[1], theta = theta))

  # -------------------------------------------------------------------------- #
  # Two items administered, but no third element added to est_history
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = resp[1], item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.2, resp = resp[2], item = ip[[2]],
         est_after = 2, se_after = -.1)
  )
  expect_equal(irt:::loglik_est_history(est_history = est_history,
                                        theta = theta, calculate_loglik = TRUE),
               resp_loglik(ip = ip[1:2], resp = resp[1:2], theta = theta))
  expect_equal(irt:::loglik_est_history(
    est_history = est_history, theta = theta, calculate_loglik = FALSE),
    resp_lik(ip = ip[1:2], resp = resp[1:2], theta = theta))

  # -------------------------------------------------------------------------- #
  # four items administered,
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = resp[1], item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.2, resp = resp[2], item = ip[[2]],
         est_after = 2, se_after = -.1),
    list(est_before = 1.2, se_before = 0.2, resp = resp[3], item = ip[[3]],
         est_after = 2, se_after = -.1),
    list(est_before = 1.2, se_before = 0.2, resp = resp[4], item = ip[[4]],
         est_after = 2, se_after = -.1),
    list(est_before = 1.2, se_before = 0.2, resp = NA, item = NULL,
         est_after = 2, se_after = -.1)
  )
  expect_equal(irt:::loglik_est_history(est_history = est_history,
                                        theta = theta, calculate_loglik = TRUE),
               resp_loglik(ip = ip[1:4], resp = resp[1:4], theta = theta))

  expect_equal(irt:::loglik_est_history(est_history = est_history,
                                     theta = theta, calculate_loglik = FALSE),
               resp_lik(ip = ip[1:4], resp = resp[1:4], theta = theta))

  # -------------------------------------------------------------------------- #
  # Error when inadmissible resp or item
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = resp[1], item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = "a", item = ip[[2]],
         est_after = NA, se_after = NA)
  )
  expect_error(
    irt:::loglik_est_history(est_history = est_history, theta = theta),
               "Inadmissable resp value")
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = resp[1], item = ip[[1]],
         est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = NA, item = ip[[2]],
         est_after = NA, se_after = NA)
  )
  expect_error(
    irt:::loglik_est_history(est_history = est_history, theta = theta),
    "Inadmissable resp value")

})


############################################################################@###
################### select_next_item_cpp #######################################
############################################################################@###

test_that("Test select_next_item_cpp Function", {
  # -------------------------------------------------------------------------- #
  # If an item from a testlet is administered, then regardless of the plan,
  # the next item will be the item from that testlet:
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")))
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")))
  ip <- c(t1, t2, itempool(b = -1:0, id = paste0("i", 1:2)))
  t1 <- ip[[1]]
  t2 <- ip[[2]]
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                         termination_rule = 'max_item',
                         termination_par = list('max_item' = 4))
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[4]], testlet = NULL, est_after = 0.3, se_after = .2),
    # As a second item administer the first item in the testlet
    list(est_before = 0.2, se_before = 0.3, resp = 0, item = t1$item_list[[1]],
         testlet = t1, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA, item = NULL,
         testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item[["item"]], t1$item_list[[2]])
  expect_equal(selected_item[["testlet"]], t1)

  # -------------------------------------------------------------------------- #
  ### Random Test Item Selection ###
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                         termination_rule = 'max_item',
                         termination_par = list('max_item' = n_ip))
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_true(selected_item$item$id %in% ip$id)

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_true(selected_item$item$id %in% ip[-1]$id)

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_true(selected_item$item$id %in% ip[-c(1:2)]$id)

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_true(selected_item$item$id %in% ip[-c(1:3)]$id)

  # -------------------------------------------------------------------------- #
  # Random Test Item Selection with Simpson-Hetter Exposure control
  # When there is an item with very low Simpson hetter value, that item will
  # not be administered:
  n_items <- 2000
  test_length <- 3
  ip <- generate_ip(n = n_items, misc = lapply(
    c(rep(1, test_length), rep(0, n_items - test_length)),
    function(x) list(sympson_hetter_k = x)))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "random",
                          exposure_control_rule = "sympson-hetter",
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = test_length))
  # Select a testlet as first item
  est_history <- list(
    list(est_before = -4, se_before = 0.3, resp = 1,
         item = ip[[1]], testlet = NULL, est_after = 0, se_after = 1),
    list(est_before = 0, se_before = 0.3, resp = 1,
         item = ip[[2]], testlet = NULL, est_after = 0, se_after = 1),
    list(est_before = 0, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  for (i in 1:5) {
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, additional_args = list())
    expect_true("set_aside_item_list" %in% names(selected_item$additional_args))
    # Even though there are hundreds of items, due to exposure control only
    # the one eligible item selected.
    expect_equal(selected_item$est_history[[3]]$item, ip[[3]])
  }

  # -------------------------------------------------------------------------- #
  # Fixed test
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(
    ip = ip,
    next_item_rule = 'fixed',
    next_item_par = list(list(item_id = "Item-2"), list(item_id = "Item-5"),
                         list(item_id = "Item-1"), list(item_id = "Item-3"),
                         list(item_id = "Item-4")),
    termination_rule = 'max_item',
    termination_par = list('max_item' = n_ip),
    exposure_control_rule = 'randomesque',
    exposure_control_par = list(num_items = 1))

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 1,
         item = ip[[4]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equivalent(selected_item$item, ip[[5]])

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 1,
         item = ip[[4]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = 0,
         item = ip[[5]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equivalent(selected_item$item, ip[[3]])

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
    next_item_par = list(item_id = c("testlet-8", "Item-1")))

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equivalent(selected_item$testlet, ip["testlet-8"][[1]])
  expect_equivalent(selected_item$item, ip["testlet-8"][[1]]$item_list[[1]])

  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = 1,
         item = ip["testlet-8"][[1]]$item_list[[1]],
         testlet = ip["testlet-8"][[1]], est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equivalent(selected_item$testlet, ip["testlet-8"][[1]])
  expect_equivalent(selected_item$item, ip["testlet-8"][[1]]$item_list[[2]])

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # Maximum Fisher Information
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  n_ip <- 5
  ip <- itempool(b = -2:2)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mfi',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[3]],
         testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[4]]$id)
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[5]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = -0.5, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = 0.3, se_after = .2)
    )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[2]]$id)

  # -------------------------------------------------------------------------- #
  # MFI with testlets should select the most infomative item or testlet
  ip <- c(item(a = 1, b = -2),
          item(a = 1, b = -3),
          testlet(itempool(item(a = 0.5, b = 0),
                           item(a = 0.75, b = -1), id = paste0("t1-i", 1:2))),
          testlet(itempool(item(a = 0.25, b = 0),
                           item(a = 0.50, b = -1),
                           item(a = 0.25, b = 1), id = paste0("t2-i", 1:3))),
          item(a = 1, b = 0)
          )
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mfi',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = 3))
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[3]][[1]],
         testlet = ip[[3]], est_after = 0.3, se_after = .2),
    list(est_before = 1, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )

  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, additional_args = list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  # Second item of the test let is selected
  expect_equal(ip[[3]][[2]], selected_item$item)
  expect_equal(ip[[3]], selected_item$testlet)

  theta <- 1
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[2]],
         testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = theta, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )

  item_info <- info(ip[-2], theta = theta)[1,]
  output <- irt:::select_next_item_fisher_max_info_cpp(
    cd = cd, est_history = est_history, additional_args = NULL)
  expect_equal(output$criteria[1], item_info[which.max(item_info)])
  expect_equal(output$criteria[4], item_info[which.min(item_info)])
  expect_equal(output$remaining_ip_list[[1]], ip[[5]])

  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, additional_args = list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item, ip[[5]])
  expect_null(selected_item$testlet)


  theta <- -1
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0, item = ip[[2]],
         testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = theta, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
  item_info <- info(ip[-2], theta = theta)[1,]
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, additional_args = list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$testlet, ip[[3]])



  # -------------------------------------------------------------------------- #
  ### Minimum Expected Posterior Variance ###
  # Item
  n_ip <- 5
  ip <- itempool(b = -2:2)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          next_item_par = list(var_calc_method = "eap"),
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = NA, item = NULL,
         est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[3]]$id)

  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 1.2, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[2]]$id)
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 1,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[5]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = -0.5, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[2]]$id)
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 1,
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[5]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0, se_before = 0.3, resp = 0,
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = -0.5, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item$id, ip[[4]]$id)

  # -------------------------------------------------------------------------- #
  # MEPV does not change
  n_ip <- 5
  ip <- itempool(a = rlnorm(n_ip, 0, .3), b = rnorm(n_ip),
                     c = runif(n_ip, 0, .3))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          next_item_par = list(var_calc_method = "owen"),
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  est_history <- list(
    list(est_before = 0, se_before = 0.3, resp = NA, item = NULL,
         est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  selected_item <- selected_item[[length(selected_item)]]
  selected_item_id <- selected_item$item$id
  for (i in 1:10) {
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item$id, selected_item_id)
  }

  # -------------------------------------------------------------------------- #
  ### Minimum Expected Posterior Variance ###
  # Item
  n_ip <- 5
  ip <- itempool(a = rlnorm(n_ip, 0, .3), b = rnorm(n_ip),
                     c = runif(n_ip, 0, .3))
  initial_ability_est <- rnorm(1)
  cd <- create_cat_design(ip = ip,
                          first_item_rule = "fixed_theta",
                          first_item_par = list(theta = initial_ability_est),
                          next_item_rule = "mepv",
                          next_item_par = list(var_calc_method = "owen"),
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0, prior_var = 4),
                          termination_rule = c("max_item"),
                          termination_par = list(max_item = n_ip))
  est_history <- list(
    list(est_before = initial_ability_est, se_before = NA, resp = NA,
         item = NULL, est_after = NA, se_after = NA)
  )
  epv <- rep(0, n_ip)
  for (i in seq_len(n_ip)) {
    item <- ip[[i]]
    for (resp in 0:1) {
      P <- prob(ip = item, theta = initial_ability_est)
      if (resp == 0) P <- 1 - P
      est <- irt:::est_ability_owen_cpp(ip = itempool(item), resp = resp,
                                        m0 = 0, v0 = 1)
      epv[i] = epv[i] + P * est$se^2
    }
  }
  expected <- ip$id[which.min(epv)]
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())$est_history
  observed <- selected_item[[length(selected_item)]]$item$id
  expect_equal(expected, observed)

  # -------------------------------------------------------------------------- #
  ### Testlet ###
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2"), D = 1.702),
                id = "t1")
  t2 <- testlet(itempool(b = 2:4, id = c("t2-i1", "t2-i2", "t2-i3"),
                          D = 1.702), id = "t2")
  i1 <- item(b = -1, id = "i1", D = 1.702)
  i2 <- item(b = 0, id = "i2", D = 1.702)
  i3 <- item(b = 1, id = "i3", D = 1.702)
  ip <- c(t1, t2, i1, i2, i3)
  expect_equal(length(ip), 5)
  for (vcm in c("eap", "owen")) {
    cd <- create_cat_design(ip = ip,
                            next_item_rule = "mepv",
                            next_item_par = list(var_calc_method = vcm),
                            termination_rule = 'min_item',
                            termination_par = list('min_item' = 4))
    # Select a testlet as first item
    est_history <- list(
      list(est_before = -4, se_before = 0.3, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, t1@item_list[[1]])
    expect_equal(selected_item$testlet, t1)

    # Select an item as first item
    est_history <- list(
      list(est_before = 0, se_before = 0.3, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, i2)
    expect_null(selected_item$testlet)

    # Select a 3-item testlet as first item
    est_history <- list(
      list(est_before = 4, se_before = 0.3, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, t2@item_list[[1]])
    expect_equal(selected_item$testlet, t2)

    # Select a second item as testlet
    est_history <- list(
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t1$item_list[[1]], testlet = t1, est_after = -1, se_after = .6),
      list(est_before = -1, se_before = 0.6, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, t1@item_list[[2]])
    expect_equal(selected_item$testlet, t1)

    # Select a third item as testlet
    est_history <- list(
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = 0, se_before = 0.3, resp = 1,
           item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = -1, se_before = 0.6, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, t2@item_list[[3]])
    expect_equal(selected_item$testlet, t2)

    # Select a new item after an administration of a testlet
    est_history <- list(
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[3]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = -1, se_before = 0.6, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    if (vcm == "owen") {
      # # Here is how mepv calculated for two items
      # owen1 <- irt:::est_ability_owen_cpp(ip = t2@item_list, resp = c(0, 0, 0),
      #                                     m0 = 0, v0 = 1)
      # # i1
      # resp_lik(ip = i1, resp = 1, theta = -1) * irt:::est_ability_owen_cpp(
      #   ip = i1, resp = 1, m0 = owen1$est, v0 = owen1$se^2)$se^2 +
      # resp_lik(ip = i1, resp = 0, theta = -1) * irt:::est_ability_owen_cpp(
      #   ip = i1, resp = 0, m0 = owen1$est, v0 = owen1$se^2)$se^2
      # # i2
      # resp_lik(ip = i2, resp = 1, theta = -1) * irt:::est_ability_owen_cpp(
      #   ip = i2, resp = 1, m0 = owen1$est, v0 = owen1$se^2)$se^2 +
      # resp_lik(ip = i2, resp = 0, theta = -1) * irt:::est_ability_owen_cpp(
      #   ip = i2, resp = 0, m0 = owen1$est, v0 = owen1$se^2)$se^2
      expect_equal(selected_item$item, t1@item_list[[1]])
      expect_equal(selected_item$testlet, t1)
    } else {
      expect_equal(selected_item$item, i1)
      expect_null(selected_item$testlet)
    }

    # Select a new item after an administration of a testlet
    est_history <- list(
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[1]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[2]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = 0, se_before = 0.3, resp = 0,
           item = t2$item_list[[3]], testlet = t2, est_after = -1, se_after = .6),
      list(est_before = -2.5, se_before = 0.6, resp = NA,
           item = NULL, testlet = NULL, est_after = NA, se_after = NA)
    )
    selected_item <- irt:::select_next_item_cpp(
      cd = cd, est_history = est_history, list())$est_history
    selected_item <- selected_item[[length(selected_item)]]
    expect_equal(selected_item$item, t1@item_list[[1]])
    expect_equal(selected_item$testlet, t1)
  }

  # -------------------------------------------------------------------------- #
  ### Sympson-Hetter ###
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2"), D = 1.702),
                id = "t1", misc = list(sympson_hetter_k = 0))
  t2 <- testlet(itempool(a = c(0.2, 0.2), b = 4:5, id = c("t2-i1", "t2-i2"),
                          D = 1.702), id = "t2",
                misc = list(sympson_hetter_k = 1))
  i1 <- item(b = -1, D = 1.702, id = "i1", misc = list(sympson_hetter_k = 1))
  i2 <- item(b = 0, D = 1.702, id = "i2", misc = list(sympson_hetter_k = 1))
  i3 <- item(b = 1, D = 1.702, id = "i3", misc = list(sympson_hetter_k = 1))
  ip <- c(t1, t2, i1, i2, i3)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "mepv",
                          next_item_par = list(var_calc_method = "eap"),
                          exposure_control_rule = "sympson-hetter",
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = 4))
  # Select a testlet as first item
  est_history <- list(
    list(est_before = -4, se_before = 0.3, resp = NA,
         item = NULL, testlet = NULL, est_after = NA, se_after = NA)
  )
  selected_item <- irt:::select_next_item_cpp(
    cd = cd, est_history = est_history, list())
  # selected_item <- irt:::select_next_item_cpp(
  #   cd = cd, est_history = est_history, additional_args = list())
  expect_true("set_aside_item_list" %in% names(selected_item$additional_args))
  expect_equal(selected_item$additional_args$set_aside_item_list[[1]], t1)
  selected_item <- selected_item$est_history
  selected_item <- selected_item[[length(selected_item)]]
  expect_equal(selected_item$item, i1)
  expect_null(selected_item$testlet)

  # skip("")
  # # -------------------------------------------------------------------------- #
  # # Performance Check of MEPV for large item pool
  # {
  #   cat("---------------------------\nStart time: ", strftime(Sys.time()), "\n")
  #   print(Sys.time())
  #   set.seed(123)
  #   n_items <- 500
  #   n_testlet <- 20
  #   ip <- itempool(c(sapply(paste0("testlet-", 1:n_testlet),
  #                            function(x) generate_testlet(id = x, n = 2)),
  #                     generate_ip(n = n_items, output = "list")))
  #   cd <- create_cat_design(ip = ip,
  #                           next_item_rule = "mepv",
  #                           next_item_par = list(var_calc_method = "owen"),
  #                           termination_rule = 'min_item',
  #                           termination_par = list('min_item' = 60))
  #   est_history <- create_est_history(num_of_steps = 50, ip = ip,
  #                                     cat_design = cd)$est_history
  #   # print_est_history(est_history)
  #   print(microbenchmark::microbenchmark(
  #     snic = irt:::select_next_item_cpp(
  #       cd = cd, est_history = est_history, additional_args = list()),
  #     bare = irt:::select_next_item_mepv_cpp(
  #       cd = cd, est_history = est_history, additional_args = list()),
  #     times = 100
  #   ))
  #   cat("---------------------------\nEnd time: ", strftime(Sys.time()), "\n")
  # }
  # # --------------------------- #
  # # Start time:  2020-07-06 15:55:36
  # # [1] "2020-07-06 15:55:36 EDT"
  # # Unit: milliseconds
  # # expr     min      lq     mean   median      uq      max neval
  # # snic 69.2849 70.2792 75.98840 71.19615 72.0817 120.5621    10
  # # bare 68.7588 69.9863 71.06902 70.62645 71.1626  77.8612    10
  # # --------------------------- #
  # # End time:  2020-07-06 15:55:42
  #
  # # --------------------------- #
  # # Start time:  2020-07-06 16:14:42
  # # [1] "2020-07-06 16:14:42 EDT"
  # # Unit: milliseconds
  # #  expr     min       lq     mean   median       uq       max neval
  # #  snic 58.6722 60.98655 63.81935 62.43820 64.10895   90.8495   100
  # #  bare 58.6808 61.84180 78.02609 63.00905 66.07295 1395.0916   100
  # # --------------------------- #
  # # End time:  2020-07-06 16:15:00
  #
  #
  #   cd_eap <- create_cat_design(ip = ip,
  #                           next_item_rule = "mepv",
  #                           next_item_par = list(var_calc_method = "eap"),
  #                           termination_rule = 'min_item',
  #                           termination_par = list('min_item' = 60))
  #   cd_owen <- create_cat_design(ip = ip,
  #                           next_item_rule = "mepv",
  #                           next_item_par = list(var_calc_method = "owen"),
  #                           termination_rule = 'min_item',
  #                           termination_par = list('min_item' = 60))
  #   est_history_eap <- create_est_history(num_of_steps = 50, ip = ip,
  #                                     cat_design = cd_eap)$est_history
  #   est_history_eap <- create_est_history(num_of_steps = 50, ip = ip,
  #                                     cat_design = cd_owen)$est_history
  #
  #   # print_est_history(est_history)
  #   print(microbenchmark::microbenchmark(
  #     eap = irt:::select_next_item_cpp(
  #       cd = cd_eap, est_history = est_history, additional_args = list()),
  #     owen = irt:::select_next_item_cpp(
  #       cd = cd_owen, est_history = est_history, additional_args = list()),
  #     times = 10
  #   ))
  #
  #
  #
  #
  #
  # # First performance check without any improvement.
  # # --------------------------- #
  # # Start time: "2020-07-06 13:49:11 EDT"
  # # Unit: seconds
  # # expr     min       lq     mean   median       uq      max neval
  # # snic 3.07022 3.074536 3.114801 3.105701 3.156388 3.190121    10
  # # bare 3.09067 3.104146 3.126132 3.127472 3.140221 3.173316    10
  # # End time: "2020-07-06 13:51:48 EDT"
  # # --------------------------- #
  #
  #
  # selected_item <- irt:::select_next_item_cpp(
  #   cd = cd, est_history = est_history, additional_args = list())
  # print_est_history(selected_item$est_history)

})

############################################################################@###
################### est_ability_cat_cpp ########################################
############################################################################@###
test_that("Test est_ability_cat_cpp Function", {
  # -------------------------------------------------------------------------- #
  n_ip <- 3
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          ability_est_rule = "eap",
                          ability_est_par = list(prior_dist = "norm",
                                                 prior_par = c(0, 1),
                                                 min_theta = -5,
                                                 max_theta = 5,
                                                 no_of_quadrature = 50
                                                 ),
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = n_ip))
  resp <- sample(0:1, n_ip, TRUE)
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = resp[1],
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[2],
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[3],
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2)
    )

  expected <- est_ability(ip, resp, method = "eap")
  observed <- irt:::est_ability_cat_cpp(true_ability = rnorm(1), cd = cd,
                                        est_history = est_history,
                                        additional_args = list())
  observed <- observed$est_history
  observed <- observed[[length(observed)]]
  expect_equal(expected$est, observed$est_after,  tol = 1e-5)
  expect_equal(expected$se, observed$se_after, tol = 1e-5)


  # -------------------------------------------------------------------------- #
  ### Owen's Bayesian Estimation ###
  n_ip <- 3
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0,
                                                 prior_var = 4),
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = n_ip))
  resp <- sample(0:1, n_ip, TRUE)
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = resp[1],
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[2],
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[3],
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2)
    )
  expected <- est_ability(ip, resp, method = "owen", prior_pars = c(0, 2))
  observed <- irt:::est_ability_cat_cpp(true_ability = rnorm(1), cd = cd,
                                        est_history = est_history,
                                        additional_args = list())
  observed <- observed$est_history
  observed <- observed[[length(observed)]]
  expect_equal(expected$est, observed$est_after,  tol = 1e-5)
  expect_equal(expected$se, observed$se_after, tol = 1e-5)

  # -------------------------------------------------------------------------- #
  # Check whether the ability estimation for Owen with all items is the same
  # as with informative prior.
  n <- 5 # number of items
  ip <- itempool(a = rlnorm(n, 0, .3), b = rnorm(n), c = runif(n, 0, .3))
  item <- item(a = rlnorm(1, 0, .3), b = rnorm(1), c = runif(1, 0, .3))
  resp <- sample(0:1, n, TRUE)
  resp_item <- sample(0:1, 1, TRUE)
  # Check for Owen
  expected <- irt:::est_ability_owen_cpp(
    ip = c(ip, item), resp = c(resp, resp_item), m0 = 0, v0 = 1)
  priors <- irt:::est_ability_owen_cpp(ip = ip, resp = resp, m0 = 0, v0 = 1)
  observed <- irt:::est_ability_owen_cpp(
    ip = itempool(item), resp = resp_item, m0 = priors$est, v0 = priors$se^2)

  expect_equal(expected$est, observed$est)
  expect_equal(expected$se, observed$se)

  # -------------------------------------------------------------------------- #
  # "sum_score" ability estimation works
  ip <- generate_ip(n = 30)
  cd <- create_cat_design(ip = ip,
                          next_item_rule = "fixed",
                          ability_est_rule = "sum_score",
                          termination_rule = 'max_item',
                          termination_par = list(max_item = 30))
  resp <- sample(0:1, 4, TRUE)
  est_history <- list(
    list(est_before = 0.2, se_before = 0.3, resp = resp[1],
         item = ip[[1]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[2],
         item = ip[[2]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[3],
         item = ip[[3]], testlet = NULL, est_after = 0.3, se_after = .2),
    list(est_before = 0.2, se_before = 0.3, resp = resp[4],
         item = ip[[4]], testlet = NULL, est_after = NA, se_after = 12)
    )
  observed <- irt:::est_ability_cat_cpp(true_ability = rnorm(1), cd = cd,
                                        est_history = est_history,
                                        additional_args = list())

  expect_equal(observed$est_history[[4]]$est_after, sum(resp))
  expect_true(is.na(observed$est_history[[4]]$se_after))

  # -------------------------------------------------------------------------- #
  ### Maximum-Likelihood Estimation using Newton-Raphson ###


  # # Benchmarking tests
  # foo1 <- function() {
  #   for (i in 1:20)
  #     irt:::est_ability_owen_cpp(
  #       ip = c(ip, item), resp = c(resp, resp_item), m0 = 0, v0 = 1)
  # }
  # foo2 <- function() {
  #   est1 <- irt:::est_ability_owen_cpp(ip = ip, resp = resp, m0 = 0, v0 = 1)
  #   for (i in 1:20)
  #     irt:::est_ability_owen_cpp(
  #       ip = itempool(item), resp = resp_item, m0 = est1$est, v0 = est1$se^2)
  # }
  # microbenchmark::microbenchmark(long = foo1(),  short = foo2())


  # # Check the same for EAP  -- NOT WORKING --
  # expected <- irt:::est_ability_eap_cpp(
  #   ip = c(ip, item), resp = matrix(c(resp, resp_item), nrow = 1),
  #   theta_range = c(-4, 4),
  #   no_of_quadrature = 50, prior_dist = "norm", prior_par = c(0, 1))
  # priors <- irt:::est_ability_eap_cpp(
  #   ip = ip, resp = matrix(resp, nrow = 1),
  #   theta_range = c(-4, 4),
  #   no_of_quadrature = 50, prior_dist = "norm", prior_par = c(0, 1))
  # observed <- irt:::est_ability_eap_cpp(
  #   ip = itempool(item), resp = matrix(resp_item, nrow = 1),
  #   theta_range = c(-4, 4), no_of_quadrature = 50, prior_dist = "norm",
  #   prior_par = c(priors$est, priors$se^2))
  #
  # expect_false(expected$est == observed$est)
  # expect_false(expected$se == observed$se)


})


############################################################################@###
################### est_ability_4pm_nr #########################################
############################################################################@###
test_that("Test est_ability_4pm_nr Function", {
  # -------------------------------------------------------------------------- #
  # Check whether the ML Newton-Raphson method is equal to ML estimate of
  for (i in 1:10) {
    ip <- generate_ip(
      model = sample(c("1PL", "2PL", "3PL", "4PL"), 1),
      n = sample(5:30, 1))
    resp <- sim_resp(ip = ip, theta = rnorm(1))

    mlnr <- irt:::est_ability_4pm_nr_itempool_cpp(
      resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 0.0001)
    mlr <- est_ability(ip = ip, resp = resp, method = "ml")
    expect_equivalent(mlnr, mlr$est, tol = 0.001)
  }
  # microbenchmark::microbenchmark(
  # cpp = irt:::est_ability_4pm_nr_itempool_cpp(resp = resp, ip = ip,
  #                                        theta_range = c(-5, 5),
  #                                        initial_est = 0,
  #                                        criterion = 0.01),
  # R = est_ability(ip = ip, resp = resp, method = "ml"))

  # -------------------------------------------------------------------------- #
  # Estimation with missing data
  for (i in 1:10) {
    ip <- generate_ip(
      model = sample(c("1PL", "2PL", "3PL", "4PL"), 1),
      n = sample(5:30, 1))
    resp <- sim_resp(ip = ip, theta = rnorm(1), prop_missing = runif(1, .1, .4))

    mlnr <- irt:::est_ability_4pm_nr_itempool_cpp(
      resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 0.0001)
    mlr <- est_ability(ip = ip, resp = resp, method = "ml")
    expect_equivalent(mlnr, mlr$est, tol = 0.001)
  }

  ip <- generate_ip(model = "3PL", n = sample(10:20, 1))

  # -------------------------------------------------------------------------- #
  # All correct and all incorrect answers.
  for (i in 0:1) {
    ip <- generate_ip(
      model = sample(c("1PL", "2PL", "3PL", "4PL"), 1),
      n = sample(5:30, 1))
    resp <- rep(i, length(ip))
    theta_range <- c(runif(1, -6, -2), runif(1, 2, 6))
    est <- irt:::est_ability_4pm_nr_itempool_cpp(
      resp = resp, ip = ip, theta_range = theta_range, criterion = 0.0001)
    expect_equal(est, theta_range[i + 1])
  }

  # -------------------------------------------------------------------------- #
  # There should not be a difference between Newton-Raphson and R's root finder
  # through optimization.
  ip <- itempool(a = c(0.5707, 0.591, 1.7997, 1.9157, 0.8327, 1.4858, 0.8798,
                       1.259, 0.6271, 0.7729, 0.7145),
                 b = c(-0.5246, -1.4423, 1.7301, -0.4387, 1.1685, 0.1699,
                       -0.5488, 0.2341, 0.6754, -0.9755, -0.8838),
                 c = c(0.0863, 0.2639, 0.2258, 0.1, 0.2839, 0.2664, 0.1728,
                       0.0304, 0.1962, 0.2932, 0.2882))
  resp <- c(0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0)
  est_opt <- est_ability(ip = ip, resp = resp, method = "ml", )$est
  est_nr <- irt:::est_ability_4pm_nr_itempool_cpp(
    resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06)
  expect_equivalent(est_opt, est_nr, tol = 1e-3)


  # Another Example
  ip <- itempool(a = c(0.712387267381568, 1.42924280545444, 0.533937500649207,
                       0.908658069797544, 1.01996576369619, 1.58730591908994,
                       0.890809537098302, 0.849616771350515, 1.33348447894585,
                       3.00943700818358, 0.856010859069746),
                 b = c(0.70438285051974, 1.46974377070636, -0.704406890404017,
                       0.0421789052992931, -0.0933942390642164, 1.5489756793278,
                       0.758683589460407, -3.40374425923407, 0.543804063596402,
                       0.397789585537004, 0.291197798755969),
                 c = c(0.225277058640495, 0.266812896472402, 0.17794270128943,
                       0.195420338166878, 0.0434793205233291, 0.0644042973406613,
                       0.00368709010072052, 0.131263936311007, 0.0388918435899541,
                       0.165485794586129, 0.264024762692861))
  resp <- c(1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 0)
  est_opt <- est_ability(ip = ip, resp = resp, method = "ml", )$est
  est_nr <- irt:::est_ability_4pm_nr_itempool_cpp(
    resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06,
    initial_estimates = c(-5 + 1e-06, 0, 5 - 1e-06))
  expect_equivalent(est_opt, est_nr, tol = 1e-3)

  # # DO NOT DELETE FOLLOWING AND UNCOMMENT IT.
  # # Another Example in favor of Newton-Raphson
  # ip <- itempool(a = c(0.501, 1.0089, 1.193, 1.0041, 1.2288, 1.0915, 1.1783,
  #                      0.8865, 0.929, 0.7253, 1.0934, 1.2778, 0.7535, 0.8765,
  #                      0.849),
  #                b = c(0.6729, 1.1055, 1.3649, -0.2446, -0.5923, -1.1431,
  #                      0.8818, -0.8183, -1.1589, 0.96, -0.8499, 1.069, -0.5085,
  #                      -0.3019, 0.4582),
  #                c = c(0.1776, 0.0553, 0.1328, 0.2787, 0.1769, 0.2688, 0.2688,
  #                      0.24, 0.0501, 0.0938, 0.1258, 0.2398, 0.2416, 0.1677,
  #                      0.2291))
  # resp <- c(0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0)
  # est_opt <- est_ability(ip = ip, resp = resp, method = "ml", )$est
  # est_nr <- irt:::est_ability_4pm_nr_itempool_cpp(
  #   resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06,
  #   initial_estimates = c(-5 + 1e-06, 0, 5 - 1e-06))
  # expect_equivalent(est_opt, est_nr, tol = 1e-3)
  #
  # # Another Example with an interesting and clear shape
  # ip <- itempool(a = c(1.7307, 1.0399, 0.7848, 0.6746),
  #                b = c(0.2322, 1.2937, -1.8468, 1.2453),
  #                c = c(4e-04, 0.1043, 0.1146, 0.1193),
  #                d = c(0.906, 0.903, 0.9741, 0.9037))
  # resp <- c(0, 1, 1, 1)
  # est_opt <- est_ability(ip = ip, resp = resp, method = "ml", )$est
  # est_nr <- irt:::est_ability_4pm_nr_itempool_cpp(
  #   resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06,
  #   initial_estimates = c(-5 + 1e-06, 0, 5 - 1e-06))
  # expect_equivalent(est_opt, est_nr, tol = 1e-3)
  # plot_resp_loglik(ip = ip, resp = resp)

  # # Another Example
  # ip <- itempool(a = c(),
  #                b = c(),
  #                c = c())
  # resp <- c()

  # writeClipboard(paste0(ip$a, collapse = ", "))
  # writeClipboard(paste0(ip$b, collapse = ", "))
  # writeClipboard(paste0(ip$c, collapse = ", "))
  # writeClipboard(paste0(ip$d, collapse = ", "))
  # writeClipboard(paste0(resp, collapse = ", "))

  # # ## ## ##
  # # # Generate item response strings and item pools that have discrepancy between
  # # # Maximum likelihood ability estimation based on R's optimization and
  # # # Newton-Raphson algorithm.
  # counter <- 0
  # while (counter < 1e7) {
  #   counter <- counter + 1
  #   ip <- generate_ip(model = "4PL", n = sample(2:4, 1))
  #   resp <- sim_resp(ip = ip, theta = rnorm(1))
  #   est_opt <- est_ability(ip = ip, resp = resp, method = "ml", )$est
  #   est_nr <- irt:::est_ability_4pm_nr_itempool_cpp(
  #     resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06)
  #   # resp_loglik(ip = ip, resp = resp, theta = est_opt)
  #   # resp_loglik(ip = ip, resp = resp, theta = est_nr)
  #   if (counter %% 1000 == 0) cat("Finished", counter, "iterations.\n")
  #   if (abs(est_opt - est_nr) > 0.01) {
  #     cat("Non-equality found!!")
  #     plot_resp_loglik(ip = ip, resp = resp)
  #     break
  #   }
  # }
  # # ## ## ##
  # ## Benchmarking between est_ability vs est_ability_nr
  # ip <- generate_ip(model = "4PL", n = sample(100:200, 1))
  # resp <- sim_resp(ip = ip, theta = rnorm(1))
  # microbenchmark::microbenchmark(
  #   est_opt = est_ability(ip = ip, resp = resp, method = "ml"),
  #   est_nr = irt:::est_ability_4pm_nr_itempool_cpp(
  #     resp = resp, ip = ip, theta_range = c(-5, 5), criterion = 1e-06)
  #   )


  # -------------------------------------------------------------------------- #


})


############################################################################@###
################### terminate_test_cpp #########################################
############################################################################@###
test_that("terminate_test_cpp Function", {
  # -------------------------------------------------------------------------- #
  ### max_item ###
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(ip = ip,
                          termination_rule = "max_item",
                          termination_par = list("max_item" = 3))
  # Test should not end before the maximum number of items administered
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))
  # Test should end when the maximum number of items administered
  co <- create_est_history(num_of_steps = 3, ip = ip, empty_last_step = FALSE)
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))

  # -------------------------------------------------------------------------- #
  ### min_item ###
  # Test min_item rule
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(ip = ip,
                          termination_rule = "min_item",
                          termination_par = list("min_item" = 3))
  # Test should not end before the minimum number of items administered
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))

  # Test should not end even minimum number of items administered
  co <- create_est_history(num_of_steps = 5, ip = ip)
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))

  # -------------------------------------------------------------------------- #
  ### min_se ###
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(
    ip = ip, termination_rule = "min_se",
    termination_par = list("min_se" = 0.75))
  # Test should not end before the minimum SE reached
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))

  # Test should end after minimum SE reached
  co <- create_est_history(num_of_steps = 12, ip = ip, empty_last_step = FALSE)
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))

  # -------------------------------------------------------------------------- #
  ### SPRT ###
  # Example on Parshall, Spray, et al. (2002) Practical considerations in
  # computer-based testing.
  ip <- itempool(
    a = c(1.80, 2.00, 1.80, 1.88, 1.71, .88, .63, .68, .92, .59, .47),
    b = c(-.64, -.16, -.83, .16, -1.51, -1.70, -1.22, -1.66, -2.09, -1.74,
          -1.94),
    c = c(.12, .15, .23, .19, .20, .31, .24, .22, .14, .18, .15),
    D = 1.7
    )
  resp <- c(1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)
  theta_0 <- -.75
  theta_1 <- -.05
  cut_score <- -0.4
  true_theta <- 1
  alpha <- 0.1
  beta <- 0.1
  A <- (1 - beta) / alpha
  B <- beta/(1-alpha)
  # -- #
  # Table 9.2  replicate:
  x <- cbind(data.frame(ip)[, -4],
        info = t(info(ip, theta = cut_score)),
        data.frame(P = t(prob(ip, c(theta_0, theta_1, true_theta)))),
        resp,
        llr_correct = sapply(1:length(ip), function(i)
          log(resp_lik(ip = ip[[i]], resp = 1, theta = theta_1)/
                resp_lik(ip = ip[[i]], resp = 1, theta = theta_0))),
        llr_incorrect = sapply(1:length(ip), function(i)
          log(resp_lik(ip = ip[[i]], resp = 0, theta = theta_1)/
                resp_lik(ip = ip[[i]], resp = 0, theta = theta_0))),
        total_ratio = sapply(1:length(ip), function(i)
          log(resp_lik(ip = ip[1:i], resp = resp[1:i], theta = theta_1)/
                resp_lik(ip = ip[1:i], resp = resp[1:i], theta = theta_0))),
        total_ratio_log_lik = sapply(1:length(ip), function(i)
          resp_loglik(ip = ip[1:i], resp = resp[1:i], theta = theta_1) -
                resp_loglik(ip = ip[1:i], resp = resp[1:i], theta = theta_0))
        )

  # This cat design is just for creating est_history.
  cd <- create_cat_design(
    ip = ip, termination_rule = c("min_item", "max_item"),
    next_item_rule = "fixed",
    next_item_par = list(item_id = x$id),
    termination_par = list("min_item" = nrow(x), "max_item" = nrow(x)))

  co <- create_est_history(num_of_steps = nrow(x), ip = ip, cat_design = cd,
                           empty_last_step = FALSE)
  for (i in 1:length(co$est_history))
    co$est_history[[i]]$resp <- x$resp[i]


  # Add a dummy item to item bank to disable max_item overrides sprt.
  ip  <- c(ip, item(b = 1))
  cd <- create_cat_design(
    ip = ip,
    next_item_rule = "fixed",
    next_item_par = list(item_id = ip$id),
    termination_rule = c("max_item", "sprt"),
    termination_par = list("sprt" = list(theta_0 = theta_0,
                                         theta_1 = theta_1,
                                         alpha = alpha,
                                         beta = beta),
                           "max_item" = length(ip)))
  # After 10 items test should not end
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history[1:10],
                                       list()))
  # After 11 items test should end
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history[1:11],
                                      list()))

  # -------------------------------------------------------------------------- #
  ## Scenario: "min_item" and "max_item"
  # A test should
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(
    ip = ip, termination_rule = c("min_item", "max_item"),
    termination_par = list("min_item" = 3, "max_item" = 5))
  # Test should not end since minimum item was not satisfied
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))

  # Test should not end since maximum item was not satisfied
  co <- create_est_history(num_of_steps = 4, ip = ip, empty_last_step = FALSE)
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))

  # Test should end since both minimum and maximum item was satisfied
  co <- create_est_history(num_of_steps = 5, ip = ip, empty_last_step = FALSE)
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))

  # -------------------------------------------------------------------------- #
  ## Scenario: "min_se" and "min_item"
  # Variable length test with minimum number of items.
  # Test should administer at least 5 items before checking whether min_se
  # is smaller than 0.75.
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(
    ip = ip, termination_rule = c("min_item", "min_se"),
    termination_par = list("min_item" = 3, "min_se" = 0.75))
  # Test should not end since minimum item was not satisfied even though min_se
  # satisfied.
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  co$est_history[[2]]$se_after <- 0.5
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))
  # Test should not end since minimum item satisfied but min_se was not
  # satisfied.
  co <- create_est_history(num_of_steps = 5, ip = ip, empty_last_step = FALSE)
  co$est_history[[5]]$se_after <- 0.85
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))
  # Test should end since both min_item and min_se satisfied.
  co <- create_est_history(num_of_steps = 5, ip = ip, empty_last_step = FALSE)
  co$est_history[[5]]$se_after <- 0.5
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))


  # -------------------------------------------------------------------------- #
  ## Scenario: "min_se" and "max_item"
  # Variable length test with maximum number of items.
  # Test should administer items until either the min_se <= 0.75 condition
  # satisfied or max_item reached.
  ip <- itempool(b = seq(-3, 3, by = 0.5))
  cd <- create_cat_design(
    ip = ip, termination_rule = c("min_se", "max_item"),
    termination_par = list("min_se" = 0.75, "max_item" = 5))
  # Test should end since min_se was satisfied even though max_item was
  # satisfied.
  co <- create_est_history(num_of_steps = 2, ip = ip, empty_last_step = FALSE)
  co$est_history[[2]]$se_after <- 0.5
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))
  # Test should not end since neither min_se or max_item satisfied.
  co <- create_est_history(num_of_steps = 4, ip = ip, empty_last_step = FALSE)
  co$est_history[[4]]$se_after <- 0.85
  expect_false(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                       est_history = co$est_history, list()))
  # Test should end since both max_item satisfied even though min_se not
  # satisfied.
  co <- create_est_history(num_of_steps = 5, ip = ip, empty_last_step = FALSE)
  co$est_history[[4]]$se_after <- 0.85
  expect_true(irt:::terminate_cat_cpp(true_ability = rnorm(1), cd = cd,
                                      est_history = co$est_history, list()))


})


############################################################################@###
################### cat_sim_single_cpp #########################################
############################################################################@###
test_that("cat_sim_single_cpp Function", {
  # A random test with three items.
  n_ip <- 3
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = n_ip))
  cs <- irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd)
  expect_is(cs, "cat_output")
  expect_equal(n_ip, length(cs$est_history))


  # A random test with min_items.
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip, next_item_rule = 'random',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  cs <- irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd)
  expect_is(cs, "cat_output")
  expect_equal(n_ip, length(cs$est_history))


  # A random test with Owen's Bayesian estimation
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'random',
                          ability_est_rule = "owen",
                          ability_est_par = list(prior_mean = 0, prior_var = 1),
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  cs <- irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd)
  expect_is(cs, "cat_output")
  expect_equal(n_ip, length(cs$est_history))


  # A long random test
  n_ip <- 100
  ip <- itempool(b = rnorm(n_ip))
  cd <- create_cat_design(
    ip = ip, next_item_rule = 'random',
    ability_est_par = list(prior_dist = "norm", prior_par = c(0, 3),
                           min_theta = -5, max_theta = 5,
                           no_of_quadrature = 50),
    termination_rule = 'max_item', termination_par = list('max_item' = n_ip))
  true_theta <- rnorm(1)
  cs <- irt:::cat_sim_single_cpp(true_ability = true_theta, cd = cd)
  # summary(cs)
  expect_is(cs, "cat_output")
  expect_equal(n_ip, length(cs$est_history))

  # Maximum Fisher Information Item Selection with EAP:
  n_ip <- 5
  ip <- itempool(a = runif(n_ip, .5, 1.5), b = rnorm(n_ip))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mfi',
                          termination_rule = 'min_item',
                          termination_par = list('min_item' = n_ip))
  cs <- irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd)
  expect_is(cs, "cat_output")
  expect_equal(n_ip, length(cs$est_history))


  # skip("Long Simulations")
  # nsim <- 100
  # output <- vector('list', nsim)
  # true_theta <- rnorm(1)
  # for (i in 1:nsim)
  #   output[[i]] <- irt:::cat_sim_single_cpp(true_ability = true_theta,
  #                                           cd = cd)
  # expected <- summary(output)
  # expect_equal(true_theta, mean(expected$est_ability), tol = 0.3)
  #
  #
  #
  # ## Micro benchmarking ##
  # n_ip <- 50
  # n_testlet <- 150
  # ip_list <- list()
  # for (i in seq_len(n_testlet)) {
  #   temp <- sample(2:4, 1)
  #   ip_list <- c(ip_list, testlet(itempool(
  #     b = rnorm(temp), id = paste0("t", i, "-i", 1:temp)),
  #     id = paste0("t", i), misc = list(sympson_hetter_k = runif(1, .5, 1))))
  # }
  # # Add individual items
  # ip <- c(itempool(ip_list), itempool(
  #   a = rlnorm(n_ip, 0, .3), b = rnorm(n_ip), c = runif(n_ip, 0, .3),
  #   misc = sapply(round(runif(n_ip, 0.5, 1), 2),
  #                 function(x) list(sympson_hetter_k = x))
  #   ))
  # cd <- create_cat_design(ip = ip,
  #                         next_item_rule = 'mepv',
  #                         exposure_control_rule = "sympson-hetter",
  #                         ability_est_rule = "owen",
  #                         ability_est_par = list(prior_mean = 0,
  #                                                prior_var = 1),
  #                         termination_rule = 'max_item',
  #                         termination_par = list(max_item = 30))
  # # x <- irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd)
  # microbenchmark::microbenchmark(
  #   irt:::cat_sim_single_cpp(true_ability = rnorm(1), cd = cd), times = 10)
})



############################################################################@###
############################################################################@###
############################################################################@###
# Test Other Functions
############################################################################@###
############################################################################@###
############################################################################@###

############################################################################@###
################### get_itempool_size #########################################
############################################################################@###

test_that("cat_sim_single_cpp Function", {
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")))
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")))
  ip <- c(t1, t2, itempool(b = -1:4))
  observed <- irt:::get_itempool_size(ip)
  expect_equivalent(observed[1], 8)
  expect_equivalent(observed[2], 2)
  expect_equivalent(observed[3], 11)
})


############################################################################@###
################### calculate_exposure_rates ###################################
############################################################################@###

test_that("calculate_exposure_rates Function", {
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")))
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")))
  ip <- c(t1, t2, itempool(b = rnorm(18)))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = 10))
  co <- cat_sim(true_ability = rnorm(3), cd = cd)
  observed <- irt:::calculate_exposure_rates_cpp(item_ids = cd$ip$id,
                                                 cat_output_list = co)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))
})



############################################################################@###
################### calculate_overlap_rates ####################################
############################################################################@###

test_that("calculate_overlap_rates Function", {
  t1 <- testlet(itempool(b = -4:-2, id = c("t1-i1", "t1-i2", "t1-i3")))
  t2 <- testlet(itempool(b = 2:3, id = c("t2-i1", "t2-i2")))
  ip <- c(t1, t2, itempool(b = rnorm(18)))
  cd <- create_cat_design(ip = ip,
                          next_item_rule = 'mepv',
                          termination_rule = 'max_item',
                          termination_par = list('max_item' = 10))
  co <- cat_sim(true_ability = rnorm(3), cd = cd)
  observed <- irt:::calculate_overlap_rates_cpp(item_ids = cd$ip$id,
                                                cat_output_list = co)
  expect_equal(length(observed), length(ip))
  expect_true(all(sapply(observed, function(x) x >=0 & x <=1)))
})
