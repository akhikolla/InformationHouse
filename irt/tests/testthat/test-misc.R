
# library(rbenchmark, testthat)


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% get_max_possible_total_score %%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

test_that("get_max_possible_total_score", {
  i1 <- sample(10:20, 1) # Number of standalone items
  i2 <- sample(3:8, 1) # number of categories of of the GRM item
  i3 <- sample(3:10, 1) # testlet size

  ip <- c(itempool(b = rnorm(i1)),
          item(a = 1, b = sort(rnorm(i2))))
  # -------------------------------------------------------------------------- #
  # No resp vector
  expect_equal(get_max_possible_total_score(ip), i1 + i2)

  # -------------------------------------------------------------------------- #
  # resp vector with one examinee
  resp <- sim_resp(ip = ip, theta = rnorm(1))
  expect_equal(get_max_possible_total_score(ip, resp = as.vector(resp)), i1 + i2)

  # -------------------------------------------------------------------------- #
  # resp vector with multiple examinees
  n_theta <- sample(10:30, 1)
  resp <- sim_resp(ip = ip, theta = rnorm(n_theta))
  expect_equivalent(get_max_possible_total_score(ip, resp = resp),
                    rep(i1 + i2, n_theta))

  # -------------------------------------------------------------------------- #
  # The size of the resp vector should be appropriate.
  expect_error(get_max_possible_total_score(ip, resp = as.vector(resp)))

  # -------------------------------------------------------------------------- #
  # item pool should be an Itempool object
  expect_error(get_max_possible_total_score(ip[[1]]))

  # -------------------------------------------------------------------------- #
  # Testlet should work
  ip <- c(
    itempool(b = rnorm(i1)),
    item(a = 1, b = sort(rnorm(i2))),
    testlet(itempool(b = rnorm(i3), id = paste0("t1-", 1:i3)))
    )
  expect_equal(get_max_possible_total_score(ip), i1 + i2 + i3)

  # -------------------------------------------------------------------------- #
  # NA works for 1pl items
  ip <- itempool(b = rnorm(i1))
  resp <- as.vector(sim_resp(ip, theta = rnorm(1, 1)))
  expect_equal(get_max_possible_total_score(ip), i1)
  expect_equal(get_max_possible_total_score(ip, resp), i1)
  resp[c(3, 5, 8)] <- NA
  expect_equal(get_max_possible_total_score(ip, resp), sum(!is.na(resp)))

  # -------------------------------------------------------------------------- #
  # NA works for 1pl items
  ip <- c(
    itempool(b = rnorm(i1)),
    item(a = 1, b = sort(rnorm(i2))),
    testlet(itempool(b = rnorm(i3), id = paste0("t1-", 1:i3)))
    )
  resp <- as.vector(sim_resp(ip, theta = rnorm(1, 3)))
  resp[c(3, 5, 8)] <- NA
  expect_equal(get_max_possible_total_score(ip, resp), i1 + i2 + i3 - 3)
  resp[i1 + 1] <- NA # set GRM item NA
  expect_equal(get_max_possible_total_score(ip, resp), i1 + i3 - 3)

  # -------------------------------------------------------------------------- #
  # Dichotomous and polytomous items
  ip <- generate_ip(model = c("3PL", "GRM", "3PL", "GRM", "GRM"),
                    n_categories = c(2, 5, 2, 4, 6))
  expect_equal(get_max_possible_total_score(ip), 1 + 4 + 1 + 3 + 5)

})




############################################################################@###
######################## Test convert_resp function ############################
############################################################################@###
test_that("Test convert_resp Function", {
  n_theta <- sample(13:24, 1)
  n_item <- sample(15:36, 1)
  resp <- data.frame(
    examinee_id = rep(paste0("Ex-", 1:n_theta), each = n_item),
    item_id = rep(paste0("Item-", 1:n_item), n_theta),
    score = sample(0:1, n_item * n_theta, replace = TRUE)
    )
  # Remove some of the columns
  resp_na <- resp
  resp_na$score[sample(1:nrow(resp), floor(0.9 * min(n_theta, n_item)))] <- NA
  # resp_wide <- convert_resp(resp_na[!is.na(resp_na$score), ])
  expect_true(TRUE)
  # for (i in sample(1:max(n_item, n_theta)))
  #   expect_equal(resp_na[i, "score"], resp_wide[resp_na[i, "examinee_id"],
  #                                               resp_na[i, "item_id"]])


  # # Check the speed. It is very slow compared to tidyr alternative.
  # # Check https://github.com/tidyverse/tidyr/blob/master/R/spread.R
  # n_theta <- 50000
  # n_item <- 100
  # resp <- data.frame(
  #   examinee_id = rep(paste0("Ex-", 1:n_theta), each = n_item),
  #   item_id = rep(paste0("Item-", 1:n_item), n_theta),
  #   score = sample(0:1, n_item * n_theta, replace = TRUE)
  #   )
  # resp_na <- resp
  # resp_na$score[sample(1:nrow(resp), floor(0.9 * min(n_theta, n_item)))] <- NA
  # resp_na <- resp_na[!is.na(resp_na$score), ]
  # microbenchmark::microbenchmark(
  # tidyr = tidyr::spread(resp_na, item_id, score) %>%
  #   tibble::column_to_rownames("examinee_id"),
  # convert_resp = convert_resp(resp_na),
  # times = 3L
  # )
  # profvis::profvis(convert_resp(resp_na))
})
