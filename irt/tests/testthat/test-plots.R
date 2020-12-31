# library(testthat)

###############################################################################@
############################# plot.Item  #######################################
###############################################################################@

test_that("Test plot.Item", {
  expect_silent(p <- plot(x = item(b = 0.3, D = 1), suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  ip <- item(a = 1.2, b = 0.3, c = .2)
  expect_silent(p <- plot(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  expect_silent(p <- plot(item(a = 1.2, b = 0.3, c = .2, d = .89, D = 1),
                          suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  ### Plot Graded Response Model ###
  ip <- item(a = 0.902, b = c(-1.411, 0.385, 1.79))
  expect_silent(p <- plot(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  # De Ayala (2009), p. 221, Fig. 8.9
  ip <- item(a = 1.5, b = c(-1, 1), D = 1)
  expect_silent(p <- plot(ip, suppress_plot = TRUE))
  plot(ip, suppress_plot = TRUE) + geom_vline(xintercept = 0)
  expect_is(p, 'ggplot')
  # De Ayala (2009), p. 222, Fig. 8.10
  ip <- item(a = 1.5, b = c(-1, 1.4, 2), D = 1)
  expect_silent(p <- plot(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  # Check category names
  expect_silent(p <- plot(ip, category_names = c("Strongly Disagree", "Disagree",
                                                 "Agree","Strongly Agree"),
                          suppress_plot = TRUE))
  expect_is(p, 'ggplot')
  # Suppress category names legend:
  expect_silent(p <- plot(ip, category_names = NULL, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  # Plot 2PM
  ip <- item(a = 0.8, b = 1, model = "GRM")
  expect_silent(p <- plot(ip, category_names = c("Incorrect", "Correct"),
                          legend_title = "Response", suppress_plot = TRUE))
  expect_is(p, 'ggplot')
})


###############################################################################@
############################# plot.Itempool  ##################################
###############################################################################@
test_that("Test plot.Itempool", {
  n <- sample(10:50,1)
  ip <- itempool(a = runif(n, .5, 2), b = rnorm(n), c = runif(n, 0, .3), D = 1)
  expect_silent(p <- plot(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  expect_silent(p <- plot(ip, size = .25, alpha = 0.3, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  # Test Characteristic Curve
  expect_silent(p <- plot(ip, tcc = TRUE, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  ip <- itempool(item(a = 0.902, b = c(-1.411, 0.385, 1.79)))
  # plot(ip)

})


###############################################################################@
############################# plot_info  #######################################
###############################################################################@
test_that("Test plot_info", {
  ip <- item(a = 1, b = 0.81, D = 1)
  expect_silent(p <- plot_info(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  n <- sample(10:20,1)
  ip <- itempool(data.frame(a = runif(n, .5, 2), b = rnorm(n),
                               c = runif(n, 0, .3), D = 1))
  expect_silent(p <- plot_info(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  expect_silent(p <- plot_info(ip, tif = TRUE, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  item1 <- item(a = runif(1, .5, 2), b = sort(runif(3, -2, 2)), model = 'GRM')
  item2 <- item(a = runif(1, .5, 2), b = sort(runif(2, -2, 2)), model = 'GRM')
  item3 <- item(a = runif(1, .5, 2), b = runif(1, -2, 2))
  ip <- c(item1, item2, item3)
  expect_silent(p <- plot_info(ip, suppress_plot = TRUE))
  expect_is(p, 'ggplot')
  expect_silent(p <- plot_info(ip, tif = TRUE, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  info(ip, theta = 3)

})


###############################################################################@
############################# plot_resp_loglik  ################################
###############################################################################@
test_that("Test plot_resp_loglik", {
  ip <- generate_ip(model = "3PL", n = sample(10:50,1))
  resp <- sim_resp(ip = ip, theta = rnorm(1))
  expect_silent(p <- plot_resp_loglik(ip, resp, suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  expect_silent(p <- plot_resp_loglik(ip, resp, likelihood = TRUE,
                                      suppress_plot = TRUE))
  expect_is(p, 'ggplot')

  ip <- generate_ip(model = "4PL", n = sample(8:10,1))
  expect_silent(p <- plot_resp_loglik(ip, resp = sim_resp(ip, theta = rnorm(1)),
                                      suppress_plot = TRUE))
  expect_is(p, 'ggplot')
})



###############################################################################@
############################# plot_empirical_icc  ##############################
###############################################################################@
test_that("Test plot_empirical_icc", {
  ip <- generate_ip(model = "3PL", n = 20)
  theta <- rnorm(100)
  resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
  # item = 1
  # bins = 10
  # type = "eicc"
  # theta = NULL
  # title = "my title"
  # suppress_plot = TRUE
  # x_axis_scale = "theta"

  # -------------------------------------------------------------------------- #
  # The default
  expect_silent(p <- plot_empirical_icc(resp = resp, item = 1,
                                       suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # Change bin number
  expect_silent(p <- plot_empirical_icc(resp = resp, type = "eicc", item = 1,
                                       bins = 3, suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # Change x axis scale
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 1, x_axis_scale = "number", suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 1, x_axis_scale = "theta", ip = ip,
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 1, x_axis_scale = "percent", suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # item can be item ID
  expect_silent(p <- plot_empirical_icc(resp = resp, item = "Item-3",
                                       suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # type = "oep"
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 2, type = "oep", ip = ip, x_axis_scale = "theta",
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # When type = "oep", only x_axis_scale = "theta" is allowed.
  expect_warning(p <- plot_empirical_icc(
    resp = resp, item = 2, type = "oep", ip = ip, x_axis_scale = "percent",
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # Expect error when graph type requires Itempool
  expect_error(p <- plot_empirical_icc(resp = resp, item = "Item-3",
                                      type = "oep", suppress_plot = TRUE))
  expect_error(p <- plot_empirical_icc(resp = resp, item = "Item-3",
                                      type = "oept", suppress_plot = TRUE))
  # Length of theta and number of rows of resp should be equal
  expect_error(p<- plot_empirical_icc(
    resp = resp, item = 1, ip = ip, type = "oept", theta = 1:2,
    suppress_plot = TRUE))
  # When x_axis_scale is "theta", either theta or an item pool should be
  # provided.
  expect_error(p<- plot_empirical_icc(
    resp = resp, item = 1, x_axis_scale = "theta", suppress_plot = TRUE))

  # expect_silent(p <- plot_empirical_icc(resp = resp, item = "Item-3", ip = ip,
  #                                     type = "oept", suppress_plot = TRUE))
  # expect_is(p, 'ggplot')

  # -------------------------------------------------------------------------- #
  # Check multiline x-axis labels
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 1, x_axis_scale = "percent", n_dodge = 2,
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # Polytomous items
  ip <- generate_ip(model = "GRM", n = 20)
  theta <- rnorm(1000)
  resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)

  expect_silent(p <- plot_empirical_icc(resp = resp, item = 1,
                                       suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # type = "oep"
  expect_silent(p <- plot_empirical_icc(
    resp = resp, item = 2, type = "oep", ip = ip, x_axis_scale = "theta",
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

})


###############################################################################@
############################# plot_distractor_icc  #############################
###############################################################################@
test_that("Test plot_distractor_icc", {
  # Without an item pool or theta
  n_item <- 40 # sample(8:12, 1)
  n_theta <- 10000 # sample(100:200, 1)
  raw_resp <- matrix(sample(LETTERS[1:4], n_item * n_theta, replace = TRUE),
                     nrow = n_theta, ncol = n_item,
                     dimnames = list(paste0("Examinee-", 1:n_theta),
                                     paste0("Item-", 1:n_item)))
  key <- sample(LETTERS[1:4], n_item, replace = TRUE)
  # item = 3
  # bins = 10
  # suppress_plot = TRUE
  # x_axis_scale = "percent"
  # title = "Distractor Graph"
  expect_silent(p <- plot_distractor_icc(raw_resp, item = 3, key, bins = 10,
                                         suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  # -------------------------------------------------------------------------- #
  # Plot with ip and theta
  n_item <- 40 # sample(8:12, 1)
  n_theta <- 4000 # sample(100:200, 1)
  theta <- rnorm(n_theta)
  ip <- generate_ip(n = n_item)
  resp <- sim_resp(ip = ip, theta = theta, prop_missing = .1)
  key <- sample(LETTERS[1:4], n_item, replace = TRUE)
  key_matrix <- matrix(rep(key, n_theta), ncol = n_item, byrow = TRUE)
  raw_resp <- ifelse(
    is.na(resp), NA,
    ifelse(resp == 1, key_matrix,
           matrix(sapply(key_matrix,
                         function(x) sample(setdiff(LETTERS[1:4], x), 1)),
                  ncol = n_item, byrow = TRUE)
           ))
  expect_silent(p <- plot_distractor_icc(raw_resp, item = 3, key, bins = 10,
                                         theta = theta, suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

  expect_silent(p <- plot_distractor_icc(
    raw_resp, item = 3, key, bins = 10, theta = theta, x_axis_scale = "percent",
    suppress_plot = TRUE))
  expect_is(p, 'ggplot'); rm(p)

})
