

test_that("rsss", {

  # -------------------------------------------------------------------------- #
  # Raw score argument cannot be larger than the maximum score of the test
  ip <- generate_ip()
  expect_error(rsss(ip = ip, raw_score = get_maximum_possible_score(ip) + 2))
  expect_error(rsss(ip = ip, raw_score = -1))
  # Error when a very small raw score presented
  expect_error(rsss(ip = ip, raw_score = 0),
               "Please provide a wider 'theta_range'")

  # -------------------------------------------------------------------------- #
  # A simple 3PL item pool
  ip <- generate_ip(model = "2PL")
  tol <- 0.0001
  theta_range <- c(-5, 5)
  n_score <- sample(2:6, 1)
  raw_score <- sort(sample(0:sum(ip$resp_max_score), n_score))
  # Remove minimum and maximum scores
  raw_score <- raw_score[!raw_score %in% c(0, ip$max_score)]
  scale_score <- rsss(ip = ip, raw_score = raw_score, theta_range = theta_range)
  # Select a random raw score which is not at the extremes
  i <- sample(which(scale_score < theta_range[2] &
                      scale_score > theta_range[1]), 1)
  expect_equivalent(raw_score[i], sum(prob(ip = ip, theta = scale_score[i],
                                           expected_value = TRUE)),
                    tolerance = tol)

  # -------------------------------------------------------------------------- #
  # Mixed dichotomous and polytomous item pool
  ip <- generate_ip(model = c(rep("2PL", 5), rep("GRM", 5)))
  tol <- 0.0001
  n_score <- sample(2:6, 1)
  theta_range <- c(-5, 5)
  raw_score <- sort(sample(0:sum(ip$resp_max_score), n_score))
  raw_score <- raw_score[!raw_score %in% c(0, ip$max_score)]
  scale_score <- rsss(ip = ip, raw_score = raw_score, theta_range = theta_range)
  # Select a random raw score which is not at the extremes
  i <- sample(which(scale_score < theta_range[2] &
                      scale_score > theta_range[1]), 1)
  expect_equivalent(raw_score[i], sum(prob(ip = ip, theta = scale_score[i],
                                           expected_value = TRUE)),
                    tolerance = tol)

  # -------------------------------------------------------------------------- #
  # Get raw scores from scale scores
  ip <- generate_ip(model = c(rep("2PL", 5), rep("GRM", 5)))
  tol <- 0.0001
  n_score <- sample(2:6, 1)
  theta_range <- c(-5, 5)
  scale_score <- sort(rnorm(n_score))
  raw_score <- rsss(ip = ip, scale_score = scale_score,
                    theta_range = theta_range)
  i <- sample(1:n_score, 1) # Select a random raw score
  p <- prob(ip = ip, theta = scale_score[i])
  expected <- sum(p * matrix(0:(ncol(p)-1), nrow = nrow(p), ncol = ncol(p),
                             byrow = TRUE), na.rm = TRUE)
  expect_equivalent(raw_score[i], expected, tolerance = tol)

  # -------------------------------------------------------------------------- #
  # Get scaled scores from raw scores - one theta
  ip <- generate_ip(model = c(rep("3PL", 5), rep("GRM", 5)))
  tol <- 0.001
  scale_score <- rnorm(1)
  p <- prob(ip = ip, theta = scale_score)
  expected <- sum(p * matrix(0:(ncol(p)-1), nrow = nrow(p), ncol = ncol(p),
                             byrow = TRUE), na.rm = TRUE)
  expect_equivalent(rsss(ip = ip, scale_score = scale_score,
                         theta_range = theta_range), expected, tolerance = tol)



})
