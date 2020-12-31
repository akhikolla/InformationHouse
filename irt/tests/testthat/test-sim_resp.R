
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% sim_resp %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# sim_resp (Item) @#################################
###############################################################################@

test_that("sim_resp - Item", {
  ##  Single theta
  expect_true(sim_resp(ip = item(b = rnorm(1)), theta = rnorm(1)) %in% c(0, 1))

  ip <- generate_item()
  expect_true(sim_resp(ip = ip, theta = rnorm(1)) %in% c(0, 1))

  # With D parameter
  ip <- generate_item(model = "2PL")
  expect_true(sim_resp(ip = ip, theta = rnorm(1)) %in% c(0, 1))

  theta <- rnorm(1)
  ip <- generate_item(model = "4PL")
  expect_true(sim_resp(ip = ip, theta = theta) %in% c(0, 1))

  ##  Multiple theta
  theta <- rnorm(5)
  ip <- generate_item(model = "4PL")
  expect_true(all(sim_resp(ip = ip, theta = theta) %in% c(0, 1)))


  # Test whether function generates approximately accurate responses
  compare <- function(v, tolerance = 0.001) all(abs(v[-1]- v[1])<tolerance)
  # Create item parameters

  # ---------------------------------------------------------------------------#
  N <- 100000 # number of examinees
  theta = rep(0.375, N)
  a = 0.904
  b = 0.178
  D = 1.702
  ip = item(a = a, b = b, D = D, model = '2PL')
  expected <- 0.575201573463713
  resp2 <- sim_resp(ip = ip, theta = theta)
  expect_equal(mean(resp2), expected, tolerance = 0.01)

  ### Graded Response Model
  # ---------------------------------------------------------------------------#
  theta <- rnorm(1)
  ip <- item(a = runif(1, .5, 2), b = c(-1, 0, 1), D = 1)
  expect_true(sim_resp(ip = ip, theta = theta) %in% 0L:3L)
  n <- 5
  theta <- rnorm(n)
  expect_true(all(sim_resp(ip = ip, theta = theta) %in% 0L:3L))
  expect_equal(length(sim_resp(ip = ip, theta = theta)), n)

  ### Generalized Partial Credit Model
  # ---------------------------------------------------------------------------#
  theta <- rnorm(1)
  ip <- item(a = runif(1, .5, 2), b = c(-1.2, 0, 1.2), D = 1, model = "GPCM")
  expect_true(sim_resp(ip = ip, theta = theta) %in% 0L:3L)
  n <- 5
  theta <- rnorm(n)
  expect_true(all(sim_resp(ip = ip, theta = theta) %in% 0L:3L))
  expect_equal(length(sim_resp(ip = ip, theta = theta)), n)

  ### A mixture of item types
  # ---------------------------------------------------------------------------#
  theta <- rnorm(1)
  ip <- c(item(a = runif(1, .5, 2), b = c(-2, -1, 0, 1, 2, 3), D = 1,
                  model = "GPCM"),
          item(a = runif(1, .5, 2), b = c(-1.1, 0, .92, 1.7), D = 1,
                  model = "GRM"),
          item(a = runif(1, .5, 2), b = runif(1, -2, 2), D = 1),
          item(a = runif(1, .5, 2), b = runif(1, -2, 2), c = .3, D = 1))
  resp <- sim_resp(ip = ip, theta = theta)
  expect_true(resp[1] %in% 0:6)
  expect_true(resp[2] %in% 0:4)
  expect_true(resp[3] %in% 0:1)
  expect_true(resp[4] %in% 0:1)

})


###############################################################################@
############################# sim_resp (Itempool) @############################
###############################################################################@

test_that("sim_resp - Itempool", {
  tol <- 1e-5 # tol

  theta <- rnorm(nTheta <- 14);     n <- 6;         D = runif(1, .5,3)
  names(theta) <- paste0("subject-", 1:nTheta)
  ip <- itempool(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3),
                 id = paste0("Item-",1:n), content = rep("Algebra", n))
  expect_true(all(  sim_resp(ip = ip, theta = theta) %in% c(0, 1)))

  theta <- rnorm(nTheta <- 1);     n <- 1  ;         D = runif(1, .5,3)
  names(theta) <- paste0("subject-", 1:nTheta)
  ipdf <- data.frame()
  ip <- item(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3),
             id = paste0("Item-",1:n), content = rep("Algebra", n))
  ip <- itempool(list(ip))
  expect_true(all(  sim_resp(ip = ip, theta = theta) %in% c(0, 1)))

  theta <- rnorm(nTheta <- 1);     n <- 8  ;         D = runif(1, .5,3)
  names(theta) <- paste0("subject-", 1:nTheta)
  ip <- itempool(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3),
                 id = paste0("Item-",1:n), content = rep("Algebra", n))
  expect_true(all(sim_resp(ip = ip, theta = theta) %in% c(0, 1)))

  theta <- rnorm(nTheta <- 5);     n <- 1  ;         D = runif(1, .5,3)
  names(theta) <- paste0("subject-", 1:nTheta)
  ip <- item(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3),
             id = paste0("Item-",1:n), content = rep("Algebra", n))
  ip <- itempool(list(ip))
  expect_true(all(sim_resp(ip = ip, theta = theta) %in% c(0, 1)))

  # ---------------------------------------------------------------------------#
  # response matrix should have the the names of the theta vector
  theta_names <- paste0("subject-", 1:5)
  theta <- rnorm(5)
  names(theta) <- theta_names
  ip <- itempool(b = rnorm(6))
  expect_equal(rownames(sim_resp(ip = ip, theta = theta)), theta_names)
  # For only one item
  expect_equal(rownames(sim_resp(ip = itempool(b = rnorm(6)),
                                 theta = theta[1])), theta_names[1])

  # ---------------------------------------------------------------------------#
  # response matrix should have default examinee names if theta vector does not
  # have a name.
  theta <- rnorm(5)
  expect_null(names(theta))
  expect_equal(rownames(sim_resp(ip = itempool(b = rnorm(6)), theta = theta)),
               paste0("S", 1:5))
  # For only one item
  expect_equal(rownames(sim_resp(ip = itempool(b = rnorm(6)), theta = theta[1])),
               paste0("S", 1))


  ### Graded Response Model
  # ---------------------------------------------------------------------------#
  theta <- rnorm(1)
  item1 <- item(a = runif(1, .5, 2), b = c(-1, 0, 1), D = 1)
  item2 <- item(a = runif(1, .5, 2), b = c(0, 2), D = 1)
  ip <- c(item1, item2)
  expect_true(all(sim_resp(ip = ip, theta = theta)[, 1] %in% 0L:3L))
  expect_true(all(sim_resp(ip = ip, theta = theta)[, 2] %in% 0L:2L))
  n <- 5
  theta <- rnorm(n)
  resp <- sim_resp(ip = ip, theta = theta)
  expect_equal(dim(resp)[1], n)
  expect_equal(dim(resp)[2], 2)
  expect_true(all(sim_resp(ip = ip, theta = theta)[, 1] %in% 0L:3L))
  expect_true(all(sim_resp(ip = ip, theta = theta)[, 2] %in% 0L:2L))

  # -------------------------------------------------------------------------- #
  # If examinee vector has names, then output response matrix will have
  # those names as row names.
  ip <- itempool(b = rnorm(10))
  row_names <- paste0("Ex-", 1:4)
  theta <- setNames(rnorm(4), row_names)
  resp <- sim_resp(ip = ip, theta = theta)
  expect_equal(rownames(resp), row_names)

  # Check for polytomous items
  ip <- itempool(data.frame(a = 1.2, b1 = rnorm(3), b2 = rnorm(3)),
                  model = "GRM")
  resp <- sim_resp(ip = ip, theta = theta)
  expect_equal(rownames(resp), row_names)

  # -------------------------------------------------------------------------- #
  # The item pool should be valid object.
  expect_error(sim_resp(ip = rnorm(5), theta = 1),
               "Cannot convert object to an 'Item' object.")

  # -------------------------------------------------------------------------- #
  # prop_missing should work.
  n_item <- sample(20:40, 1)
  n_theta <- sample(30:100, 1) * 10
  prop <- round(runif(1, .1, .94), 1)
  ip <- generate_ip(model = "GRM", n = n_item)
  theta <- rnorm(n_theta)
  resp <- sim_resp(ip = ip, theta = theta, prop_missing = prop)
  expect_equal(sum(is.na(resp)), prop * n_item * n_theta)
})

###############################################################################@
############################# sim_resp (Testlet) @##############################
###############################################################################@

test_that("sim_resp - REST", {
  # Ability Estimation with the testlets
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2")), id = "t1")
  t2 <- testlet(itempool(b = 2:4, id = c("t2-i1", "t2-i2", "t2-i3")),
                   id = "t2")
  i1 <- item(b = -1, id = "i1")
  i2 <- item(b = 0, id = "i2")
  i3 <- item(b = 1, id = "i3")
  ip <- c(t1, t2, i1, i2, i3)
  n_examinee <- 4
  resp <- sim_resp(ip = ip, theta = rnorm(n_examinee, 0, .4))
  expect_is(resp, "matrix")
  expect_equal(ncol(resp), length(t1) + length(t2) + 3)
  expect_equal(nrow(resp), n_examinee)

  # -------------------------------------------------------------------------- #
  # Collection of models
  testlet <- testlet(c(item(b = 1), item(a = .8, b = 3.1),
                     item(b = -1:1, model = "PCM")))
  expect_is(sim_resp(ip = testlet, theta = rnorm(1)), 'matrix')
})



###############################################################################@
############################# sim_resp (REST) @#################################
###############################################################################@

test_that("sim_resp - REST", {
  # Try character, 1 theta
  expect_error(sim_resp("1", theta = 2),
               "Cannot convert object to an 'Item' or an 'Itempool' object.")
  tol <- 1e-5 # tol

  # Try numeric, 1 theta
  theta <- rnorm(1);     n <- 1;        D = 1
  ipdf <- data.frame(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_true(all(  sim_resp(ip = ipdf, theta = theta)    %in% c(0, 1)))

  # Try numeric, multiple theta
  theta <- rnorm(4);     n <- 1;        D = 1
  ipdf <- data.frame(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_true(all(  sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))

  # Try matrix, 1 theta
  theta <- rnorm(1);     n <- 1;        D = 1
  ipdf <- matrix(c(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3)),
                 nrow = n, dimnames = list(NULL, c("a", "b", "c")))
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))

  # Try numeric, multiple theta
  theta <- rnorm(4);     n <- 5;        D = 1
  names(theta) <- paste0("subject-", 1:4)
  ipdf <- matrix(c(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3)),
                 nrow = n, dimnames = list(NULL, c("a", "b", "c")))
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))

  # Try data.frame, 1 theta
  theta <- rnorm(1);     n <- 1;        D = 1
  ipdf <- data.frame(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))

  # Try data.frame, multiple theta
  theta <- rnorm(7);     n <- 5;        D = 1
  names(theta) <- paste0("subject-", 1:7)
  ipdf <- data.frame(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))


  # Try data.frame, 1 theta
  theta <- rnorm(1);     n <- 1;        D = 1
  ipdf <- item(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))

  # Try data.frame, multiple theta
  theta <- rnorm(nTheta <- 14);     n <- 5;        D = 1
  names(theta) <- paste0("subject-", 1:nTheta)
  ipdf <- itempool(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3),
                   d = runif(n, .85, 1))
  ipdf <- ipdf$item_list
  expect_true(all(sim_resp(ip = ipdf, theta = theta)   %in% c(0, 1)))
})

