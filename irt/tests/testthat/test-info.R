

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

###############################################################################@
############################# info (Item) @#####################################
###############################################################################@

test_that("info - Item", {

  # ---------------------------------------------------------------------------#
  # All models runs using single theta
  expect_is(info(ip = generate_item(model = "1PL"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "2PL"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "3PL"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "4PL"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "GRM"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "PCM"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "GPCM"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "GPCM2"), theta = rnorm(1)),
            "numeric")
  expect_is(info(ip = generate_item(model = "Rasch"), theta = rnorm(1)),
            "numeric")

  # ---------------------------------------------------------------------------#
  # All models runs using multiple theta
  expect_is(info(ip = generate_item(model = "1PL"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "2PL"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "3PL"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "4PL"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "GRM"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "PCM"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "GPCM"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "GPCM2"),
                 theta = rnorm(sample(2:5, 1))), "numeric")
  expect_is(info(ip = generate_item(model = "Rasch"),
                 theta = rnorm(sample(2:5, 1))), "numeric")

  ##  Single theta
  # ---------------------------------------------------------------------------#
  theta = 1.386
  b = -1.252
  D = 1.702
  ip = item(b = b, D = D, model = '1PL', id = paste0('i', 1:1))
  expected <- 0.0317905243367295
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # ---------------------------------------------------------------------------#
  theta = 0.663
  b = -1.391
  D = 1.702
  ip = item(b = b, D = D, model = '1PL', id = paste0('i', 1:1))
  expected <- 0.0827409079413328
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = -1.764
  a = 1.603
  b = 1.863
  c = 0.277
  D = 1.702
  ip = item(a = a, b = b, c = c, D = D, model = '3PL', id = paste0('i', 1:1))
  expected <- 4.93335767520636e-08
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = 1.211
  a = 0.623
  b = -0.837
  c = 0.136
  d = 0.954
  D = 1.702
  ip = item(a = a, b = b, c = c, d = d, D = D, model = '4PL',
            id = paste0('i', 1:1))
  expected <- 0.0562367939887916
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # Multiple theta
  # ---------------------------------------------------------------------------#
  theta = c(1.753, -1.77, -0.45)
  a = 1.239
  b = 1.465
  D = 1.702
  ip = item(a = a, b = b, D = D, model = '2PL', id = paste0('i', 1:1))
  expected <- c(1.01520882483516, 0.00483507962344057, 0.0756952424600891)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  # Observed Information
  theta = 0.832
  b = -0.682
  D = 1.702
  ip = item(b = b, D = D, model = '1PL', ID = paste0('i', 1:1))
  resp = c(TRUE)
  expected <- 0.19018685
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  theta = 0.452
  ip = item(a = 1.714, b = 0.898, c = 0.116, D = 1.702, model = '3PL',
            ID = paste0('i', 1:1))
  resp = c(1)
  expected <- -0.35163524
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  expect_equivalent(object, expected, tolerance = 1e-6)
  expect_equivalent(object, -resp_loglik(ip = ip, resp = resp, theta = theta,
                                         derivative = 2), tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  theta = -1.101
  a = 1.683
  b = 1.808
  c = 0.078
  D = 1.702
  ip = item(a = a, b = b, c = c, D = D, model = '3PL', ID = paste0('i', 1:1))
  resp = c(0)
  expected <- 0.00197251
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  theta = 1.301
  a = 1.283
  b = -1.586
  c = 0.234
  d = 0.97
  D = 1.702
  ip = item(a = a, b = b, c = c, d = d, D = D, model = '4PL',
               ID = paste0('i', 1:1))
  resp = c(1)
  expected <- 0.00658582
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  expect_equivalent(object, expected, tolerance = 1e-6)
  expect_equivalent(object, -resp_loglik(ip = ip, resp = resp, theta = theta,
                                         derivative = 2), tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  theta = -1.438
  D = 1.702
  a = 1.22
  b = c(-0.184, -2.862, -1.469, 0.169, -0.157)
  ip = item(a = a, b = b, D = D, model = 'GPCM')
  resp = 0
  expected <- 2.11764253
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  # TODO: Check the following test:
  # expect_equivalent(object, expected, tolerance = 1e-6)
  expect_equivalent(object, -resp_loglik(ip = ip, resp = resp, theta = theta,
                                         derivative = 2), tolerance = 1e-6)

  ### Graded Response Model ###
  # -------------------------------------------------------------------------- #
  theta = -1.822
  a = 0.812
  b = c(-1.545, -1.503, 0.536)
  D = 1.702
  ip = item(a = a, b = b, D = D, model = 'GRM', id = paste0('i', 1:1))
  expected <- 0.474452438901097
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  # Multiple Thetas
  theta = c(1.618, 1.172)
  a = 1.174
  b = c(-0.009, 0.426, 0.634)
  D = 1.702
  ip = item(a = a, b = b, D = D, model = 'GRM', id = paste0('i', 1:1))
  expected <- c(0.432263612036406, 0.776788276793933)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ### Generalized Partial Credit Model ###
  # -------------------------------------------------------------------------- #
  theta = -1.572
  D = 1.702
  a = 0.859
  b = c(-2.402, -1.106, 1.56)
  ip = item(a = a, b = b, D = D, model = 'GPCM')
  expected <- 0.32534918
  object <- info(ip = ip, theta = theta)
  # TODO: double check the information of Generalized partial credit model
  # expect_equivalent(object, expected, tolerance = 1e-6)
  expect_true(inherits(object, 'numeric'))

})


###############################################################################@
############################# info (Itempool) @################################
###############################################################################@

test_that("info - Itempool", {

  # ---------------------------------------------------------------------------#
  # All models runs using single theta
  expect_is(info(ip = generate_ip(model = "1PL"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "2PL"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "3PL"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "4PL"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "GRM"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "PCM"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "GPCM"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "GPCM2"), theta = rnorm(1)),
            "matrix")
  expect_is(info(ip = generate_ip(model = "Rasch"), theta = rnorm(1)),
            "matrix")

  # ---------------------------------------------------------------------------#
  # All models runs using multiple theta
  expect_is(info(ip = generate_ip(model = "1PL"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "2PL"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "3PL"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "4PL"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "GRM"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "PCM"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "GPCM"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "GPCM2"),
                 theta = rnorm(sample(2:5, 1))), "matrix")
  expect_is(info(ip = generate_ip(model = "Rasch"),
                 theta = rnorm(sample(2:5, 1))), "matrix")


  ##  Single theta
  # ---------------------------------------------------------------------------#
  theta = 0.226
  b = c(-1.529, -1.38, -0.593, -0.448, -0.08)
  D = 1.702
  ip = itempool(data.frame(b = b), D = D, model = '1PL',
                    id = paste0('i', 1:5))
  expected <- c(0.132414725735083, 0.166003207528995, 0.461363556709579,
                0.529896093199218, 0.677229701666625)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = -0.017
  a = c(1.875, 1.537, 0.961, 1.073, 1.79)
  b = c(-1.8, -1.335, 0.013, 1.653, 1.974)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b), D = D, model = '2PL',
                    id = paste0('i', 1:5))
  expected <- c(0.0341866121473768, 0.204493870492437, 0.668412411266363,
                0.144011156870809, 0.0214434042065179)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = -1.575
  a = c(1.489, 1.235, 0.93, 1.555, 1.999)
  b = c(-1.163, -0.977, -0.79, -0.284, 0.704)
  c = c(0.118, 0.289, 0.218, 0.184, 0.243)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c), D = D, model = '3PL',
                    id = paste0('i', 1:5))
  expected <- c(0.816980470555487, 0.268723586337344, 0.194023741905176,
                0.0266167442914817, 6.62213613273249e-06)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = -1.479
  a = c(0.913, 1.291, 1.18, 1.551, 1.766)
  b = c(-1.54, -0.4, -0.4, 1.091, 1.253)
  c = c(0.289, 0.116, 0.269, 0.343, 0.293)
  d = c(0.983, 0.924, 0.937, 0.936, 0.908)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c, d), D = D, model = '4PL',
                    id = paste0('i', 1:5))
  expected <- c(0.319106179781748, 0.127585971221441, 0.0684054134021274,
                1.3844056224489e-05, 1.2139206400796e-06)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ## Multiple theta and multiple item
  # ---------------------------------------------------------------------------#
  theta = c(1.425, -0.831, 0.041)
  a = c(1.871, 1.421, 1.183, 0.59, 1.433)
  b = c(-1.28, -0.164, -0.12, -0.033, 0.856)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b), D = D, model = '2PL',
                    id = paste0('i', 1:5))
  expected <- matrix(c(
    0.0018404879628147, 0.120135319887308, 0.165581637094441, 0.153835359014424,
    0.950925240905621, 1.58021391753661, 0.810390054168032, 0.631056873500936,
    0.215592201376169, 0.0940630086562428, 0.146655871578628, 1.37602654672357,
    0.987344633236002, 0.251746679726459, 0.630400798654426),
                     nrow = length(theta), byrow = TRUE)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ## TIF - Single Theta
  # ---------------------------------------------------------------------------#
  theta = 0.855
  a = c(1.685, 1.209, 1.804, 0.913, 1.061)
  b = c(-1.986, -1.082, -0.968, -0.835, 1.378)
  c = c(0.065, 0.178, 0.132, 0.188, 0.111)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c), D = D, model = '3PL',
                    id = paste0('i', 1:5))
  expected <- 0.67083317449735
  object <- info(ip = ip, theta = theta, tif = TRUE)
  expect_equivalent(object, expected, tolerance = 1e-6)
  # ---------------------------------------------------------------------------#
  theta = -0.711
  a = c(0.676, 1.713, 0.798, 1.592, 1.961)
  b = c(-1.688, -0.736, -0.209, 1.415, 1.446)
  c = c(0.211, 0.227, 0.094, 0.167, 0.179)
  d = c(0.938, 0.929, 0.939, 0.983, 0.923)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c, d), D = D, model = '4PL',
                    id = paste0('i', 1:5))
  expected <- sum(c(0.131316413890116, 1.08007543820367, 0.278790686278441,
                    0.000340028742097405, 2.32962070134738e-05))
  object <- info(ip = ip, theta = theta, tif = TRUE)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ## TIF - Multiple Theta
  # ---------------------------------------------------------------------------#
  theta = c(-0.864, -0.373, 1.384)
  a = c(0.546, 0.926, 1.465, 1.605)
  b = c(-1.284, -0.078, 0.735, 1.738)
  c = c(0.294, 0.082, 0.318, 0.329)
  D = 1.702
  ip = itempool(data.frame(a = a, b = b, c), D = D, model = '3PL',
                    id = paste0('i', 1:4))
  expected <- c(0.436162812117237, 0.631040513977452, 1.31562844559115)
  object <- info(ip = ip, theta = theta, tif = TRUE)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ### Graded Response Model ###
  # -------------------------------------------------------------------------- #
  theta = 1.516
  D = 1.702
  a = 1.032
  b = c(-0.968, -0.516, 1.57)
  item1 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.259
  b = c(-0.783, 0.407, 0.457)
  item2 = item(a = a, b = b, D = D, model = 'GRM')
  ip = itempool(c(item1, item2), id = paste0('i', 1:2))
  expected <- c(0.791551907180495, 0.3904719458343)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  # -------------------------------------------------------------------------- #
  # Different number of categories
  theta = 0.721
  D = 1.702
  a = 1.059
  b = c(0.259, 1.026)
  item1 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.071
  b = c(-0.558, 0.576, 1.902)
  item2 = item(a = a, b = b, D = D, model = 'GRM')
  a = 1.152
  b = c(-0.215, 1.167, 1.201, 1.624)
  item3 = item(a = a, b = b, D = D, model = 'GRM')
  ip = itempool(c(item1, item2, item3), id = paste0('i', 1:3))
  expected <- c(0.960425549030691, 0.950834524804291, 1.03609957317723)
  object <- info(ip = ip, theta = theta)
  expect_equivalent(object, expected, tolerance = 1e-6)

  ### Mixture of Models ###
  # -------------------------------------------------------------------------- #
  # Function can deal with the information of mixture of item models
  # Here the item pool consist of two 2PL, three GRM and two GPCM items.
  D = 1.702
  ip = itempool(data.frame(a = c(0.983, 1.03), b = c(-0.581, 1.643)),
                   D = D, model = '2PL', id = paste0('irt', 1:2))
  item1 = item(a = 1.027, b = c(-2.439, -0.21, 0.693), D = D, model = 'GRM')
  item2 = item(a = 1.024, b = c(0.241, 0.257, 0.311), D = D, model = 'GRM')
  item3 = item(a = 0.896, b = c(-0.225, 0.458, 0.764), D = D, model = 'GRM')
  ip = c(ip, itempool(c(item1, item2, item3), id = paste0('grm', 1:3)))
  item1 = item(a = 0.927, b = c(2.086, 0.107), D = D, model = 'GPCM')
  item2 = item(a = 1.201, b = c(-0.349, 1.162), D = D, model = 'GPCM')
  ip = c(ip, itempool(c(item1, item2), id = paste0('gpcm', 1:2)))
  expected <- info(ip = ip, theta = rnorm(1), tif = TRUE)
  expect_is(expected, 'matrix')

  # Observed information
  # -------------------------------------------------------------------------- #
  theta = 1.243
  b = c(-0.286, 0.839, 1, 1.094, 1.998)
  D = 1.702
  ip = itempool(data.frame(b = b), D = D, model = '1PL',
                    ID = paste0('i', 1:5))
  resp = sim_resp(ip = ip, theta)
  expected <- c(0.18605487, 0.64491951, 0.69409371, 0.71268098, 0.49170446)
  object <- info(ip = ip, theta = theta, observed = TRUE, resp = resp)
  expect_equivalent(object, expected, tolerance = 1e-6)
  object <- info(ip = ip, theta = theta, tif = TRUE, observed = TRUE,
                 resp = resp)
  expect_equivalent(object, sum(expected))

  # -------------------------------------------------------------------------- #
  # Information of an item pool with testlets
  ip1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4),
                     c = rlnorm(4, -2, .4), id = paste0("t1-i", 1:4))
  t1 <- testlet(ip1)
  ip2 <- itempool(a = rlnorm(3, 0, .3), b = rnorm(3),
                     c = rlnorm(3, -2, .4), id = paste0("t2-i", 1:3))
  t2 <- testlet(ip2)
  ip3 <- itempool(a = rlnorm(10, 0, .3), b = rnorm(10),
                     c = rlnorm(10, -2, .4))
  ip <- c(t1, t2, ip3)
  theta <- rnorm(1)
  expect_is(info(ip, theta), "matrix")
  expect_is(info(ip, theta, tif = TRUE), "matrix")
  expect_equivalent(info(ip, theta, tif = TRUE),
                    info(ip1, theta, tif = TRUE) +
                      info(ip2, theta, tif = TRUE) +
                      info(ip3, theta, tif = TRUE) )


  # -------------------------------------------------------------------------- #
  # Information of item pool mixed with models and testlets
  t1 <- testlet(itempool(b = rnorm(2), id = c("t1-i1", "t1-i2")),
                   id = "t1")
  t2 <- testlet(itempool(a = rlnorm(3, 0, .3), b = rnorm(3),
                                id = c("t2-i1", "t2-i2", "t2-i3")), id = "t2")
  i1 <- item(b = rnorm(1), id = "i1")
  i2 <- item(a = rlnorm(1, 0, .3), b = rnorm(1), c = .2, id = "i2")
  i3 <- item(a = rlnorm(1, 0, .3), b = sort(runif(3)), id = "i3")
  ip <- c(t1, t2, i1, i2, i3)
  theta <- rnorm(1)
  resp <- sim_resp(ip = ip, theta = theta)[1, ]

  resp1 <- resp
  resp1["i3"] <- NA
  ip1 <- c(t1, t2, i1, i2)
  expected <-  info(theta = theta, ip = ip1, tif = TRUE, observed = FALSE,
                    resp = NULL)
  observed <-  info(theta = theta, ip = ip, tif = TRUE, observed = FALSE,
                    resp = resp1)
  expect_equivalent(observed, expected)

  resp1 <- resp
  resp1[c("i3", "t1-i1", "t1-i2")] <- NA
  ip1 <- c(t2, i1, i2)
  expected <-  irt:::info_itempool_bare_cpp(
    theta = theta, ip = ip1, tif = TRUE, observed = FALSE, resp = NULL)
  observed <-  irt:::info_itempool_bare_cpp(
    theta = theta, ip = ip, tif = TRUE, observed = FALSE, resp = resp1)
  expect_equivalent(observed, expected)

  resp1 <- resp
  resp1[c("i3", "t1-i1")] <- NA
  ip1 <- c(t1@item_list[[2]], t2, i1, i2)
  expected <-  irt:::info_itempool_bare_cpp(
    theta = theta, ip = ip1, tif = TRUE, observed = FALSE, resp = NULL)
  observed <-  irt:::info_itempool_bare_cpp(
    theta = theta, ip = ip, tif = TRUE, observed = FALSE, resp = resp1)
  expect_equivalent(observed, expected)
})


###############################################################################@
############################# info (Testlet) @##################################
###############################################################################@
test_that("info - REST", {
  theta <- rnorm(1)
  ip <- itempool(a = rlnorm(10, 0, .3), b = rnorm(10),
                     c = rlnorm(10, -2, .4))
  testlet <- testlet(ip)
  expect_equivalent(info(ip, theta, tif = TRUE),
                    info(testlet, theta, tif = TRUE))


  # Ability Estimation with the testlets
  t1 <- testlet(itempool(b = rnorm(2), id = c("t1-i1", "t1-i2")), id = "t1")
  t2 <- testlet(itempool(b = rnorm(3), id = c("t2-i1", "t2-i2", "t2-i3")),
                   id = "t2")
  i1 <- item(b = -1, id = "i1")
  i2 <- item(b = 0, id = "i2")
  i3 <- item(b = 1, id = "i3")
  theta <- rnorm(1)
  ip <- c(t1, t2, i1, i2, i3)
  observed <- info(ip = ip, theta = theta, tif = TRUE)
  ip <- c(i1, i2, i3, t1@item_list, t2@item_list)
  expected <- info(ip = ip, theta = theta, tif = TRUE)
  expect_equal(observed, expected)
})

###############################################################################@
############################# info (REST) @#####################################
###############################################################################@

test_that("info - REST", {

  expect_error(info("1", theta = 2), regexp = "Cannot convert object to an")

  # ---------------------------------------------------------------------------#
  theta = 1.83
  ipdf = data.frame(a = c(0.722, 1.975, 0.725, 0.89),
                    b = c(-1.303, -1.202, -0.535, 1.662),
                    c = c(0.206, 0.159, 0.15, 0.294),
                    D = 1.702)
  expected <- 0.411672740009343
  object <- info(ip = ipdf, theta = theta, tif = TRUE)
  expect_equivalent(object, expected, tolerance = 1e-6)
})


