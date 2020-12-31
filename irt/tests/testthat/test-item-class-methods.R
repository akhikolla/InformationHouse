# library(testthat)

###############################################################################@
############################# item  ############################################
###############################################################################@
test_that("Test item", {
  # -------------------------------------------------------------------------- #
  ### model + parameters  ###
  # D parameter may be missing.
  expect_is(item(model = "2PL", parameters = list(a = 1.9, b = 1)), "Item")
  expect_is(item(model = "2PL", parameters = list(a = 1.9, b = 1, D = 1.7)),
            "Item")

  ### model + without parameters ###
  expect_is(item(model = "2PL", a = 1.9, b = 1, D = 1.7), "Item")
  expect_is(item(model = "2PL", a = 1.9, b = 1), "Item")

  # Rasch model
  item <- item(b = 1,  model = "Rasch")
  expect_is(item, "Item")
  expect_null(item$D)
  expect_null(item$a)

  # Create a GPCM model
  item <- item(parameters = list(a = 1.9, b = c(-1, 2), D = 1), model = "GPCM")
  expect_is(item, "Item")
  expect_equal(item$model, "GPCM")
  # GPCM without D parameter
  item <- item(parameters = list(a = 1.9, b = c(-1, 2)), model = "GPCM")
  expect_is(item, "Item")
  expect_equal(item$model, "GPCM")

  # Create a GPCM2 model
  item <- item(parameters = list(a = 1.9, b = 0.2, d = c(-1, 2), D = 1),
               model = "GPCM2")
  expect_is(item, "Item")
  expect_equal(item$model, "GPCM2")

  # Redundant parameters are ignored
  item <- item(a = 1.9, b = 1, c = .2,  model = "2PL")
  expect_is(item, "Item")
  expect_null(item$c)

  # Names of the previous parameters should be removed
  item <- item(a = 1.9, b = c(a = -1, b = 1.2, c = 2.2),  model = "GPCM")
  expect_null(names(item$b))

  # Additional parameters successfully added
  item <- item(a = 1.9, b = c(-1.1, 0.8, 1.6),  model = "GRM", id = "abc",
               content = "Geometry", misc = list(sympson_hetter_k = 1),
               se_parameters = list(a = .8, b = c(.2, .3, .4)))
  expect_is(item, "Item")
  expect_equal(item@id, "abc")
  expect_equal(item@content, "Geometry")
  expect_equal(item@misc, list(sympson_hetter_k = 1))
  expect_equal(item@se_parameters, list(a = .8, b = c(.2, .3, .4)))


  ### No model but with parameters ###
  # Create a 2PM model
  item <- item(parameters = list(a = 1.9, b = 1, c = .2))
  expect_is(item, "Item")
  expect_equal(item$model, "3PL")
  expect_equal(item$a, 1.9)
  expect_equal(item$b, 1)
  expect_equal(item$c, .2)
  # -------------------------------------------------------------------------- #
  # Create a GRM model
  item <- item(parameters = list(a = 1.9, b = c(-1, 2)),
               id = "i1", content = "Algebra")
  expect_is(item, "Item")
  expect_equal(item$model, "GRM")
  expect_equal(item$a, 1.9)
  expect_equal(item$b, c(-1, 2))
  expect_null(names(item$b))
  expect_equal(item$id, "i1")
  expect_equal(item$content, "Algebra")

  # -------------------------------------------------------------------------- #
  # Create a PCM model
  item <- item(parameters = list(b = c(x = -1, y = 2)))
  expect_is(item, "Item")
  expect_equal(item$model, "PCM")
  expect_null(names(item$b))
  expect_equal(item$b, c(-1, 2))


  ### No model or parameters ###
  # Create a 2PM model
  item <- item(a = 1.9, b = 1, c = .2)
  expect_is(item, "Item")
  expect_equal(item$model, "3PL")
  expect_equal(item$a, 1.9)
  expect_equal(item$b, 1)
  expect_equal(item$c, .2)
  # -------------------------------------------------------------------------- #
  # Unconventional vectors with only b parameters
  item <- item(b = list(1))
  expect_is(item, "Item")
  expect_equal(item$model, "Rasch")
  item <- item(b = matrix(1))
  expect_is(item, "Item")
  expect_equal(item$model, "Rasch")
  # -------------------------------------------------------------------------- #
  # Unnecessary parameters ("c") ignored
  item <- item(b = 1, c = .2)
  expect_is(item, "Item")
  expect_equal(item$model, "Rasch")
  expect_equal(item$b, 1)
  expect_null(item$c)
  # -------------------------------------------------------------------------- #
  # Create a GRM model
  item <- item(a = 1.9, b = c(-1, 2), id = "i1", content = "Algebra")
  expect_is(item, "Item")
  expect_equal(item$model, "GRM")
  expect_equal(item$a, 1.9)
  expect_equal(item$b, c(-1, 2))
  expect_equal(item$id, "i1")
  expect_equal(item$content, "Algebra")
  # -------------------------------------------------------------------------- #
  # Create a PCM model
  item <- item(b = c(-1, 2))
  expect_is(item, "Item")
  expect_equal(item$model, "PCM")
  expect_equal(item$b, c(-1, 2))
  # -------------------------------------------------------------------------- #
  # Create a M3PL model
  item <- item(a = c(1.9, 1.1), d = 1, c = .2)
  expect_is(item, "Item")
  expect_equal(item$model, "M3PL")
  # -------------------------------------------------------------------------- #
  # Create a M2PL model
  item <- item(a = c(1.9, 1.1), d = 1)
  expect_is(item, "Item")
  expect_equal(item$model, "M2PL")

  ### Other Cases ###
  # -------------------------------------------------------------------------- #
  # A row of a data frame can be used for creating items:
  dtf <- data.frame(a = c(1, 2), b1 = c(-2, 1), b2 = c(.2, .9))
  expect_is(item(a = dtf$a[1], b = dtf[1, -1], model = "GRM"), 'Item')
  expect_is(item(a = dtf$a[1], b = dtf[1, -1]), 'Item')

  # -------------------------------------------------------------------------- #
  # Can process an Item class object within item() function
  item <- item(a = 1.9, b = c(-1.1, 0.8, 1.6),  model = "GRM")
  expect_is(item(item), "Item")
  expect_null(item$id)
  expect_null(item$content)
  expect_null(item$misc)
  item <- item(item, id = "abc",
               content = "Geometry", misc = list(sympson_hetter_k = 1),
               se_parameters = list(a = .8, b = c(.2, .3, .4)))
  expect_equal(item$id, "abc")
  expect_equal(item$content, "Geometry")

  # -------------------------------------------------------------------------- #
  # item() function should accept "id" "ID", "iD"
  id <- "xyz123"
  item <- item(a = 1.9, b = 1, id = id)
  expect_equal(item$id, id)
  item <- item(a = 1.9, b = 1, Id = id)
  expect_equal(item$id, id)
  item <- item(a = 1.9, b = 1, iD = id)
  expect_equal(item$id, id)
  item <- item(a = 1.9, b = 1, ID = id)
  expect_equal(item$id, id)

  # -------------------------------------------------------------------------- #
  ##### Errors #####
  ## model ##
  # Model should be valid
  expect_error(item(model = "abc"), "Invalid 'model' specification.")

  # Model without parameters arguments but incomplete parameters
  expect_error(item(model = "2PL", a = 1, D = 1), "Incomplete parameters.")
  expect_error(item(model = "2PL", a = 1), "Incomplete parameters.")

  ## parameters ##
  expect_error(item(model = "2PL", parameters = list()),
               "Invalid item parameters")
  expect_error(item(model = "2PL", parameters = list(b = 1)),
               "Invalid item parameters")

  # Ambiguous parameters should raise error:
  expect_error(item(a = c(1, 2), b = c(-1, .2)), "Invalid parameters.")
  expect_error(item(a = c(1, 2, 3), b = c(-1, .2, 2)), "Invalid parameters.")
  expect_error(item(a = c(1, 2), b = c(-1, .2, 2)), "Invalid parameters.")
})


###############################################################################@
############################# is.Item  #########################################
###############################################################################@
# Test is.Item
test_that("Test is.Item", {
  expect_true(object = is.Item(new(
    "Item", model =  "2PL", parameters = list(a = 2, b = -.22, D = 1.702))))
  expect_true(object = is.Item(new(
    "Item", model = "M2PL", parameters = list(a = c(1,2), d = -.22, D = 1))))
})

###############################################################################@
############################# print.Item  ######################################
###############################################################################@
# Test print.Item and show.Item
test_that("Test show.Item", {
  # IRT-2PM
  item <- new("Item", model =  "2PL",
              parameters = list(b = 2.7, a = 1.3, D = 1))
  expect_output(print(item), regexp = "2PL")
  expect_output(print(item), regexp = "An object of class 'Item'.")
  # IRT-3PM
  item <- new("Item", model =  "3PL",
              parameters = list(b = 2, c = .12, a = 1.2, D = 1))
  expect_output(print(item), regexp = "1.20 2.00 0.12")
  # IRT-3PM - Print Content
  item <- new("Item", model =  "3PL", content = "C-1.2", id = "MyItem122",
              parameters = list(b = 2, c = .123, a = 1.2, D = 1.7))
  expect_output(print(item), regexp = "ID:")
  expect_output(print(item), regexp = "Content:")
  expect_output(print(item), regexp = " C-1.2", fixed = TRUE)
  # MIRT-3PM
  item <- new("Item", model =  "M3PL", content = "C-12",
              parameters = list(d = 2, c = .12, a = c(1,1.2), D = 1.7))
  expect_output(print(item), regexp = "a1   a2    c    d")
  # GRM
  item <- item(a = 1.2, b = c(-2, 1, 2), model = "GRM")
  expect_output(print(item), regexp = "GRM")
  expect_output(print(item), regexp = "a   b1   b2   b3")

  # misc field printed properly:
  item <- item(a = 1.2, b = -0.94, id = "item1", content = "Earth Science",
               misc = list(key = "C",  operational = TRUE, type = "MC",
                           enemies = c("i2", "i3"), ep = c(1.2, -2.11, 3.1),
                           np = list(x = 1, y = 2)))
  expect_output(print(item), regexp = "key")
  expect_output(print(item), regexp = "operational")
  expect_output(print(item), regexp = "type")
  expect_output(print(item), regexp = "enemies")
})


###############################################################################@
############################# item (legacy from as.Item)  ######################
###############################################################################@
# Test as.Item
test_that("Test item (legacy from as.Item)", {
  # When there is no argument, an error should be raised
  expect_error(object = item())

  ### Test numeric part ####################################################@###
  # 'item()' returns an "Item" class when an "Item" class entered.
  expect_is(object = item(new(
    "Item", model =  "2PL", parameters = list(a = 2, b = -.22, D = 1.702))),
    class = "Item")
  expect_is(object = item(a = 2, b = -.22, c = .24, d = .99), class = "Item")
  expect_is(object = item(a = 1, b = 1.2), class = "Item")
  item1 <- item(a = 2, b = -.22, c = .24, d = .99)
  expect_true(item1@parameters[2] == -.22)
  expect_is(object = item(b = c(2.11), id = "Item 12"), class = "Item")
  # -------------------------------------------------------------------------- #
  # 1PL
  item <- item(b = 1.2, D = 1)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "1PL")
  expect_equal(item@parameters$D, 1)
  item <- item(b = 1.2)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "Rasch")
  item <- item(b = 1.2, D = 1)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "1PL")

  item <- item(b = c(-1, 2, 1))
  expect_is(object = item, class = "Item")
  expect_true(item$model == 'PCM')

  # -------------------------------------------------------------------------- #
  # 2PL
  item <- item(a = 1, b = 1.2, D = 1)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "2PL")
  expect_equal(item@parameters$D, 1)
  item <- item(a = 1, b = 1.2)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "2PL")

  # TODO: The following gives error, solve this.
  item <- item(a = 1.38, b = -.81, D = 1.2, model = "2PL")
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "2PL")
  expect_equal(item$a, 1.38)
  expect_equal(item$b, -0.81)
  expect_equal(item@parameters$D, 1.2)

  # -------------------------------------------------------------------------- #
  # 3PL
  item <- item(a = 1, b = 2, c = .2, D = 1.7)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "3PL")
  #
  item <- item(a = 1, b = 2, c = .2)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "3PL")
  expect_equal(item@parameters$D, 1)
  #
  item <- item(a = 1, b = 0.02173983, c = 0, D = 0.7636914, model = "3PL")
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "3PL")
  expect_equal(item@parameters$D, 0.7636914)

  # -------------------------------------------------------------------------- #
  # 4PL
  item <- item(a = 1, b = 2, c = .2, d=.99, D = 1.1)
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "4PL")
  expect_equal(item@parameters$b, 2)
  expect_equal(item@parameters$d, .99)
  expect_equal(item@parameters$D, 1.1)

  # Default value of D
  item <- item(a = 1, b = 1.2)
  expect_equal(item@parameters$D, 1)

  # Test the IDs are "Item-xx"
  ip <- item(b = rnorm(10))
  expect_null(ip$id)

  # -------------------------------------------------------------------------- #
  # Create Items for Graded Response Model
  # The default method is GRM when no model presented and there are multiple
  # b parameter values.
  item <- item(parameters = list(a = 1.2, b = c(-2.3, -1.2, 0.4, 1.3), D = 1))
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "GRM")
  item <- item(parameters = list(a = 1.2, b = c(-2.3, -1.2, 0.4, 1.3), D = 1),
               model = "GRM")
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "GRM")

  # ignore D
  item <- item(parameters = list(a = 1.2, b = c(-2.3, -1.2, 0.4, 1.3)))
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "GRM")
  expect_equal(item@parameters$D, 1)

  # -------------------------------------------------------------------------- #
  # Create Items for Generalized Partial Credit Model
  item <- item(parameters = list(a = 1.2, b = c(-2.3, -1.2, 0.4, 1.3),
                                 D = 1.702),
                  model = "GPCM")
  expect_is(object = item, class = "Item")
  expect_equal(item@model, "GPCM")

  # -------------------------------------------------------------------------- #
  item1 <- item(b = 2.11, c = .22, a = 1, id = "Item 11",
                   content = "Math-12",
                   se_parameters = list(a = .22, b = .32, c = .1))
  expect_is(object = item1, class = "Item")
  expect_true(all(item1@id == "Item 11", item1@model == "3PL",
                  item1@parameters[2] == 2.11, item1@se_parameters[3] == .1,
                  item1@content == "Math-12"))

  # -------------------------------------------------------------------------- #
  expect_error(item(b = 1, model = "abc"))

  # -------------------------------------------------------------------------- #
  # Another thing to pass
  item1 <- item(a = 1, b = c(-1, 0, 1), D = 1)
  expect_equivalent(item1$model, 'GRM')

  # -------------------------------------------------------------------------- #
  # Correctly specify an incorrectly ordered parameter vector
  item1 <- item(b = 2.11, c = .22, a = 1)
  expect_true(item1@parameters[2] == 2.11)

  # -------------------------------------------------------------------------- #
  # Issue error if model is incorrectly specified
  expect_error(object = item(x = c(1,2,.2), model = "2PL"),
               "Incomplete parameters. Make sure to provide")
  expect_error(object = item(x = c(1,2,.2), model = "WeirdoModel"),
               "Invalid 'model' specification.")

  # -------------------------------------------------------------------------- #
  # Issue error if parameters are NULL or NA
  expect_error(object = item(NULL))
  expect_error(object = item(NA))

  # -------------------------------------------------------------------------- #
  # Check whether the names of parameters removed
  dtf <- data.frame(a = c(1, 2), xxb1 = c(-2, 1), xxb2 = c(.2, .9))
  i = 1
  item <- item(a = dtf$a[i], b = unlist(dtf[i, paste0('xxb', 1:2)]))
  expect_null(names(item$parameters$a))
  expect_null(names(item$parameters$b))

  # ################# Test matrix part #######################################@###
  # ## Test a matrix with one column
  # ip <- matrix(c(a = .12, b = 2.21), nrow = 1, dimnames = list(NULL, c("a", "b")),byrow = FALSE)
  # expect_is(object = item(ip), class = "Item")
  # ip <- item(ip)
  # expect_true(ip@parameters[2] == 2.21)
  # expect_true(is.null(ip@id))
  # ip <- item(matrix(c(a = .12, b = 2.21),
  #              nrow = 1, byrow = FALSE), id = "Geo-Item-1")
  # expect_true(ip@id == "Geo-Item-1")
  #
  # # Issue error when there is more columns than 4
  # ip <- matrix(c(a = .12, b = 2.21, c = .43, d = 1, e = 2),
  #             nrow = 1, byrow = FALSE)
  # expect_error(object = item(ip),
  #              regexp = "Model cannot have more than 4 parameters.")
  #
  # # A matrix without a name
  # ip <- item(matrix(runif(3*10), ncol = 3))
  # expect_is(object = ip, class = "list")
  # expect_true(all(sapply(ip, class) == 'Item'))
  #
  # ## Multiple rows
  # expect_error(object = item(matrix(c(1:12), nrow = 2, byrow = FALSE)),
  #              regexp = "Matrix cannot have more than 5 columns.")
  #
  # ip <- matrix(c(1,-2,.2,1, 1.3, 1.2, .35, .93), nrow = 2, byrow = TRUE,
  #              dimnames = list(NULL, c("a", "b", "c", "d")))
  # expect_is(object = item(ip), class = "list")
  #
  # expect_is(object = item(
  #   matrix(c(a = runif(10, .5, 1.5), b = runif(10, -3,3), D = rep(1.7, 10)),
  #          ncol = 3, byrow = FALSE, dimnames = list(NULL, c('a', 'b', 'D')))),
  #   class = "list")
  # ip <- item(ip)
  # expect_true(ip[[2]]@parameters$a == 1.3)
  # expect_true(is.null(ip[[1]]@id))
  #
  # expect_is(object = item(
  #   as.matrix(data.frame(b = runif(10, -3,3), D = rep(1.7, 10))),
  #   id = paste0("My-Item-",1:10), content = rep("Algebra", 10)), class = "list")
  #
  # ip <- as.matrix(data.frame(a = runif(10, .5, 1.5), b = runif(10, -3,3),
  #                            c = runif(10, 0,.3)))
  # rownames(ip) <- paste0("My-Items-",1:10)
  # expect_is(object = item(ip), class = "list")
  # rownames(ip) <- rep(NA, 10)
  # expect_is(object = item(ip), class = "list")

  # ################# Test data.frame part ###################################@###
  # ## Single line data.frames.
  # expect_is(item(data.frame(a = .12, b = 2.21)), class = "Item")
  #
  # ip <- item(data.frame(a = .12, b = 2.21), id = "Ka-2", content = "Alg1")
  # expect_is(ip, class = "Item")
  # expect_true(ip@id == "Ka-2")
  # expect_true(ip@content == "Alg1")
  #
  # ## Multi-line data.frames.
  # ip <- data.frame(a = runif(10, .5, 1.5), b = runif(10, -3,3),
  #                  c = runif(10, 0, .2))
  # expect_is(item(ip), class = "list")
  # ip1 <- item(ip, id = paste0(1:10), content = rep("Alg1",10))
  # expect_true(ip1[[9]]@id == "9")
  # expect_true(ip1[[3]]@content == "Alg1")
  #
  #
  # # Test creation of GRM objects
  # n <- sample(8:13, 1)
  # ip_dtf <- data.frame(a = rlnorm(n, 0, .3), b1 = rnorm(n, -1))
  # ip_dtf$b2 <- ip_dtf$b1 + runif(n, 0.1)
  # ip_dtf$b3 <- ip_dtf$b2 + runif(n, 0.1)
  # # args <- list(ip_dtf, model = "GRM")
  # # Test Data Frame
  # item_list <- item(ip_dtf, model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # i <- sample(1:n, 1)
  # expect_equal(item_list[[i]]$a, ip_dtf[i, "a"])
  # # Test Matrix
  # item_list <- item(as.matrix(ip_dtf), model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # i <- sample(1:n, 1)
  # expect_equal(item_list[[i]]$a, ip_dtf[i, "a"])
  # # # Test Partial Credit Model
  # # # This test stop working on 2020-09-06
  # # item_list <- item(as.matrix(ip_dtf), model = "PCM")
  # # expect_true(all(sapply(item_list, is.Item)))
  # # expect_true(all(sapply(item_list, function(x) x$model == "PCM")))
  # # expect_null(item_list[[i]]$a)
  # # Test "GPCM"
  # item_list <- item(as.matrix(ip_dtf), model = "GPCM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GPCM")))
  # i <- sample(1:n, 1)
  # expect_equal(item_list[[i]]$a, ip_dtf[i, "a"])
  #
  # # Test D parameter
  # ip_dtf$D <- 1.1
  # # Test Data Frame
  # item_list <- item(ip_dtf, model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # expect_equal(item_list[[sample(1:n, 1)]]$D, 1.1)
  # # Test Matrix
  # item_list <- item(as.matrix(ip_dtf), model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # expect_equal(item_list[[sample(1:n, 1)]]$D, 1.1)
  #
  # # Test Content, only valid for data.frame
  # item_list <- item(cbind(ip_dtf, content = rep(c("G", "M"), len = n)),
  #                      model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # expect_equal(item_list[[sample(1:n, 1)]]$D, 1.1)
  # expect_equal(item_list[[3]]$content, "G")
  # expect_equal(item_list[[4]]$content, "M")
  #
  # # Test id, only valid for data.frame
  # item_list <- item(cbind(ip_dtf, content = rep(c("G", "M"), len = n),
  #                            id = paste0("ii", 1:n)),
  #                      model = "GRM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GRM")))
  # expect_equal(item_list[[sample(1:n, 1)]]$D, 1.1)
  # expect_equal(item_list[[3]]$content, "G")
  # expect_equal(item_list[[4]]$content, "M")
  # i <- sample(1:n, 1)
  # expect_equal(item_list[[i]]$id, paste0("ii", i))
  #
  # # Test GPCM
  # item_list <- item(cbind(ip_dtf, content = rep(c("G", "M"), len = n),
  #                            id = paste0("ii", 1:n)),
  #                      model = "GPCM")
  # expect_true(all(sapply(item_list, is.Item)))
  # expect_true(all(sapply(item_list, function(x) x$model == "GPCM")))
  # expect_equal(item_list[[sample(1:n, 1)]]$D, 1.1)
  # expect_equal(item_list[[3]]$content, "G")
  # expect_equal(item_list[[4]]$content, "M")
  # i <- sample(1:n, 1)
  # expect_equal(item_list[[i]]$id, paste0("ii", i))
  #
  # # # Test PCM
  # # # This test stop working on 2020-09-06
  # # item_list <- item(cbind(ip_dtf, content = rep(c("G", "M"), len = n),
  # #                            id = paste0("ii", 1:n)),
  # #                      model = "PCM")
  # # expect_true(all(sapply(item_list, is.Item)))
  # # expect_true(all(sapply(item_list, function(x) x$model == "PCM")))
  # # expect_null(item_list[[sample(1:n, 1)]]$D)
  # # expect_null(item_list[[sample(1:n, 1)]]$a)
  # # expect_equal(item_list[[3]]$content, "G")
  # # expect_equal(item_list[[4]]$content, "M")
  # # i <- sample(1:n, 1)
  # # expect_equal(item_list[[i]]$id, paste0("ii", i))

})


###############################################################################@
############################# '$' method #######################################
###############################################################################@
# Test '$' method
test_that("Test $ method ", {
  item <- item(a = 1, b = .2, c = .3, id = 'myid', content = 'Algebra')
  expect_equivalent(item$a, 1)
  expect_equivalent(item$b, .2)
  expect_equivalent(item$c, .3)
  expect_equivalent(item$D, 1)
  expect_null(item$d) # Non-exsitent parameter
  expect_equivalent(item$id, 'myid')
  expect_equivalent(item$parameters, item@parameters)
  expect_equivalent(item$content, 'Algebra')
  expect_equivalent(item$model, '3PL')
  expect_equal(item$max_score, 1)

  # -------------------------------------------------------------------------- #
  item <- generate_item(model = "GRM", n_categories = 5)
  expect_equal(item$max_score, 4)

})

###############################################################################@
############################# '$<-' method #####################################
###############################################################################@
# Test '$<-' method
test_that("Test $<- method ", {
  item <- item(a = 1, b = .2, c = .3, id = 'myid', content = 'Algebra')
  expect_equivalent(item$a, 1)
  expect_equivalent(item$b, .2)
  expect_equivalent(item$c, .3)
  expect_equivalent(item$D, 1)
  expect_equivalent(item$id, 'myid')
  expect_equivalent(item$parameters, item@parameters)
  expect_equivalent(item$content, 'Algebra')
  expect_equivalent(item$model, '3PL')
  item$a <- 2
  item$b <- -1
  item$c <- .41
  item$D <- 1.1
  item$d <- .99   # Non-exsitent parameter
  item$id <- "New_ID"
  item$content <- "Geometry"
  expect_equivalent(item$a, 2)
  expect_equivalent(item$b, -1)
  expect_equivalent(item$c, .41)
  expect_equivalent(item$D, 1.1)
  expect_null(item$d)  # Non-exsitent parameter
  expect_equivalent(item$id, 'New_ID')
  expect_equivalent(item$content, 'Geometry')
  # expect_equivalent(item$model, '3PL')

  # -------------------------------------------------------------------------- #
  # One cannot change the item's model or parameters
  item <- item(a = 1, b = .2, c = .3, id = 'myid', content = 'Algebra')
  expect_equal(item$model, '3PL')
  item$model <- "Rasch"
  expect_equal(item$model, "3PL")

  # -------------------------------------------------------------------------- #
  # Item parameters cannot be updated if they are not in agreement with the
  # underlying model.
  expect_warning(item$parameters <- list(b = 1), "Invalid item parameters.")
  # Parameters should be completed
  expect_warning(item$parameters <- list(a = .5, b = .4, c = .2),
                 "Invalid item parameters.")
  # Parameters should be valid
  expect_error(item$parameters <- list(a = "x", b = "y", c = "z", D = "d"),
               "Invalid parameters")
  # -------------------------------------------------------------------------- #
  # Item parameters can be updated with acceptable list of parameters.
  item$parameters <- list(a = .5, b = .4, c = .2, D = 1.99)
  expect_equivalent(item$a, .5)
  expect_equivalent(item$b, .4)
  expect_equivalent(item$c, .2)
  expect_equivalent(item$D, 1.99)
  })
