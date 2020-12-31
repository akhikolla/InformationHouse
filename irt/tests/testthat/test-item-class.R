# library(testthat)
# Test setValidity of "Item" class

test_that("Test setValidity of Item class", {
  # irt items
  expect_is(object = new("Item", model = "2PL",
                         parameters = list(a = 2, b = -.22, D = 1.702)),
            class = "Item")
  expect_is(object = new("Item", model = "3PL",
                         parameters = list(a = 2, b = -.22, c = 0, D = 1.702)),
            class = "Item")
  expect_is(object = new("Item", model = "3PL",
                         parameters = list(a = 2, b = -.22, c = 0, D = 1.702),
                         id = "MyItem12", content = "2"), class = "Item")
  expect_is(object = new("Item", model = "3PL",
                         parameters = list(a = 2, b = -.22, c = 0, D = 1.702),
                         id = "MyItem12", content = "2"), class = "Item")
  expect_is(object = new("Item", model = "1PL",
                         parameters = list(b = 2, D = 1.702),
                         id = "MyItem12", content = "1.2"), class = "Item")
  expect_is(object = new("Item", model = "1PL",
                         parameters = list(b = 2, D = 1.702),
                         se_parameters = list(b = .12),
                         id = "MyItem12", content = "1.2"), class = "Item")
  expect_is(object = new("Item", model = "3PL",
                         parameters = list(b = 2, c = .12, a = 1.2, D = 1.702),
                         se_parameters = list(c = .12, a = .15, b = .17),
                         id = "MyItem12", content = "1.2"), class = "Item")

  # Test the default values of "Item"
  item <- new("Item")
  expect_is(item, "Item")
  expect_equal(item@model, "2PL")
  expect_equal(item@parameters$a, 1)
  expect_equal(item@parameters$b, 0)
  expect_null(item@se_parameters)
  expect_null(item@content)
  expect_null(item@id)
  expect_null(item@misc)

  # -------------------------------------------------------------------------- #
  # mirt items
  expect_is(object = new(Class = "Item", model = "M1PL",
                         parameters = list(d = 1.2, D = 1)), class = "Item")
  expect_is(object = new(Class = "Item", model = "M2PL",
                         parameters = list(d = 1.2, a = c(1, 1.2), D = 1)),
            class = "Item")
  expect_is(object = new(Class = "Item", model = "M3PL",
                         parameters = list(d = 1.2, a = c(1, 1.2), c = .1,
                                           D = 1)),
            class = "Item")

  # -------------------------------------------------------------------------- #
  # Generalized Partial Credit items
  expect_is(object = new(Class = "Item", model = "GPCM",
                         parameters = list(b = c(-3, -.5, 1.5), a = 1, D = 1)),
            class = "Item")

  # -------------------------------------------------------------------------- #
  # Reparametrized Generalized Partial Credit items
  expect_is(object = new(Class = "Item", model = "GPCM2",
                         parameters = list(d = c(-3, -.5, 1.5), b = 0.4, a = 1,
                                           D = 1)),
            class = "Item")
  # An error will be raised if user switches b with d (b should have length 1)
  expect_error(
    object = new(Class = "Item", model = "GPCM2",
                 parameters = list(b = c(-3, -.5, 1.5), d = 0.4, a = 1, D = 1)),
    "Invalid parameters.")


  # -------------------------------------------------------------------------- #
  # The item parameters should saved in ascending order by name
  a1 <- rlnorm(1, 0, .3)
  b1 <- rnorm(1)
  c1 <- runif(1, 0, .3)
  D1 <- 1.7
  i1 <- new("Item", parameters = sample(list(b = b1, a = a1, D = D1, c = c1)),
            model = "3PL")
  expect_equal(i1@parameters[[1]], a1)
  expect_equal(i1@parameters[[2]], b1)
  expect_equal(i1@parameters[[3]], c1)
  expect_equal(i1@parameters[[4]], D1)

  ################# ERROR CASES ###############################################@
  # -------------------------------------------------------------------------- #
  ###### id ######@
  # id should have length 1 or NULL
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(a = 2, b = -.22, D = 1.702),
                            id = c("Item1", "Item2")),
               regexp = "Invalid item id.")


  # -------------------------------------------------------------------------- #
  ###### model ######@
  # Model should be specified, otherwise it will be 2PL which will not be
  # valid
  expect_error(object = new("Item", parameters = list(a = 1,b = 0, c = 0,
                                                      d = 1)),
               regexp = "Invalid parameter names.")
  expect_error(object = new(Class = "Item",
                            parameters = list(d = 1.2, a = c(1, 1.2))),
               regexp = "Invalid parameter names.")
  expect_error(object = new("Item", model = "IRTXYZ"),
               regexp = "Invalid model value.")

  # -------------------------------------------------------------------------- #
  # Incorrect model specification
  expect_error(object = new("Item", model =  "4PL",
                            parameters = list(a = 2, b = -.22, D = 1.702)),
               regexp = paste0("Invalid parameter names. Parameter names of ",
                               "Item class should be unique and complete. "))
  expect_error(object = new(Class = "Item", model = "M2PL",
                         parameters = list(d = 1.2, a = c(1, 1.2), c = .1)))

  # -------------------------------------------------------------------------- #
  # Model length should be 1
  expect_error(object = new("Item", model =  c("2PL", "3PL"),
                            parameters = list(a = 2, b = -.22)),
               regexp = "Invalid model")


  # -------------------------------------------------------------------------- #
  ###### parameters ######@
  # Parameters cannot be NULL or NA
  expect_error(object = new("Item", parameters = NULL),
               regexp = "Invalid parameters. Parameters should be a 'list'.")
  expect_error(
    object = new("Item", model = '1PL',
                 parameters = list(a = 2, b = -.22, c = NA)),
    regexp = "Invalid parameter. Item parameters cannot be NULL or NA.")
  # Item parameters cannot have NA in their values
  expect_error(
    object = new("Item", model =  "3PL",
                 parameters = list(a = 2, b = -.22, c = NA)),
    regexp = "Invalid parameter. Item parameters cannot be NULL or NA.")
  # Item parameter values cannot be non-numeric
  expect_error(object = new("Item", model =  "1PL",
                            parameters = list(b = "ABC", D = 1.1)),
               regexp = "Invalid parameters. All parameters should be numeric.")


  # -------------------------------------------------------------------------- #
  # c Parameter cannot be larger than 1 or smaller than 0.
  expect_error(object = new("Item",
                            parameters = list(a = 1, b = 2, c = 3, D = 1.702),
                            model = "3PL"),
               regexp = "Invalid parameters.")
  expect_error(
    object = new("Item", parameters = list(a = 1, b = 2, c = -0.1, D = 1.702),
                 model = "3PL"), regexp = "Invalid parameters.")
  expect_error(
    object = new("Item", parameters = list(a = 1, b = 2, c = -0.1, D = 1.702),
                 model = "3PL"), regexp = "Invalid parameters.")
  expect_error(object = new(Class = "Item", model = "M3PL",
                            parameters = list(d = 1.2, a = c(1, 1.2), c = 1.1,
                                              D = 1)),
               regexp = "Invalid parameters.")

  # -------------------------------------------------------------------------- #
  # d Parameter cannot be larger than 1 or smaller than 0. Also cannot be
  # smaller than c parameter.
  expect_error(object = new(
    "Item", parameters = list(a = 1, b = 2, c = .2, d = 3, D = 1.702),
    model = "4PL"),
    regexp = "Invalid parameters.")
  expect_error(object = new(
    "Item",
    parameters = list(a = 1, b = 2, c = .2, d = -0.3, D = 1.702),
    model = "4PL"),
    regexp = "Invalid parameters.")

  # -------------------------------------------------------------------------- #
  # Problematic naming of parameters
  expect_error(object = new("Item", model = "3PL",
                            parameters = list(2.11, c = .22, a = 1, D = 1.702)),
               regexp = "Invalid parameter names.")
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(c = .22, a = 1, D = 1.702)),
               regexp = "Invalid parameter names.")
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(b = .22, b = 1, D = 1.702)),
               regexp = "Invalid parameter names.")
  expect_error(object = new(Class = "Item", model = "M2PL",
                            parameters = list(e = 1.2, a = c(1, 1.2), D = 1)),
               regexp = "Invalid parameter names.")
  expect_error(
    object = new("Item", model =  "4PL",
                 parameters = list(c = .22, a = 1, d = .9, D = 1.702)),
    regexp = paste0("Invalid parameter names. Parameter names of ",
                    "Item class should be unique and complete. "))
  expect_error(object = new("Item", model =  "M3PL",
                            parameters = list(a = c(1, 1.2) , b = 2.11, c = .22,
                                              D = 1)),
               regexp = "Invalid parameter names.")

  # -------------------------------------------------------------------------- #
  # Length of parameters should be 1.
  expect_error(
    object = new("Item", model =  "3PL",
                 parameters = list(a = c(1, 1.2) , b = 2.11, c = .22,
                                   D = 1.702)),
    regexp = paste0("Invalid parameters."))
  expect_error(
    object = new("Item", model =  "3PL",
                 parameters = list(a = 1 , b = 2.11, c = c(.22, 1.2),
                                   D = 1.702)),
    regexp = paste0("Invalid parameters."))
  expect_error(
    object = new("Item", model =  "M3PL",
                 parameters = list(a = c(1, 1.2) , d = 2.11, c = c(.22, 1.2),
                                   D = 1.702)),
    regexp = "Invalid parameters.")


  # -------------------------------------------------------------------------- #
  # Length of each parameter element should correspond to the size argument of
  # Pmodels.
  expect_error(
    object = new("Item", model =  "GRM",
                 parameters = list(a = c(1, 2.2) , b = c(-1, 0, 2.11),
                                   D = 1.702)),
    regexp = paste0("Invalid parameters."))
  expect_error(
    object = new("Item", model =  "GPCM",
                 parameters = list(a = c(0.2, 1, 2.2) , b = c(-1, 0, 2.11, 3),
                                   D = 1.702)),
    regexp = paste0("Invalid parameters."))


  # -------------------------------------------------------------------------- #
  ###### se_parameters ######@

  # se_parameters should be numeric
  expect_is(object = new(
    "Item", model =  "4PL",
    parameters = list(a = 2, b = -.22, c = 0.2, d = .99, D = 1.702),
    se_parameters = list(a = Inf, b = 0, c = .2, d = .99)),
    class = "Item")

  # NA is an acceptable se_parameter
  expect_is(new("Item", model =  "2PL",
                parameters = list(a = 2, b = -.22, D = 1.702),
                se_parameters = list(a = Inf, b = as.numeric(NA))),
            "Item")

  expect_is(new("Item", model =  "GRM",
                parameters = list(a = 2, b = c(-1, 0, 1), D = 1.702),
                se_parameters = list(a = Inf, b = rep(as.numeric(NA), 3))),
            "Item")

  # Even logical NA is acceptable
  expect_is(new("Item", model =  "GRM",
                parameters = list(a = 2, b = c(-1, 0, 1), D = 1.702),
                se_parameters = list(a = Inf, b = rep(NA, 3))),
            "Item")


  # se_paramters for polytomous items
  expect_is(new("Item", model =  "GPCM",
                parameters = list(a = 2, b = c(-.22, 1.2, 2.3), D = 1.702),
                se_parameters = list(a = 1, b = c(.13, .2, .3))),
            "Item")
  # The order of se_parameters can be different
  expect_is(new("Item", model =  "GPCM",
                parameters = list(a = 2, b = c(-.22, 1.2, 2.3), D = 1.702),
                se_parameters = list(b = c(.13, .2, .3), a = .3)),
            "Item")

  # se_parameters should be larger than 0
  expect_error(
    object = new("Item", model =  "2PL",
                 parameters = list(a = 2, b = -.22, D = 1.702),
                 se_parameters = list(a = .14, b = -.1)),
    regexp = paste0("Invalid 'se_parameters' values. Standard error values of ",
                    "item parameters ('se_parameters') cannot be smaller than ",
                    "0."),
    fixed = TRUE)

  # se_parameters should match item parameters
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(a = 2, b = -.22, D = 1.702),
                            se_parameters = list(a = .14, b = .13, c = .24)),
               regexp = paste0("Invalid 'se_parameters' values.\n",
                               "'se_parameters' should be"),
               fixed = TRUE)

  # se_parameters should match item parameters
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(a = 2, b = -.22, D = 1.702),
                            se_parameters = list(b = .13)),
               regexp = paste0("Invalid 'se_parameters' values."),
                               fixed = TRUE)

  # se_parameters cannot be an empty list
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(a = 2, b = -.22, D = 1.702),
                            se_parameters = list()),
               regexp = paste0("Invalid 'se_parameters' values."),
                               fixed = TRUE)

  # se_parameters elements cannot be empty, and should have
  expect_error(object = new("Item", model =  "2PL",
                            parameters = list(a = 2, b = -.22, D = 1.702),
                            se_parameters = list(a = Inf, b = numeric(0))),
               regexp = paste0("Invalid 'se_parameters' values."),
                               fixed = TRUE)

  # Item parameter values cannot be non-numeric
  expect_error(
    object = new("Item", model =  "1PL", parameters = list(b = 1, D = 1),
                 se_parameters = list(b = "ABC")),
    regexp = paste0("Invalid 'se_parameters' values. All standard error values",
                    " of item parameters ('se_parameters') should be numeric."),
    fixed = TRUE)


  # None of the se_parameter values can be NULL individually
  expect_error(new("Item", model =  "2PL",
                parameters = list(a = 2, b = -.22, D = 1.702),
                se_parameters = list(a = .2, b = NULL)),
            paste0("Invalid 'se_parameters' values. Individual elements of ",
                   "'se_parameters' cannot be NULL."))

  # For polytomous items, the standard errors should have the same number of
  # elements as the number of parameters
  expect_error(object = new("Item", model =  "GPCM",
                            parameters = list(a = 2, b = c(-.22, 1.2, 2.3),
                                              D = 1.702),
                            se_parameters = list(a = 1, b = .13)),
               regexp = paste0("Invalid 'se_parameters' values."),
               fixed = TRUE)
})

