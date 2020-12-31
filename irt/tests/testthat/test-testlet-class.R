
# library(testthat)
# Test setValidity of "Item" class
test_that("Test setValidity of 'Testlet' class", {
  # Each element of the "Testlet" should be 'Item' class
  item1 <- item(a = 1.12, b = -2.1, c = 0.28)
  item2 <- item(a = 2, b = 3.2, c = 0.21)
  item3 <- c(a = 1.2, b = 2.8, c = 0.12) # This is not 'Item' class
  item4 <- item(a = 1.12, b = -1.23, c = .2, id = "I-21")
  item5 <- item(a = 0.84, b = 2.23, c = .25, id = "I-22")
  # Create a new Itempool
  expect_is(object = new("Testlet", item_list = itempool(item4, item5)),
            class = "Testlet")
  expect_is(object = new("Testlet", item_list = itempool(item4)),
            class = "Testlet")



  ########################### id ##############################################@
  # Testlet id should have length 1 or NULL
  expect_error(object = new("Testlet", model = "BTM",
                            item_list = itempool(item4, item5),
                            parameters = NULL,
                            id = c("Item1", "Item2")),
               regexp = "Invalid testlet 'id'.")
  expect_error(
    object = new("Testlet", model = "BTM", item_list = itempool(item4, item5),
                 parameters = NULL, id = NA),
    regexp = "assignment of an object of class \"logical\" is not valid for")
  expect_error(object = new("Testlet", model = "BTM",
                            item_list = itempool(item4, item5),
                            parameters = NULL,
                            id = as.character(NA)),
               regexp = "Invalid testlet 'id'.")

  ########################### model ###########################################@
  # Default model name is "BTM":
  t1 <- new("Testlet", item_list = itempool(item4))
  expect_equal(t1@model, "BTM")
  # Incorrect model name:
  expect_error(object = new("Testlet",
                            item_list = itempool(item4, item5),
                            model = "XYZ_MODEL"),
               regexp = "Invalid model value.")

  ########################### item_list #######################################@
  # # All of the elements of Itempool class should be 'Item' class.
  # expect_error(new("Testlet", item_list = itempool(item1, item2, item3, 12)),
  #              regexp = "Invalid elements in item_list.")

  # Incorrect item_list Class
  expect_error(new("Testlet", item_list = item4))
  expect_error(new("Testlet", item_list = list()),
               regexp = "'item_list' should be an 'Itempool' object.")

  ########################### parameters ######################################@
  # Parameters cannot be NULL or NA if model is other than 'BTM'.
  expect_error(object = new("Testlet", parameters = NULL, model = 'RTM',
                            item_list = itempool(item4, item5)),
               regexp = "Invalid testlet parameter.")
  expect_error(object = new("Testlet", parameters = list(mean = 12),
                            model = 'RTM', item_list = itempool(item4, item5)),
               regexp = "Invalid testlet parameter.")
  expect_error(object = new("Testlet", parameters = list(x = 12, y = 1),
                            model = 'RTM', item_list = itempool(item4, item5)),
               regexp = "Invalid parameter names. Parameter names for RTM")

  # expect_error(object = new("Item", model = '1PL',
  #                           parameters = list(a = 2, b = -.22, c = NA)),
  #              regexp = "Invalid parameter. Item parameters cannot be NULL or NA.")





})
