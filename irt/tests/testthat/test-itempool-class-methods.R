# library(testthat)

###############################################################################@
############################# itempool #########################################
###############################################################################@

test_that("Test itempool", {
  ################# Test item_pool part ####################################@###
  ### Itempool object
  ip <- itempool(list(item(b = 1, id = "K1"), item(a = 1, b = 2, id = "K2")))
  expect_is(itempool(ip), "Itempool")
  ip1 <- itempool(ip, id = c('i1', 'i2'))
  expect_equal(ip1$id, c('i1', 'i2'))

  ip1 <- itempool(ip, id = c('a1', 'a2'), content = c("A", "G"))
  expect_equal(ip1$id, c('a1', 'a2'))
  expect_equivalent(ip1$content, c('A', 'G'))
  expect_equal(names(ip1$content)[2], "a2")
  # When there is no named argument return the Itempool, ignore the arguments
  ip1 <- itempool(ip, 3, "b", c(1, 2,3))
  expect_equal(ip, ip1)
  # Change the D parameter
  ip <- itempool(item_list = list(item(a = .9, b = -1.7, id = "K1"),
                                  item(a = 1, b = 2, id = "K2")))
  expect_equivalent(ip$D, c(1, 1))
  ip1 <- itempool(ip, D = 1.2)
  expect_equivalent(ip1$D, c(1.2, 1.2))
  ip1 <- itempool(ip, D = c(3.4, 5.5))
  expect_equivalent(ip1$D, c(3.4, 5.5))

  # -------------------------------------------------------------------------- #
  ## item
  ip <- item(a = 1, b = 2)
  # Item do not have an ID automatically
  expect_null(ip$id)
  ip1 <- itempool(ip)
  expect_is(ip1, "Itempool")
  # Default id is 'Item-1, itempool() function adds an ID
  expect_equal(ip1$id, 'Item-1')

  # -------------------------------------------------------------------------- #
  # Note: 2020/02/03: following feature deprecated
  ip1 <- itempool(ip, id = 'i1', content = 'algebra')
  expect_is(ip1, "Itempool")
  expect_equal(ip1$id, 'i1')
  expect_equivalent(ip1$content, 'algebra')
  # The following would be great if it works.
  expect_equivalent(ip$D, 1)
  ip1 <- itempool(ip, D = 3.2)
  expect_equivalent(ip1$D, 3.2)

  # -------------------------------------------------------------------------- #
  ip <- itempool(b = c(-1, 0.2, 1.1), model = "Rasch")
  expect_equal(ip[[1]]$model, "Rasch")

  # -------------------------------------------------------------------------- #
  ## numeric or integer
  ip = itempool(a =  1, b = .2, D = 1.1, content = 'dd')
  expect_is(ip, "Itempool")
  expect_equal(length(ip), 1)
  expect_equivalent(ip$content, 'dd')
  # Another example with multiple items
  ip <- itempool(a = c(1, 2), b = c(-1, .2), content = "Algebra", D = 3)
  expect_is(ip, "Itempool")
  expect_equal(length(ip), 2)
  expect_equivalent(ip$content, c("Algebra", "Algebra"))
  expect_equal(names(ip$content), ip$id)
  expect_equivalent(ip$D, c(3, 3))

  # -------------------------------------------------------------------------- #
  ## List
  ip_list <- generate_ip(n = 3, output = "list")
  expect_is(itempool(ip_list), "Itempool")
  expect_equal(length(itempool(ip_list)), 3)
  # Error:
  expect_error(itempool(list(b = c(1,2,3))),
               regexp = paste0("Invalid elements. All elements of the list ",
                               "should be an 'Item' or 'Testlet' object."))

  # An item without an id will be automatically assigned an id.
  item_no_id <- generate_item()
  item_no_id@id <- NULL
  expect_true(validObject(item_no_id))
  ip_list <- list(generate_item(id = "abc1"), item_no_id)
  ip <- itempool(ip_list)
  expect_is(ip, "Itempool")
  expect_equal(ip[[2]]$id, "Item-1")


  # -------------------------------------------------------------------------- #
  # A list of items in which some of items don't have 'id's, the function
  # should automatically assign names
  ip_list <- vector('list', 3)
  ip_list[[1]] <- item(b = 1)
  ip_list[[2]] <- item(a = .99, b = 1)
  ip_list[[3]] <- item(a = .99, b = 1, c = .2, id = 'my_item1')
  expect_true(is.null(ip_list[[1]]$id))
  expect_true(is.null(ip_list[[2]]$id))
  expect_false(is.null(ip_list[[3]]$id))
  expect_is(ip_list, "list")
  ip <- itempool(ip_list)
  expect_is(ip, "Itempool")
  expect_false(any(is.null(ip$id)))
  expect_equal(ip$id[3], ip_list[[3]]$id)

  # -------------------------------------------------------------------------- #
  ## data.frame or matrix
  n <- 10
  ip <- itempool(data.frame(a = runif(n, .5, 1.5), b = rnorm(n)))
  expect_is(ip, "Itempool")
  expect_equivalent(ip$D, rep(1, n))
  expect_equivalent(ip$model, rep('2PL', n))

  # -------------------------------------------------------------------------- #
  n <- 10
  ip <- itempool(data.frame(a = runif(n, .5, 1.5), b = rnorm(n), D = 1))
  expect_is(ip, "Itempool")
  expect_equivalent(ip$D, rep(1, n))
  expect_equivalent(ip$model, rep('2PL', n))

  # -------------------------------------------------------------------------- #
  ## Matrix
  n <- 10
  ip_matrix <- matrix(runif(3*n), ncol = 3, dimnames = list(NULL, c("a", "b", "c")))
  ip <- itempool(ip_matrix)
  expect_is(ip, "Itempool")
  expect_equivalent(ip$model[sample(1:n, 1)], "3PL")


  # Test the IDs are "Item-xx"
  ip <- itempool(b = rnorm(10))
  expect_true(all(ip$id == paste0("Item-", 1:10)))
  ip <- itempool(a = runif(10, 1, 2), b = rnorm(10))
  expect_true(all(ip$id == paste0("Item-", 1:10)))

  #############@###
  expect_equal(itempool(new(
    "Itempool",
    item_list = list(K1 = item(b = 1, id = "K1"),
                     K2 = item(a = 1, b = 2, id = "K2"))),
    id = c("E1", "E2"))[[1]]@id, "E1")
  expect_equal(itempool(new(
    "Itempool",
    item_list = list(K1 = item(b = 1, id = "K1"),
                     K2 = item(a = 1, b = 2, id = "K2"))),
    content = c("E1", "E2"))[[1]]@content, "E1")

  # -------------------------------------------------------------------------- #
  ## Testlet
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2")), id = "t1")
  t2 <- testlet(itempool(b = 2:4, id = c("t2-i1", "t2-i2", "t2-i3")),
                   id = "t2")
  i1 <- item(b = -1, id = "i1")
  i2 <- item(b = 0, id = "i2")
  expect_is(itempool(i1, t1), "Itempool")
  expect_is(itempool(t1, i1), "Itempool")
  expect_is(itempool(t1, i1, t2), "Itempool")
  expect_is(itempool(i1, t1, t2), "Itempool")
  expect_is(itempool(t1), "Itempool")
  expect_is(itempool(t1, t2), "Itempool")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    item_list = list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  expect_is(itempool(ip), "Itempool")
  expect_is(itempool(ip@item_list), "Itempool")

  # -------------------------------------------------------------------------- #
  # misc fields correctly specified.
  ip <- itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                 misc = list(list(sympson_hetter_k = .8, form = "b3"),
                             list(sympson_hetter_k = .9)))
  expect_equal(ip[[1]]@misc$sympson_hetter_k, 0.8)
  expect_equal(ip[[1]]@misc$form, "b3")
  expect_equal(ip[[2]]@misc$sympson_hetter_k, 0.9)

  # Add the same misc value to all of the items
  ip <- itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                 misc = list(sympson_hetter_k = .8))
  expect_equal(ip[[1]]@misc$sympson_hetter_k, 0.8)
  expect_equal(ip[[2]]@misc$sympson_hetter_k, 0.8)

  # Add the same misc value to all of the items
  ip <- itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                 misc = list(list(sympson_hetter_k = .8, form = "Form1")))
  expect_equal(ip[[1]]@misc$sympson_hetter_k, 0.8)
  expect_equal(ip[[2]]@misc$sympson_hetter_k, 0.8)

  ################# Test item part #########################################@###
  # Put name if there is no name.
  expect_is(itempool(item(b = 1)), "Itempool")
  expect_is(itempool(item(a = 1, b = 1.2, id = "It-12")), "Itempool")
  expect_is(itempool(item(a = 1, b = 1.2, id = "It-12", content = "A24")),
            "Itempool")
  expect_is(itempool(item(a = 1, b = 1.2, c = .23)), "Itempool")

  expect_is(itempool(item(a = 1, b = 1.2, c = .23, id = "jeb1")), "Itempool")

  expect_equal(itempool(item(a = 1, b = 1.2, c = .23, id = "jeb1"))$id, "jeb1")
  expect_equivalent(itempool(item(a = 1, b = 1.2, c = .23),
                             content = "A1")$content, "A1")

  ################# Test numeric part ######################################@###
  expect_error(itempool(1, model = "1PL"))
  expect_error(itempool(c(1, 1.2), id = "It-12"),
               "Itempool object cannot be created for '' model.")
  expect_error(itempool(c(1, 1.2,.2,.93), id = "It-12", content = "A24"),
               "Itempool object cannot be created for '' model.")

  # When only b parameter entered without a model, it is 1PL
  ip <- itempool(b = rnorm(5))
  expect_equal(ip[[1]]$model, "1PL")

  ip <- itempool(a = rnorm(3), b = rnorm(3))
  expect_equal(ip[[1]]$model, "2PL")

  ip <- itempool(a = rlnorm(3), b = rnorm(3), c = runif(3, 0, .3))
  expect_equal(ip[[1]]$model, "3PL")

  ip <- itempool(a = rlnorm(3), b = rnorm(3), c = runif(3, 0, .3),
                 d = runif(3, .8, 1))
  expect_equal(ip[[1]]$model, "4PL")

  ################# Test matrix part #######################################@###
  ip <- matrix(c(a = .12, b = 2.21), nrow = 1,
               dimnames = list(NULL, c("a", "b")), byrow = FALSE)
  expect_is(itempool(ip), "Itempool")
  expect_is(itempool(ip, id = "Kan12", content = "kn1"), "Itempool")
  ip <- matrix(c(1,-2,.2,1, 1.3, 1.2, .35, .93), nrow = 2, byrow = TRUE,
               dimnames = list(NULL, c("a", "b", "c", "d")))
  expect_is(itempool(ip), "Itempool")
  expect_equal(itempool(ip, content = c("c1", "c2"))[[1]]@content, "c1")

  ip <- matrix(1:10)
  colnames(ip) <- "b"
  expect_equal(itempool(ip, content = rep("Alg1", 10))[[4]]@content, "Alg1")
  ################# Test data.frame part ###################################@###

  n <- sample(5:20,1)
  ip <- data.frame(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0,.3))
  expect_is(itempool(ip), "Itempool")
  expect_is(itempool(ip), "Itempool")


  # -------------------------------------------------------------------------- #

  n <- sample(10:20, 1)
  ipdf <- data.frame(a = rlnorm(n, 0, .3), b1 = rnorm(n, -1))
  ipdf$b2 <- ipdf$b1 + runif(n, 0.1)
  ipdf$b3 <- ipdf$b2 + runif(n, 0.1)
  ipdf$D <- 1.7
  ipdf$content = rep(c("G", "M"), len = n)
  rownames(ipdf) <- paste0("item--", 1:n)

  ip <- itempool(ipdf, model = "GRM")
  expect_is(ip, "Itempool")
  expect_true(all(sapply(ip$item_list, function(x) x$model == "GRM")))
  i <- sample(1:n, 1)
  expect_equal(ip$item_list[[i]]$a, ipdf[i, "a"])
  i <- sample(1:n, 1)
  expect_equal(ip$item_list[[i]]$b[2], ipdf[i, "b2"])
  expect_null(names(ip$b))
  expect_equal(names(ip$a), ip$id)
  expect_equal(names(ip$D), ip$id)
  expect_equal(names(ip$content), ip$id)
  expect_equivalent(ip$D[i], 1.7)
  expect_equal(ip$item_list[[3]]$content, "G")
  expect_equal(ip$item_list[[4]]$content, "M")
  i <- sample(1:n, 1)
  # The following feature (data.frame row names as item ids) discontinued on
  # 2020-09-06
  # expect_equal(ip$id[i], paste0("item--", i))

  # Check whether id is updated when id is explicitly put in the data.frame
  ipdf$id = paste0("ii", 1:n)
  ip <- itempool(ipdf, model = "GRM")
  i <- sample(1:n, 1)
  expect_equal(ip$id[i], paste0("ii", i))

  # -------------------------------------------------------------------------- #
  # Both id, Id, ID, iD should be acceptable as data frame column name, but if
  # multiple of them are present, "ID" or "Id" or "iD" should be ignored.
  n <- sample(4:7, 1)
  ids <- sample(letters, n)
  ip <- itempool(data.frame(id = ids, b = rnorm(n)))
  expect_equal(ip$id, ids)
  ip <- itempool(data.frame(ID = ids, b = rnorm(n)))
  expect_equal(ip$id, ids)
  ip <- itempool(data.frame(Id = ids, b = rnorm(n)))
  expect_equal(ip$id, ids)
  ip <- itempool(data.frame(iD = ids, b = rnorm(n)))
  expect_equal(ip$id, ids)


  ip <- itempool(data.frame(b = rnorm(n)), ID = ids)
  expect_equal(ip$id, ids)
  ip <- itempool(data.frame(b = rnorm(n)), iD = ids)
  expect_equal(ip$id, ids)

  # -------------------------------------------------------------------------- #
  # The following data frame or tibble cannot be created when id's are
  # duplicated. A more informative error should be issued.
  n <- sample(4:7, 1)
  ipdf <- data.frame(id = sample(letters[1:3], n, TRUE),
                       a = rlnorm(n, 0, .3), b = rnorm(n),
                       content = sample(c("Geo", "Alg"), n, TRUE))
  expect_error(itempool(ipdf), "Invalid id's. There are duplicated item id's.")
  ipdf <- data.frame(ID = sample(letters[1:3], n, TRUE),
                       a = rlnorm(n, 0, .3), b = rnorm(n),
                       content = sample(c("Geo", "Alg"), n, TRUE))
  expect_error(itempool(ipdf), "Invalid id's. There are duplicated item id's.")
  # ID or id as argument
  expect_error(itempool(data.frame(b = rnorm(n)),
                         ID = sample(letters[1:3], n, TRUE)),
               "Item ID cannot be duplicated.")
  expect_error(itempool(data.frame(b = rnorm(n)),
                         id = sample(letters[1:3], n, TRUE)),
               "Item ID cannot be duplicated.")

  ################# First element Numeric ##################################@###
  ip <- itempool(a = 1:2, b = 2:3, model = "2PL")
  ip <- itempool(a = 1:2, b = 2:3, c = c(.2, .3), model = "2PL")

  ################# se_parameters ##########################################@###
  # Test a simple list for se_parameters
  item_list <- generate_ip(n = 2, model = "2PL", output = "list")
  ip <- itempool(item_list,
                 se_parameters = list(list(a = .2, b = .3),
                                      list(a = .4, b = .5)))
  expect_equal(ip[[1]]@se_parameters$a, .2)
  expect_equal(ip[[1]]@se_parameters$b, .3)
  expect_equal(ip[[2]]@se_parameters$a, .4)
  expect_equal(ip[[2]]@se_parameters$b, .5)

  # -------------------------------------------------------------------------- #
  # Entering se_parametrs as data.frame
  ipdf <- as.data.frame(generate_ip(se_parameters = TRUE))
  ip <- itempool(ipdf[, c("a", "b", "c", "id")],
                 se_parameters = ipdf[, c("a_se", "b_se", "c_se")])
  i <- sample(1:nrow(ipdf), 1)
  expect_equal(ip[[i]]@se_parameters$a, ipdf$a_se[i])
  expect_equal(ip[[i]]@se_parameters$b, ipdf$b_se[i])
  expect_equal(ip[[i]]@se_parameters$c, ipdf$c_se[i])

})





###############################################################################@
############################# concatenation of 'Item' objects ##################
###############################################################################@

# Test concatenation function "c" of "Itempool" class
test_that("Test concatenation function 'c' of 'Itempool' class", {
  # Each element of the "Itempool" should be "Item" class
  item1 <- item(a = 1.12, b = -2.1, c = 0.28)
  item2 <- item(a = 2, b = 3.2, c = 0.21)
  item3 <- c(a = 1.2, b = 2.8, c = 0.12) # This is not 'Item' class
  item4 <- item(a = 1.12, b = -1.23, c = .2, id = "I-21")
  item5 <- item(a = 0.84, b = 2.23, c = .25, id = "I-22")
  # Create a new Itempool
  expect_is(object = c(item4, item5), class = "Itempool")
  # Creat a new Itempool object even the 'Item's don't have names
  expect_is(object = c(item1, item2), class = "Itempool")
  expect_is(object = c(item1, item2, item4), class = "Itempool")

  # -------------------------------------------------------------------------- #
  # mirt items
  item4 <- item(a = 1.12, b = -1.23, c = .2, id = "I-1")
  item5 <- item(a = 0.84, b = 2.23, c = .25, id = "I-2")
  item6 <- new("Item", model = "M2PL",
               parameters = list(a = c(2.1, 1.23, 1.3), d = -1.2, D = 1.7),
               id = "MI-1")
  item7 <- new("Item", model = "M2PL", parameters = list(
    a = c(1.1, 1.23), d = -1.2, D = 1.7), id = "MI-2")
  item8 <- new("Item", model = "M3PL", id = "MI-3", parameters = list(
    a = c(2.1, 1.23, 1.3), d = -1.2, c = .2, D = 1.7))
  item9 <- new("Item", model = "M3PL", content = "Algebra",
               parameters = list(a = c(2.1, 1.23, 1.3), d = -1.2, c = .2,
                                 D = 1.7))

  ip <- itempool(list(item4 ,item6, item7, item8))
  expect_is(object = ip, class = "Itempool")
  expect_is(object = c(item4 ,item6, item7, item8), class = "Itempool")
  expect_is(object = c(item4 ,item6, item7, item8, item9), class = "Itempool")

  # -------------------------------------------------------------------------- #
  # Concatenate Testlet and Item objects
  i1 <- item(b=rnorm(1))
  i2 <- item(b=rnorm(1))
  i3 <- item(b=rnorm(1))
  i4 <- item(b=rnorm(1))
  i5 <- item(b=rnorm(1), id = 'ii')
  i6 <- item(b=rnorm(1), id = 'ii')
  i7 <- item(b=rnorm(1), id = 'i7')
  t1 <- testlet(i3, i4)
  t2 <- testlet(i6, i7)
  ip <- c(i1, i2, t1)
  expect_is(ip, 'Itempool')
  expect_equal(length(ip), 3)
  # Combine two items and two testlets
  ip <- c(i1, i2, t1, t2)
  expect_is(ip, 'Itempool')
  expect_equal(length(ip), 4)

  t3 <- testlet(i3, i5)
  expect_error(c(i6, t3), "Invalid id's.")


  # -------------------------------------------------------------------------- #
  ### Errors ###
  # Give error when one of the item is not 'Item' class
  expect_error(c(item4, item3),
               regexp = "All of the elements should be 'Item' class.")
  expect_error(c(item4, item3, item6),
               regexp = "All of the elements should be 'Item' class.")
})

###############################################################################@
################### Subsetting 'Itempool' objects with "[" ####################
###############################################################################@
test_that("Test concatenation with [", {
  ip <- itempool(data.frame(a = runif(8, .5, 1.5), b = rnorm(8),
                               c = runif(8, 0,.3)), id = paste0("Item-",1:8),
                    content = c(rep("Algebra", 3), rep("Geometry", 2),
                                rep("Arithmetic", 3)))
  ip_new <- ip[1]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 1)
  expect_equal(ip_new[[1]], ip[[1]])

  ip_new <- ip[c(1:3)]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 3)
  expect_equal(ip_new[[2]], ip[[2]])

  ip_new <- ip[c(T, F, T, T, F, F, F, T)]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 4)
  expect_equal(ip_new[[2]], ip[[3]])

  ip_new <- ip[-2]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), length(ip) - 1)
  expect_equal(ip_new[[2]], ip[[3]])

  # If the logical vector is shorter than the vector being subsetted,
  # it will be recycled to be the same length.
  ip_new <- ip[c(T, F)]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 4)
  expect_equal(ip_new[[2]], ip[[3]])

  ip_new <- ip[c("Item-1", "Item-2")]

  # -------------------------------------------------------------------------- #
  # The order of numbers matter
  ip_new <- ip[c(6, 3, 7, 1)]
  expect_equal(ip_new[[1]], ip[[6]])
  expect_equal(ip_new[[2]], ip[[3]])
  expect_equal(ip_new[[3]], ip[[7]])
  expect_equal(ip_new[[4]], ip[[1]])

  # -------------------------------------------------------------------------- #
  # Missing indieces
  ip_new <- ip[c(7, NA, 2)]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 2)
  expect_equal(ip_new[[1]], ip[[7]])

  ip_new <- ip[c(NA, F, T, T, F, F, F, NA)]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 2)
  expect_equal(ip_new[[1]], ip[[3]])

  # -------------------------------------------------------------------------- #
  # Subsetting by IDs
  ip_new <- ip[c("Item-6", "Item-3", "Item-4")]
  expect_is(ip_new, "Itempool")
  expect_equal(length(ip_new), 3)
  expect_equal(ip_new[[1]], ip[[6]])
  expect_equal(ip_new[[2]], ip[[3]])
  expect_equal(ip_new[[3]], ip[[4]])

  # -------------------------------------------------------------------------- #
  # Errors
  # Cannot subset using an invalid item id.
  expect_error(ip_new <- ip[c("abc", "Item-3", "Item-4")], "Failed to subset")

  # -------------------------------------------------------------------------- #
  # misc is also transferred to the new Itempool object
  ip <- itempool(b = rnorm(5))
  ip$misc <- list(form_id = "F12")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  # The following tests are crucial. For example, in ability estimation using
  # maximum likelihood, the NA items are removed for standard error calculation
  # using ip[!is.na(resp)].
  expect_is(ip[1], "Itempool")
  expect_is(ip[c(1:3)], "Itempool")
  expect_is(ip[c(T, F, T, T, F)], "Itempool")


  # -------------------------------------------------------------------------- #
  # When subsetting an item pool, if the result is an empty item pool, return
  # NULL
  ip <- generate_ip(model = "2PL")
  expect_error(ip[ip$model == "GPC"],
               "The selection did not match any Item/Testlet object")

})

###############################################################################@
################### Subsetting 'Itempool' objects with "[[" ###################
###############################################################################@

test_that("Test concatenation with [[", {
  # Each element of the "Itempool" should be "Item" class
  ip <- itempool(data.frame(a = runif(8, .5, 1.5), b = rnorm(8),
                               c = runif(8, 0,.3)), id = paste0("Item-",1:8),
                    content = c(rep("Algebra", 3), rep("Geometry", 2),
                                rep("Arithmetic", 3)))

  expect_is(ip[[1]], "Item")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    item_list = list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  expect_is(ip[[1]], "Item")

  # -------------------------------------------------------------------------- #
  # When out of bounds, a custom error will be shown
  ip <- generate_ip(n = 10)
  expect_error(ip[[11]], paste0("Subscript out of bounds. Please use an index ",
                                "between 1 and 10."))
})

###############################################################################@
############################# [[<-  method (Itempool) #########################
###############################################################################@
test_that("Test Setting 'Itempool' objects with '[[<-'", {
  ip <- itempool(b = rnorm(5))
  item <- item(b = 1, a = 1.5)
  expect_equal(ip[[1]]$model, "1PL")
  ip[[1]] <- item
  expect_equal(ip[[1]]$model, "2PL")
  # An invalid assignment
  expect_error(ip[[1]] <- 12, "Invalid assignment.")
})

###############################################################################@
############################# $<-  method (Itempool) ##########################
###############################################################################@
test_that("Test $<- method (Itempool)", {
  # Setting 'misc'
  ip <- itempool(b = rnorm(5))
  expect_null(ip$misc)
  expect_null(ip@misc)
  ip$misc <- list(form_id = "Form-123")
  expect_equal(ip@misc$form_id, "Form-123")

  # -------------------------------------------------------------------------- #
  # Setting 'id'
  ip <- itempool(b = rnorm(2))
  expect_equal(ip$id[1], "Item-1")
  expect_equal(ip$id[2], "Item-2")
  ip$id <- c("a", "b")
  expect_equal(ip$id[1], "a")
  expect_equal(ip$id[2], "b")

  # -------------------------------------------------------------------------- #
  # Errors in setting 'id'
  ip <- itempool(b = rnorm(2))
  # id length and ip length should be the same
  expect_error(ip$id <- c("a", "b", "c"),
               "'id' should be a character vector with length equal to 2")
  expect_error(ip$id <- c(2, 1),
               "'id' should be a character vector with length equal to 2")
  expect_error(ip$id <- c("a", "a"),
               "'id' should not have any duplicated values")

  # -------------------------------------------------------------------------- #
  # Setting 'content'
  ip <- itempool(b = rnorm(2))
  expect_null(ip$content[1])
  expect_null(ip$content[2])
  ip$content <- c("a", "b")
  expect_equivalent(ip$content[1], "a")
  expect_equivalent(ip$content[2], "b")

  # -------------------------------------------------------------------------- #
  # One can set content to NULL
  ip <- itempool(b = rnorm(2), content = c("Algebra", "Geometry"))
  expect_equivalent(ip$content, c("Algebra", "Geometry"))
  ip$content <- NULL
  expect_null(ip$content)

  # -------------------------------------------------------------------------- #
  # Errors in setting 'content'
  ip <- itempool(b = rnorm(2))
  # id length and ip length should be the same
  expect_error(ip$content <- c("a", "b", "c"),
               "'content' should be a character vector with length equal to 2")
  expect_error(ip$content <- c(2, 1),
               "'content' should be a character vector with length equal to 2")

  # -------------------------------------------------------------------------- #
  # A single value can be assigned to all content
  ip <- itempool(b = rnorm(4))
  ip$content <- "Algebra"
  expect_equivalent(ip$content, rep("Algebra", 4))

  # -------------------------------------------------------------------------- #
  # Setting 'item_list'
  ip <- itempool(b = rnorm(2))
  item_list <- list(item(a = 1, b = 2), item(a = 2, b = 3, c = .2),
                    item(b = 1, model = 'Rasch'),
                    testlet(itempool(b = rnorm(2))))
  names(item_list) <- sapply(item_list, function(i) i$id)
  ip$item_list <- item_list
  expect_equal(length(ip), 4)
  expect_equal(ip[[3]]$model, 'Rasch')

  # -------------------------------------------------------------------------- #
  # Parameter values can be assigned to item pools through $<-
  ip <- itempool(b = rnorm(3))
  new_b <- rnorm(3)
  new_D <- rnorm(3, 6)
  ip$D <- new_D
  ip$b <- new_b
  expect_equivalent(ip$D, new_D)
  expect_equivalent(ip$b, new_b)

  # -------------------------------------------------------------------------- #
  # A single value can be assigned as a single parameters to all elements
  ip <- itempool(b = rnorm(3))
  ip$D <- 1.9
  ip$b <- 1
  expect_equivalent(ip$D, rep(1.9, 3))
  expect_equivalent(ip$b, rep(1, 3))

  # -------------------------------------------------------------------------- #
  # Errors in setting 'item_list'
  ip <- itempool(b = rnorm(2))
  # item_list should be a list of testlets or length and ip length should be the same
  expect_error(ip$item_list <- c("a", "b", "c"),
               "'item_list' should be a list of 'Item' or 'Testlet' objects.")
  expect_error(ip$item_list <- itempool(b = rnorm(3)),
               "'item_list' should be a list of 'Item' or 'Testlet' objects.")

  # -------------------------------------------------------------------------- #
  # Errors in setting unknown name
  ip <- itempool(b = rnorm(2))
  expect_error(ip$xyz <- "abc", "is not a valid name.")

  # -------------------------------------------------------------------------- #
  # Parameter values cannot be assigned to item pools with mixed models
  ip <- itempool(item(a = 1, b = 2), item(b = 1, model = "Rasch"))
  expect_error(ip$b <- 1, "is not a valid name.")
  expect_error(ip$D <- 1, "is not a valid name.")

})



###############################################################################@
############################# as.list (Itempool) ##############################
###############################################################################@

test_that("Test as.list( Itempool)", {
  # Each element of the "Itempool" should be "Item" class
  ip <- itempool(a = runif(4, .5, 1.5), b = rnorm(4), c = runif(4, 0,.3),
                 id = paste0("Item-",1:4), content = rep("Algebra", 4))
  expect_is(as.list(ip), "list")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    item_list = list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  expect_is(as.list(ip), "list")
})


###############################################################################@
############################# as.data.frame (Itempool) ########################
###############################################################################@

test_that("Test as.data.frame (Itempool)", {
  # Each element of the "Itempool" should be "Item" class
  ip <- generate_ip(model = "3PL", n = 4, content = rep("Algebra", 4))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_true(all(c("id", "model", "a", "b", "c", "D", "content") %in%
                    colnames(ipdf)))
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$id, ip$id)
  expect_equivalent(ipdf$a, ip$a)
  expect_equivalent(ipdf$b, ip$b)
  expect_equivalent(ipdf$c, ip$c)
  expect_equivalent(ipdf$D, ip$D)
  expect_equivalent(ipdf$model, ip$model)
  expect_equivalent(ipdf$content, ip$content)
  expect_is(ipdf$content, 'character')

  # -------------------------------------------------------------------------- #
  # Single item:
  ip <- itempool(data.frame(a = 7.4, b = rnorm(1), c = .3, D = 1.7),
                     id = paste0("Item-",1), content = rep("Algebra", 1))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equivalent(ipdf$id, ip$id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$model)
  expect_equivalent(ipdf$content, ip$content)

  # -------------------------------------------------------------------------- #
  # A mixture of Dichotomous items without a testlet.
  i3 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     id = paste0("i3-", 1:2))
  i4 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), id = paste0("i4", 1:2))
  i5 <- itempool(b = rnorm(3), id = paste0("i5-", 1:3))

  ip <- c(i3, i4, i5)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_true(all(c("id", "model", "a", "b", "c", "D") %in% colnames(ipdf)))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # A mixture of Dichotomous items without a testlet.
  ip <- generate_ip(model = c("2PL", "1PL", "4PL", "4PL", "1PL", "3PL", "Rasch"))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, 'data.frame')
  expect_true("model" %in% colnames(as.data.frame(ip)))
  expect_true(all(c("id", "model", "a", "b", "c", "d", "D") %in% colnames(ipdf)))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # A mixture of testlets and items but all of the standalone items and the
  # testlet items has the same model, All D's are equal
  n_t1 <- sample(3:9, 1)
  n_t2 <- sample(3:9, 1)
  n_t3 <- 1
  n_i1 <- sample(3:9, 1)
  n_i2 <- sample(3:9, 1)
  t1 <- testlet(itempool(
    a = rlnorm(n_t1, 0, .3), b = rnorm(n_t1), c = runif(n_t1, 0, .3),
    d = runif(n_t1, .95, 1), id = paste0("t1-i", 1:n_t1)), id = "t1")
  t2 <- testlet(itempool(
    a = rlnorm(n_t2, 0, .3), b = rnorm(n_t2), c = runif(n_t2, 0, .3),
    d = runif(n_t2, .95, 1), id = paste0("t2-i", 1:n_t2)), id = "t2")
  t3 <- testlet(itempool(
    a = rlnorm(n_t3, 0, .3), b = rnorm(n_t3), c = runif(n_t3, 0, .3),
    d = runif(n_t3, .95, 1), id = paste0("t3-i", 1:n_t3)), id = "t3")
  i1 <- generate_ip(n = n_i1, model = "4PL", id = paste0("i1-", 1:n_i1))
  i2 <- generate_ip(n = n_i2, model = "4PL", id = paste0("i2-", 1:n_i2))
  ip <- c(t1, i1, t2, i2, t3)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  # Testlet column is correctly specified
  expect_true(all(t1$id == ipdf[1:n_t1, "testlet"]))
  expect_true(is.na(ipdf[n_t1 + 1, "testlet"]))
  # Check parameters
  expect_equal(t2[[2]]$a, ipdf[n_t1 + n_i1 + 2, "a"])
  expect_equal(i1[[3]]$b, ipdf[n_t1 + 3, "b"])
  # Check model
  expect_true(all(ipdf$model == "4PL"))

  # -------------------------------------------------------------------------- #
  # A mixture of testlets and items but all of the standalone items and the
  # testlet items has the same model, D's are different
  n_t1 <- sample(3:6, 1)
  n_t2 <- sample(3:6, 1)
  n_t3 <- 1
  n_i1 <- sample(3:9, 1)
  n_i2 <- sample(3:9, 1)
  weird_D <- 92.8174
  t1 <- testlet(itempool(
    a = rlnorm(n_t1, 0, .3), b = rnorm(n_t1), c = runif(n_t1, 0, .3),
    d = runif(n_t1, .95, 1), id = paste0("t1-i", 1:n_t1), D = 1), id = "t1")
  t2 <- testlet(itempool(
    a = rlnorm(n_t2, 0, .3), b = rnorm(n_t2), c = runif(n_t2, 0, .3),
    d = runif(n_t2, .95, 1), id = paste0("t2-i", 1:n_t2), D = weird_D),
    id = "t2")
  t3 <- testlet(itempool(
    a = rlnorm(n_t3, 0, .3), b = rnorm(n_t3), c = runif(n_t3, 0, .3),
    d = runif(n_t3, .95, 1), id = paste0("t3-i", 1:n_t3)), id = "t3")
  i1 <- generate_ip(n = n_i1, model = "4PL", id = paste0("i1-", 1:n_i1))
  i2 <- generate_ip(n = n_i2, model = "4PL", id = paste0("i2-", 1:n_i2))
  ip <- c(t1, c(i1, i2), t2, t3)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equal(ipdf$D[n_t1 + n_i1 + n_i2 + 1], weird_D)
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # PCM with different number of thresholds
  ip <- itempool(
            list(
              new("Item", model = "PCM", id = "PCM-1", content = "Stress",
                  parameters = list(b = c(-0.837, 0.529))),
              new("Item", model = "PCM", id = "PCM-2", content = "Stress",
                  parameters = list(b = c(-0.837, 0.529, 1.2)))
            ))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equal(colnames(ipdf), c("id", "model", "b1", "b2", "b3", "content"))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # GPCM with different number of thresholds
  ip <- itempool(
    list(new("Item", model = "GRM", id = "GRM-1", content = "Stress",
             parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382,
                                                1.649), D = 1.702)),
         new("Item", model = "GRM", id = "GRM-2", content = "Depression",
             parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567),
                               D = 1.702))
            ))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equal(colnames(ipdf), c("id", "model", "a", "b1", "b2", "b3",
                                 "b4", "D", "content"))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # mirt with different number of dimensions
  ip <- itempool(list(
    new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
        parameters = list(a = c(runif(2,1,2)), d = rnorm(1), D = 1.7)),
    new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
        parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
    new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
        parameters = list(a = c(runif(4,1,2)), d = rnorm(1), D = 1.7))
    ))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equal(colnames(ipdf), c("id", "model", "a1", "a2", "a3", "a4",
                                 "d", "D", "content"))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # GRM item with only two categories and a 3PL item
  # TODO: Check whether the following "Item-1" with only 2 categories converted to
  # the data frame nicely. It the threshold parameter is not labelled as "b1"
  # but as "b"
  ip <- c(
    generate_item(model = "GRM", n_categories = 2),
    generate_item(model = "GRM", n_categories = 3),
    generate_item(model = "3PL", n = 4)
    )
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")

  # -------------------------------------------------------------------------- #
  # A mixture of Polytomous items without a testlet.
  ip <- itempool(list(
    new("Item", model = "GRM", id = "GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "GRM-2", content = "Stress",
        parameters = list(a = 1.41, b = c(-0.837, 0.529), D = 1.702)),
    new("Item", model = "GRM", id = "GRM-3", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702)),
    new("Item", model = "GPCM", id = "GPCM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, 0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GPCM", id = "GPCM-2", content = "Stress",
        parameters = list(a = 1.41, b = c(-0.837, 0.529), D = 1.702)),
    new("Item", model = "GPCM", id = "GPCM-3", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702)),
    new("Item", model = "PCM", id = "PCM-1", content = "Stress",
        parameters = list(b = c(-0.837, 0.529))),
    new("Item", model = "PCM", id = "PCM-2", content = "Stress",
        parameters = list(b = c(-0.837, 0.529, 1.2)))
            ))
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_true(all(c("id", "model", "a", "b1", "b2", "b3", "b4", "D",
                    "content") %in% colnames(ipdf)))
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(ipdf$content, ip$resp_content)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # Polytomous items with a testlet.
  t1 <- testlet(itempool(list(
    new("Item", model = "GRM", id = "t1GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "t1GRM-2",
        parameters = list(a = 1.41, b = c(-0.837, 0.529), D = 1.702)))),
    id = "t1")
  i2 <- itempool(list(
    new("Item", model = "GRM", id = "GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "GRM-2", content = "Depression",
        parameters = list(a = 1.41, b = c(-0.837, 0.529), D = 1.702))))

  ip <- c(t1, i2)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_true(all(c("id", "testlet", "model", "a", "b1", "b2", "b3", "b4", "D",
                    "content") %in% colnames(ipdf)))
  expect_equivalent(ipdf$content, ip$resp_content)
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # A mixture of Dichotmous and Polytomous items without a testlet.
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1), id = paste0("i1-", 1:4))
  i2 <- itempool(list(
    new("Item", model = "GRM", id = "GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "GRM-2", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702))
            ))
  i3 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     id = paste0("i3-", 1:2))
  i4 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), id = paste0("i4-",
                                                                     1:2))
  i5 <- itempool(b = rnorm(3), id = paste0("i5-", 1:3))

  ip <- c(i1, i2, i3, i4, i5)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equivalent(ipdf$content, ip$resp_content)
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # A mixture of items: Dichotomous and Polytomous items scattered among
  # testlets
  t1 <- generate_testlet(n = 3, item_models = "3PL", id = "t1")
  i1 <- generate_ip(model = "3PL", n = 4, id = paste("i1-", 1:4))
  t2 <- generate_testlet(item_models = c("GRM", "3PL", "3PL", "3PL", "GRM"),
                         id = "t2", item_id_preamble = "t2-")
  i2 <- c(generate_item(model = "GRM", n_categories = 3, id = "grm-1",
                        content = "Stress"),
          generate_item(model = "GRM", n_categories = 4, id = "grm-2",
                        content = "Depression"))
  t3 <- generate_testlet(item_models = c("GRM", "2PL", "3PL", "Rasch", "GPCM"),
                         id = "t3", item_id_preamble = "t3-")
  ip <- c(t1, i1, t2, i2, t3)
  ipdf <- as.data.frame(ip)
  expect_is(ipdf, "data.frame")
  expect_equivalent(ipdf$id, ip$resp_id)
  expect_equivalent(rownames(ipdf), paste0(1:ip$n$items))
  expect_equivalent(ipdf$model, ip$resp_model)

  # -------------------------------------------------------------------------- #
  # If a misc field has one printable value, it should be added as a column.
  ip <- c(generate_item(misc = list(form = "Form-1", key = "A")),
          generate_item(misc = list(form = "Form-1", key = "C", level = "C-")),
          # Add an item without any misc field
          generate_item(),
          # add an item with a misc field length more than 1
          generate_item(misc = list(form = "Form-1", key = "C",
                                    enemies = c("i1", "i2"))),
          generate_item(misc = list(form = "Form-1", key = "B", level = 12))
          )
  ipdf <- as.data.frame(ip)
  # a misc field with length more than 1 is not included in the data frame
  expect_false("enemies" %in% colnames(ipdf))
  expect_true(all(c("form", "key", "level") %in% colnames(ipdf)))
  expect_equal(ipdf$form[5], ip[[5]]$misc$form)
  expect_equal(ipdf$level[2], ip[[2]]$misc$level)
  # numeric misc field converted to character
  expect_equal(ipdf$level[5], as.character(ip[[5]]$misc$level))
  expect_true(is.na(ipdf$form[3]))

  # -------------------------------------------------------------------------- #
  # misc field when item pool has testlets
  i1 <- item(b = 1, id = "I-1", content = "bec",
             misc = list(sympson_hetter_k = 0.1,
                         form = "A1"))
  t1 <- testlet(itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                         misc = list(list(sympson_hetter_k = .2, form = "b3"),
                                     list(sympson_hetter_k = .3))),
                   id = "t1")
  ip <- c(i1, t1)
  ipdf <- as.data.frame(ip)
  expect_equal(ipdf$form[1], i1@misc$form)
  expect_equal(ipdf$sympson_hetter_k[1], i1@misc$sympson_hetter_k)
  expect_equal(ipdf$form[2], t1[[1]]@misc$form)
  expect_equal(ipdf$sympson_hetter_k[2], t1[[1]]@misc$sympson_hetter_k)
  expect_equal(ipdf$form[3], as.character(NA))
  expect_equal(ipdf$sympson_hetter_k[3], t1[[2]]@misc$sympson_hetter_k)

  # -------------------------------------------------------------------------- #
  # "2PL" item pool with standard errors
  ip <- generate_ip(n = 10, model= "2PL", se_parameters = TRUE)
  ipdf <- as.data.frame(ip)
  expect_true("a_se" %in% colnames(ipdf))
  expect_true("b_se" %in% colnames(ipdf))

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # - se_parameters                                                          - #
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # All same dichotomous models where all items have se_parameters
  ip <- generate_ip(model = "2PL", se_parameters = TRUE)
  ipdf <- as.data.frame(ip)
  expect_true("a_se" %in% colnames(ipdf))
  expect_true("b_se" %in% colnames(ipdf))
  i <- sample(1:length(ip), 1)
  expect_equal(ipdf$a_se[i], ip[[i]]@se_parameters$a)
  expect_equal(ipdf$b_se[i], ip[[i]]@se_parameters$b)
  # When 'include_se = FALSE' se_parameters are not included in data frame
  ipdf <- as.data.frame(ip, include_se = FALSE)
  expect_false("a_se" %in% colnames(ipdf))
  expect_false("b_se" %in% colnames(ipdf))

  # -------------------------------------------------------------------------- #
  # All same polytomous models with same number of categories where all items
  # have se_parameters
  ip <- generate_ip(model = "GPCM", se_parameters = TRUE)
  ipdf <- as.data.frame(ip)
  expect_true("a_se" %in% colnames(ipdf))
  expect_true("b1_se" %in% colnames(ipdf))
  i <- sample(1:length(ip), 1)
  expect_equal(ipdf$a_se[i], ip[[i]]@se_parameters$a)
  expect_equal(ipdf$b1_se[i], ip[[i]]@se_parameters$b[1])

  # -------------------------------------------------------------------------- #
  # All same polytomous models but with different number of categories where
  # all items have se_parameters
  n_items <- sample(5:10, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  ip <- generate_ip(model = "GPCM", n = n_items, n_categories = n_categories,
                    se_parameters = TRUE)
  ipdf <- as.data.frame(ip)
  expect_true("a_se" %in% colnames(ipdf))
  expect_true("b1_se" %in% colnames(ipdf))
  i <- sample(1:length(ip), 1)
  expect_equal(ipdf$a_se[i], ip[[i]]@se_parameters$a)
  expect_equal(ipdf$b1_se[i], ip[[i]]@se_parameters$b[1])

  # -------------------------------------------------------------------------- #
  # Mixture of polytomous and dichotomous items with different number of
  # categories where all items have se_parameters
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = TRUE)
  ipdf <- as.data.frame(ip)
  expect_true("a_se" %in% colnames(ipdf))
  expect_true("b_se" %in% colnames(ipdf))
  i <- sample(1:length(ip), 1)
  expect_equal(ipdf$a_se[i], ip[[i]]@se_parameters$a)
  expect_equal(ipdf$b_se[2], ip[[2]]@se_parameters$b)
  expect_equal(ipdf$b2_se[3], ip[[3]]@se_parameters$b[2])
  expect_equal(ipdf$b_se[4], ip[[4]]@se_parameters$b)
  expect_equal(ipdf$d2_se[4], ip[[4]]@se_parameters$d[2])

  # -------------------------------------------------------------------------- #
  # An item pool without any se_parameters should not return any se_parameters
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = NULL)
  ipdf <- as.data.frame(ip)
  expect_false("a_se" %in% colnames(ipdf))
  expect_false("b_se" %in% colnames(ipdf))
})


###############################################################################@
############################# is.Itempool ######################################
###############################################################################@

test_that("Test is.Itempool", {
  # Each element of the "Itempool" should be "Item" class
  item1 <- item(a = 1.12, b = -2.1, c = 0.28)
  item2 <- item(a = 2, b = 3.2, c = 0.21)
  expect_true(is.Itempool(c(item1, item2)))

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
    ))
  expect_true(is.Itempool(ip))
})


###############################################################################@
############################# length (Itempool) ###############################
###############################################################################@

test_that("Test length (Itempool)", {
  n <- sample(10:100,1)
  ip <- itempool(a = runif(n, .5, 1.5), b = rnorm(n),
                 c = runif(n, 0,.3), id = paste0("Item-",1:n))
  expect_true(length(ip) == n)

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))))
  expect_equivalent(length(ip), 5)
})

###############################################################################@
################### $ (Itempool)  #############################################
###############################################################################@

test_that("Test $ operator", {
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # Overall Tests
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # Each element of the "Itempool" should be "Item" class
  n <- sample(10:20, 1)
  ip <- itempool(a = runif(n, .5, 1.5), b = rnorm(n), c = runif(n, 0, .3),
                 id = paste0("Item_",1:n), content = rep("Algebra", n))
  ip_list <- ip$item_list
  i <- sample(1:n, 1)
  expect_equivalent(ip$model[i], "3PL")
  expect_equal(names(ip$model)[i], ip$id[i]) # vectors should be named
  expect_equivalent(ip$content[i], "Algebra")
  expect_equal(names(ip$content)[i], ip$id[i]) # vectors should be named
  expect_equal(ip$id[i], paste0("Item_", i))
  expect_is(ip$parameters, "matrix")
  expect_equal(ip$parameters[3,1], ip[[3]]@parameters$a)
  expect_equivalent(ip$parameters[3, "a"] , ip[[3]]@parameters$a)
  # Extract single parameter
  expect_equivalent(ip$a[i], ip_list[[i]]@parameters$a)
  expect_equivalent(ip$b[i], ip_list[[i]]@parameters$b)
  expect_equivalent(ip$D[i], ip_list[[i]]@parameters$D)
  expect_equal(names(ip$a)[i], ip$id[i]) # vectors should be named
  expect_equal(names(ip$b)[i], ip$id[i]) # vectors should be named
  expect_equal(names(ip$D)[i], ip$id[i]) # vectors should be named


  # Benchmarking
  # n <- sample(1000:2000, 1)
  # ip <- generate_ip(n = n, id = paste0("Item_",1:n),
  #                   content = rep("Algebra", n))
  #
  # library(Rcpp)
  # cppFunction('
  # Rcpp::Nullable<Rcpp::StringVector> gsic(Rcpp::S4 ip, std::string slotName)
  # {
  #   Rcpp::S4 tempS4;
  #   List item_list = as<List>(ip.slot("item_list"));
  #   int noItem = item_list.size();
  #   Rcpp::StringVector output(noItem);
  #   int count_empty_str = 0;
  #   for (int i = 0; i < noItem; i++)
  #   {
  #     tempS4 = as<Rcpp::S4> (item_list(i));
  #     if (Rf_isNull(tempS4.slot(slotName))) {
  #       output[i] = StringVector::get_na();
  #       count_empty_str++;
  #     } else output[i] = as<std::string>(tempS4.slot(slotName));
  #   }
  #   if (count_empty_str == noItem)
  #   {
  #     return R_NilValue;
  #   }
  #   return(output);
  # }
  # ')
  #
  # microbenchmark::microbenchmark(
  #   cpp = gsic(ip, "model"),
  #   R = sapply(ip@item_list, slot, "model"),
  #   times = 1e2
  #   )


  # -------------------------------------------------------------------------- #
  # Test with various item types
  ip <- itempool(list(
    item(b = 1, D = 1, id = "I-1", content = "bec"),
    item(a = 1, b = 2, id = "I-2", content = "cab"),
    item(a = 1, b = 2, c = .12, id = "I-3", content = "dac")))
  expect_equivalent(ip$model[1], "1PL")
  expect_equivalent(ip$content[2], "cab")
  expect_equal(ip$id[3], "I-3")
  expect_is(ip$parameters, "list")
  expect_equal(ip$parameters[[3]]$c, 0.12)

  # -------------------------------------------------------------------------- #
  # Test 1PL
  ip <- itempool(data.frame(b = rnorm(10)))
  expect_equivalent(ip$model[1], "Rasch")
  expect_true(is.null(ip$content))
  expect_equal(ip$id[3], "Item-3")
  expect_is(ip$parameters, "matrix")
  expect_is(ip$parameters[,"b"], "numeric")

  # -------------------------------------------------------------------------- #
  # For GPCM with different number of categories, "$" still can return parameter
  ip <- itempool(lapply(2:5, function(x)
    generate_item(model = "GPCM", n_categories = x)))
  a <- sapply(ip$item_list, function(x) x$a)
  expect_equivalent(a, ip$a)
  expect_equal(names(ip$a), ip$id)

  # -------------------------------------------------------------------------- #
  # One should not extract an unknown parameter
  ip <- itempool(b = rnorm(4))
  expect_null(ip$kkk)
  expect_null(ip$a)

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))))
  expect_equivalent(ip$model[1], "M2PL")
  expect_equal(ip$id[3], "Item 3")
  expect_is(ip$parameters, "matrix")
  expect_equal(ip$parameters[3, 2], ip@item_list[[3]]@parameters$a[2])
  expect_equivalent(ip$parameters[,'d'][4],
                    ip@item_list[[4]]@parameters$d)
  expect_equivalent(ip$content[4], "Algebra")
  expect_equivalent(which(ip$content %in% c("Algebra"))[2], 4)
  expect_equivalent(which(ip$content %in% c("Algebra", "Arithmetic"))[3], 5)
  expect_equivalent(which(ip$content %in% c("Algebra")),
                    which(ip$content %in% c("Algebra", "Reading")))


  # -------------------------------------------------------------------------- #
  # Graded Response Model (GRM)
  ip <- itempool(
    list(
      new("Item", model = "GRM", id = "Item 1", content = "Stress",
          parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382,  1.649),
                            D = 1.702)),
      new("Item", model = "GRM", id = "Item 2", content = "Depression",
          parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567,  0.646),
                            D = 1.702))))
  expect_equivalent(ip$model[1], "GRM")
  expect_equal(names(ip$model)[1], "Item 1")
  expect_equal(ip$id[2], "Item 2")
  expect_is(ip$parameters, "matrix")

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $content
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ip <- itempool(data.frame(a = runif(8, .5, 1.5), b = rnorm(8),
                               c = runif(8, 0,.3)), id = paste0("MyItem-",1:8),
                    content = c(rep("Algebra", 3), rep("Geometry", 2),
                                rep("Arithmetic", 3)))
  expect_equivalent(ip$content[4], "Geometry")
  expect_equivalent(ip$id[3], "MyItem-3")
  expect_equivalent(which(ip$content == "Algebra")[2], 2)
  expect_equivalent(which(ip$content %in% c("Algebra", "Arithmetic"))[4], 6)
  expect_equivalent(which(ip$content %in% "Algebra"),
                    which(ip$content %in% c("Algebra", "Reading")))

  # -------------------------------------------------------------------------- #
  # Test content where some content is missing
  ip <- itempool(item(b = 0, content = "A"),
                  item(b = 1, id = "i2"),
                  item(b = 2, content = "C"))
  expect_true(is.na(ip$content[2]))

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $se_parameters
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # All same dichotomous models where all items have se_parameters
  ip <- generate_ip(model = "2PL", se_parameters = TRUE)
  se_pars <- ip$se_parameters
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$b[i], ip[[i]]@se_parameters$b)

  # -------------------------------------------------------------------------- #
  # All same dichotomous models where all items have se_parameters except one
  # is NULL
  ip <- generate_ip(model = "2PL", se_parameters = TRUE)
  ip[[1]]@se_parameters <- NULL
  se_pars <- ip$se_parameters
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$b[i], ip[[i]]@se_parameters$b)
  expect_true(is.na(se_pars$a[1]))
  expect_true(is.na(se_pars$b[1]))

  # -------------------------------------------------------------------------- #
  # All same polytomous models with same number of categories where all items
  # have se_parameters
  ip <- generate_ip(model = "GPCM", se_parameters = TRUE)
  se_pars <- ip$se_parameters
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$b1[i], ip[[i]]@se_parameters$b[1])


  # -------------------------------------------------------------------------- #
  # All same polytomous models with same number of categories where all items
  # have se_parameters except one is NULL
  ip <- generate_ip(model = "GPCM", se_parameters = TRUE)
  ip[[1]]@se_parameters <- NULL
  se_pars <- ip$se_parameters
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$b1[i], ip[[i]]@se_parameters$b[1])
  expect_true(is.na(se_pars$a[1]))
  expect_true(is.na(se_pars$b1[1]))

  # -------------------------------------------------------------------------- #
  # All same polytomous models but with different number of categories where
  # all items have se_parameters
  n_items <- sample(5:10, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  ip <- generate_ip(model = "GPCM", n = n_items, n_categories = n_categories,
                    se_parameters = TRUE)
  se_pars <- ip$se_parameters
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(1:length(ip), 1)
  expect_equal(se_pars$b1[i], ip[[i]]@se_parameters$b[1])

  # -------------------------------------------------------------------------- #
  # All same polytomous models but with different number of categories where
  # all items have se_parameters except one is NULL
  n_items <- sample(5:10, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  ip <- generate_ip(model = "GPCM", n = n_items, n_categories = n_categories,
                    se_parameters = TRUE)
  ip[[1]]@se_parameters <- NULL
  se_pars <- ip$se_parameters
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$a[i], ip[[i]]@se_parameters$a)
  i <- sample(2:length(ip), 1)
  expect_equal(se_pars$b1[i], ip[[i]]@se_parameters$b[1])
  expect_true(is.na(se_pars$a[1]))
  expect_true(is.na(se_pars$b1[1]))

  # -------------------------------------------------------------------------- #
  # Mixture of polytomous and dichotomous items with different number of
  # categories where all items have se_parameters
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = TRUE)
  se_pars <- ip$se_parameters
  expect_equal(se_pars$a[1], ip[[1]]@se_parameters$a)
  expect_equal(se_pars$b[1], ip[[1]]@se_parameters$b)
  expect_equal(se_pars$c[2], ip[[2]]@se_parameters$c)
  expect_equal(se_pars$a[3], ip[[3]]@se_parameters$a)
  expect_equal(se_pars$b[4], ip[[4]]@se_parameters$b)
  expect_equal(se_pars$d2[4], ip[[4]]@se_parameters$d[2])

  # -------------------------------------------------------------------------- #
  # Mixture of polytomous and dichotomous items with different number of
  # categories where all items have se_parameters except one is NULL
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = TRUE)
  ip[[length(ip)]]@se_parameters <- NULL
  se_pars <- ip$se_parameters

  expect_equal(se_pars$a[1], ip[[1]]@se_parameters$a)
  expect_equal(se_pars$b[1], ip[[1]]@se_parameters$b)
  expect_equal(se_pars$c[2], ip[[2]]@se_parameters$c)
  expect_equal(se_pars$a[3], ip[[3]]@se_parameters$a)
  expect_equal(se_pars$b[4], ip[[4]]@se_parameters$b)
  expect_equal(se_pars$d2[4], ip[[4]]@se_parameters$d[2])
  expect_true(all(is.na(unlist(se_pars[length(ip), ]))))

  # -------------------------------------------------------------------------- #
  # Mixture of polytomous and dichotomous items with different number of
  # categories where all items have se_parameters except one is NULL and other
  # item has one parameter SE which is NA
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = TRUE)
  ip[[length(ip)]]@se_parameters <- NULL
  ip[[length(ip)-1]]@se_parameters$a <- NA
  se_pars <- ip$se_parameters
  expect_equal(se_pars$a[1], ip[[1]]@se_parameters$a)
  expect_equal(se_pars$b[1], ip[[1]]@se_parameters$b)
  expect_equal(se_pars$c[2], ip[[2]]@se_parameters$c)
  expect_equal(se_pars$a[3], ip[[3]]@se_parameters$a)
  expect_equal(se_pars$b[4], ip[[4]]@se_parameters$b)
  expect_equal(se_pars$d2[4], ip[[4]]@se_parameters$d[2])
  expect_true(all(is.na(unlist(se_pars[length(ip), ]))))
  expect_true(is.na(se_pars$a[length(ip)-1]))

  # -------------------------------------------------------------------------- #
  # An item pool without any se_parameters should return NULL
  n_items <- sample(9:19, 1)
  n_categories <- sample(3:7, n_items, replace = TRUE)
  models <- c("2PL", "3PL", "GPCM", "GPCM2",
              sample(c("2PL", "3PL", "GPCM", "GPCM2"), n_items - 4,
                     replace = TRUE))
  ip <- generate_ip(model = models, n_categories = n_categories,
                    se_parameters = NULL)
  expect_null(ip$se_parameters)


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $a, $b, $c, $D
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # The parameter order should not be important
  a1 <- rlnorm(1, 0, .3)
  b1 <- rnorm(1)
  c1 <- runif(1, 0, .3)
  D1 <- 1.7
  a2 <- rlnorm(1, 0, .3)
  b2 <- rnorm(1)
  c2 <- runif(1, 0, .3)
  D2 <- 1

  i1 <- new("Item", parameters = sample(list(b = b1, a = a1, D = D1, c = c1)),
            model = "3PL", misc = list(form = "F1"))
  i2 <- new("Item", parameters = sample(list(c = c2, D = D2, a = a2, b = b2)),
            model = "3PL", misc = list(form = "F2"))
  ip <- c(i1, i2)
  expect_equal(i1@parameters$a, a1)
  expect_equal(i1@parameters[[1]], a1)
  expect_equal(i2@parameters$a, a2)
  expect_equal(i1@parameters$b, b1)
  expect_equal(i1@parameters[[2]], b1)
  expect_equal(i2@parameters$b, b2)
  expect_equal(i1@parameters$c, c1)
  expect_equal(i1@parameters[[3]], c1)
  expect_equal(i2@parameters$c, c2)
  expect_equal(i1@parameters$D, D1)
  expect_equal(i1@parameters[[4]], D1)
  expect_equal(i2@parameters$D, D2)

  expect_equivalent(ip$a, c(a1, a2))
  expect_equivalent(ip$b, c(b1, b2))
  expect_equivalent(ip$c, c(c1, c2))
  expect_equivalent(ip$D, c(D1, D2))


  a1 <- rlnorm(1, 0, .3)
  b1 <- rnorm(1)
  c1 <- runif(1, 0, .3)
  D1 <- 1.7
  i1 <- item(b = b1, a = a1, D = D1, c = c1)
  expect_equal(i1@parameters$a, a1)
  expect_equal(i1@parameters$b, b1)
  expect_equal(i1@parameters$c, c1)
  expect_equal(i1@parameters$D, D1)

  i1 <- item(parameters = sample(list(b = b1, a = a1, D = D1, c = c1)))
  expect_equal(i1@parameters$a, a1)
  expect_equal(i1@parameters$b, b1)
  expect_equal(i1@parameters$c, c1)
  expect_equal(i1@parameters$D, D1)

  # -------------------------------------------------------------------------- #
  # Pulling parameters when there are multiple models:
  ip <- generate_ip(model = c("2PL", "GPCM", "GPCM2"))
  expect_false(is.null(ip$a))
  expect_null(ip$c)
  expect_equivalent(ip$a[1], ip[[1]]@parameters$a)
  expect_equivalent(ip$a[2], ip[[2]]@parameters$a)
  expect_equivalent(ip$a[3], ip[[3]]@parameters$a)
  expect_equivalent(ip$b[1], ip[[1]]@parameters$b)
  expect_true(is.na(ip$b[2]))
  expect_equivalent(ip$b[3], ip[[3]]@parameters$b)
  expect_true(is.na(ip$d1[1]))
  expect_true(is.na(ip$d1[2]))
  expect_equivalent(ip$d1[3], ip[[3]]@parameters$d[1])


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $resp_ip, $resp_content, $resp_model, $resp_item_list
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  i1 <- item(b = 1, id = "i1", content = "A")
  i2 <- item(a = 1, b = 1, id = "i2", content = "B")
  i3 <- item(a = 1, b = 1, c = 1, id = "i3", content = "C")
  i4 <- item(a = 1, b = c(1, 2), id = "i4", content = "D", model = "GRM")
  i5 <- item(a = 1, b = c(1, 2), id = "i5", content = "E", model = "GPCM")
  t1 <- testlet(i2, i5, id = "testlet1", content = "TZM")
  t2 <- generate_testlet(n = 3, item_id_preamble = "t2-", id = "t2")
  ip <- c(i1, i3, t1, i4, t2)
  # There is a warning when printing this item pool:
  expect_output(show(ip), "Itempool")
  expect_equal(ip$id, c("i1", "i3", "testlet1", "i4", "t2"))
  expect_equal(ip$resp_id, c("i1", "i3", "i2", "i5", "i4",
                             paste0("t2-Item-", 1:3)))
  expect_equivalent(ip$content, c("A", "C", "TZM", "D", NA))
  expect_equivalent(ip$resp_content, c("A", "C", "B", "E", "D", NA, NA, NA))
  expect_equal(names(ip$resp_content), c("i1", "i3", "i2", "i5", "i4",
                                         paste0("t2-Item-", 1:3)))
  expect_equivalent(ip$model, c("Rasch", "3PL", "BTM", "GRM", "BTM"))
  expect_equivalent(ip$resp_model, c("Rasch", "3PL", "2PL", "GPCM", "GRM",
                                     "3PL", "3PL", "3PL"))

  # $resp_item_list: A list of standalone items
  expect_true("testlet1" %in% ip$id)
  sa_ip_list <- ip$resp_item_list
  expect_false("testlet1" %in% names(sa_ip_list))
  expect_false("i2" %in% ip$id)
  expect_true("i2" %in% names(sa_ip_list))
  expect_equal(i4, sa_ip_list[[5]])
  expect_equal(t2[[2]], sa_ip_list[[7]])

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $n
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  expect_is(ip$n, "list")
  expect_equal(ip$n$testlets, 2)
  expect_equal(ip$n$elements, 5)
  expect_equal(ip$n$items, 8)


  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $items
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  ip <- generate_ip(n = 7)
  expect_is(ip$items, "list")
  expect_equal(length(ip$items), 7)
  expect_equal(names(ip$items), ip$id)

  # With testlets
  ip <- c(generate_ip(n = 2), generate_testlet(n = 4))
  ip_list <- ip$items
  expect_is(ip_list, "list")
  expect_equal(length(ip_list), 6)
  expect_equal(ip_list[[4]], ip[[3]]@item_list[[2]])
  expect_equal(names(ip_list), ip$resp_id)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $max_Score
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  ip <- c(generate_ip(model = c("GPCM", "GRM", "3PL"),
                      n_categories = c(7, 5, 2)),
          generate_testlet(n = 4, item_models = c("GRM", "2PL", "GRM", "GPCM"),
                           n_categories = 3))
  expect_equal(ip$max_score, 18)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $resp_misc
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  i1 <- item(b = 1, id = "I-1", content = "bec",
             misc = list(sympson_hetter_k = 1,
                         form = "A1"))
  t1 <- testlet(itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                         misc = list(list(sympson_hetter_k = .8, form = "B1"),
                                     list(sympson_hetter_k = .9))),
                   id = "t1")
  t2 <- generate_testlet(n = 4, id = "t2")
  t2[[2]]@misc <- list(form = "B2")
  t3 <- generate_testlet(n = 3, id = "t3")
  i2 <- generate_item(misc = list(form = "A1"))
  i3 <- generate_item()
  ip <- c(i1, t1, i2, t2, i3, t3)
  misc <- ip$resp_misc
  expect_equal(i1$misc, misc[[1]])
  expect_equal(t1@item_list$misc[[1]], misc[[2]])
  expect_equal(t1@item_list$misc[[2]], misc[[3]])
  expect_equal(t2@item_list$misc[[2]], misc[[6]])
  expect_equal(t2@item_list$misc[[3]], misc[[7]])
  expect_equal(i2$misc, misc[[4]])
  expect_equal(i3$misc, misc[[9]])
  expect_equal(t3@item_list$misc[[3]], misc[[ip$n$items]])

  # -------------------------------------------------------------------------- #
  # $resp_misc returns NULL when no items or testlet has misc field
  ip <- c(generate_ip(), generate_testlet())
  expect_null(ip$resp_misc)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $misc
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  # misc field should return a list
  ip <- c(generate_item(misc = list(form = "F1")),
          generate_item(misc = list(form = "F2")))
  expect_is(ip$misc, "list")
  expect_equal(length(ip$misc), 2)
  expect_equal(names(ip$misc), ip$id)
  expect_true(all(sapply(ip$misc, names) == "form"))

  # -------------------------------------------------------------------------- #
  # misc: when no misc field function should return NULL
  ip <- generate_ip()
  expect_null(ip$misc)

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  # $resp_max_score
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #

  ip <- generate_ip(n = 7)
  ip_max_score <- ip$resp_max_score
  expect_equivalent(ip_max_score, rep(1, 7))
  expect_equal(names(ip_max_score), ip$id)

  # With testlets
  ip <- c(generate_ip(n = 2),
          generate_testlet(n = 4, item_models = c("GRM", "2PL", "GRM"),
                           n_categories = 3))
  ip_max_score <- ip$resp_max_score
  expect_equivalent(ip_max_score, c(1, 1, 2, 1, 2, 2))
  expect_equal(names(ip_max_score), ip$resp_id)

  # With testlets and polytomous items
  ip <- c(generate_ip(n = 2),
          generate_testlet(
            id = "t1",
            item_models = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
            n_categories = c(2, 3, 6, 7, 2)),
          generate_testlet(n = 3, id = "t2"))
  ip_max_score <- ip$resp_max_score
  expect_equal(names(ip_max_score), ip$resp_id)
  expect_equivalent(ip_max_score, c(1, 1, 1, 2, 5, 6, 1, 1, 1, 1))
})

###############################################################################@
############################# print.Itempool ##################################
###############################################################################@

test_that("Test print.Itempool", {
  # Each element of the "Itempool" should be "Item" class
  ip <- generate_ip(n = 4, id = paste0("Item-", 1:4),
                    content = rep("Algebra", 4))
  expect_output(print(ip), "An object of class 'Itempool'")
  expect_true(grepl("An object of class 'Itempool'", capture_output(print(ip))))
  expect_output(print(ip), "Model of items: (.*)3PL(.*)D = ")
  expect_output(print(ip), "D = 1")


  ip <- itempool(data.frame(a = runif(4, .5, 1.5), b = 1:4, D = 1:4),
                    id = paste0("Item-",1:4),
                    content = rep("Algebra", 4))
  expect_output(print(ip), "Model of items: ")
  expect_output(print(ip), "a b D")
  expect_output(print(ip), "Content = Algebra")

  # -------------------------------------------------------------------------- #
  # Test with various item types
  ip <- itempool(list(item(b = 1, id = "I-1", content = "bec"),
                      item(a = 1, b = 2, id = "I-2", content = "cab"),
                      item(a = 1, b = 2, c = .12, id = "I-3", content = "dac")))
  expect_output(print(ip), "An object of class 'Itempool'")

  ip <- itempool(data.frame(b = rnorm(10)))
  expect_output(print(ip), "An object of class 'Itempool'")

  # -------------------------------------------------------------------------- #
  # Test Itempool with a misc field printed correctly
  i1 <- item(b = 1, id = "I-1", content = "bec",
             misc = list(sympson_hetter_k = 1))
  t1 <- testlet(itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                         misc = list(list(sympson_hetter_k = .8),
                                     list(sympson_hetter_k = .9))),
                   id = "t1")
  ip <- c(i1, t1)
  expect_output(print(ip), "An object of class 'Itempool'")

  # -------------------------------------------------------------------------- #
  # mirt
  ip <- itempool(
    list(
      new("Item", model = "M2PL", id = "Item 1", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 2", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 3", content = "Geometry",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 4", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7)),
      new("Item", model = "M2PL", id = "Item 5", content = "Algebra",
          parameters = list(a = c(runif(3,1,2)), d = rnorm(1), D = 1.7))
      ))
  expect_output(print(ip), "        d  content")
  expect_output(print(ip), "D = 1.7")

  # -------------------------------------------------------------------------- #
  # GRM
  ip <- itempool(list(
    new("Item", model = "GRM", id = "Item 1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "Item 2", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567, 0.646),
                          D = 1.702))
            ))
  expect_output(print(ip), "id     a     b1     b2     b3    b4    content")
  expect_output(print(ip), "Model of items: ")
  expect_output(print(ip), "D = 1.702")

  # -------------------------------------------------------------------------- #
  # An item pool consist of Testlet can be printed nicely.
  t1 <- testlet(itempool(b = rnorm(5), id = paste0("t1-i", 1:5)),
                   id = "t1")
  t2 <- testlet(itempool(b = rnorm(4), id = paste0("t2-i", 1:4)),
                   id = "t2")
  ip <- c(t1, t2)
  expect_output(object = print(t2), regexp = "t2-i2")
  expect_output(object = print(ip), regexp = "t2-i2")

  # -------------------------------------------------------------------------- #
  # If argument n is presented only that number of items printed to the console
  ip <- generate_ip(n = 50)
  expect_output(object = print(ip), regexp = "# ... with 40 more items")
  expect_output(object = print(ip, n = 45), regexp = "# ... with 5 more items")
  expect_output(object = print(ip), regexp = "Total number of items = 50")

  # -------------------------------------------------------------------------- #
  # If n is between 11 and 20 all of the items are printed
  ip <- generate_ip(n = sample(11:19, 1))
  text <- capture.output(print(ip))
  expect_false(any(grepl("with (.*) more items", text)))
  expect_false(any(grepl("Total number of items ", text)))

  # If n is between 1 and 10 all of the items are printed
  ip <- generate_ip(n = sample(1:10, 1))
  text <- capture.output(print(ip))
  expect_false(any(grepl("with (.*) more items", text)))
  expect_false(any(grepl("Total number of items ", text)))

  # -------------------------------------------------------------------------- #
  # For PCM, there shouldn't be a "D = " in the printed output
  ip <- itempool(list(
    new("Item", model = "PCM", id = "PCM-1", content = "Stress",
        parameters = list(b = c(-0.837, 0.529))),
    new("Item", model = "PCM", id = "PCM-2", content = "Stress",
        parameters = list(b = c(-0.837, 0.529, 1.2)))
            ))
  expect_output(print(ip), "PCM-1")
  expect_output(print(ip), "b1")
  expect_false(any(grepl(pattern = "D = ", capture.output(ip))))

  # -------------------------------------------------------------------------- #
  # A mixture of testlets and items but all of the standalone items and the
  # testlet items has the same model, All D's are equal
  n_t1 <- sample(5:9, 1)
  n_t2 <- sample(4:7, 1)
  n_t3 <- 1
  t1 <- testlet(itempool(
    a = rlnorm(n_t1, 0, .3), b = rnorm(n_t1), c = runif(n_t1, 0, .3),
    d = runif(n_t1, .95, 1), id = paste0("t1-i", 1:n_t1)), id = "t1")
  t2 <- testlet(itempool(
    a = rlnorm(n_t2, 0, .3), b = rnorm(n_t2), c = runif(n_t2, 0, .3),
    d = runif(n_t2, .95, 1), id = paste0("t2-i", 1:n_t2)), id = "t2")
  t3 <- testlet(itempool(
    a = rlnorm(n_t3, 0, .3), b = rnorm(n_t3), c = runif(n_t3, 0, .3),
    d = runif(n_t3, .95, 1), id = paste0("t3-i", 1:n_t3)), id = "t3")
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1), id = paste0("i1-", 1:4))
  i2 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     d = runif(2, .95, 1), id = paste0("i2-", 1:2))
  ip <- c(t1, i1, t2, i2, t3)
  expect_output(print(ip), "t1-i1")
  expect_output(print(ip), "D = 1")

  # If item pools with testlets has more than 20 items, additional info will
  # be printed
  n_items <- ip$n$items
  i3 <- generate_ip(n = 15, model = "4PL", id = paste0("i3-", 1:15))
  ip <- c(ip, i3)
  expect_output(print(ip), paste0("with ", n_items + 15 - 10," more items"))
  expect_output(print(ip), paste0("Number of Testlets =  ", ip$n$testlets))
  expect_output(print(ip), paste0("Number of items within Testlets = ",
                                  ip$n$items-(ip$n$elements-ip$n$testlets)))
  expect_output(print(ip), paste0("Number of standalone items = ",
                                  ip$n$elements-ip$n$testlets))
  expect_output(print(ip), paste0("Total number of items = ", ip$n$items))

  # Change the "n = " argument for the print() function
  expect_output(print(ip, n = 13), paste0("with ", n_items + 15 - 13," more items"))
  expect_output(print(ip, n = 13), paste0("Number of Testlets =  ", ip$n$testlets))
  expect_output(print(ip, n = 13),
                paste0("Number of items within Testlets = ",
                       ip$n$items-(ip$n$elements-ip$n$testlets)))
  expect_output(print(ip, n = 13), paste0("Number of standalone items = ",
                                          ip$n$elements-ip$n$testlets))
  expect_output(print(ip, n = 13),
                paste0("Total number of items = ", ip$n$items))

  # -------------------------------------------------------------------------- #
  # A mixture of testlets and items but all of the standalone items and the
  # testlet items has the same model, D's are different
  n_t1 <- sample(3:6, 1)
  n_t2 <- sample(3:6, 1)
  n_t3 <- 1
  weird_D <- 92.81
  t1 <- testlet(itempool(
    a = rlnorm(n_t1, 0, .3), b = rnorm(n_t1), c = runif(n_t1, 0, .3),
    d = runif(n_t1, .95, 1), id = paste0("t1-i", 1:n_t1), D = 1), id = "t1")
  t2 <- testlet(itempool(
    a = rlnorm(n_t2, 0, .3), b = rnorm(n_t2), c = runif(n_t2, 0, .3),
    d = runif(n_t2, .95, 1), id = paste0("t2-i", 1:n_t2), D = weird_D),
    id = "t2")
  t3 <- testlet(itempool(
    a = rlnorm(n_t3, 0, .3), b = rnorm(n_t3), c = runif(n_t3, 0, .3),
    d = runif(n_t3, .95, 1), id = paste0("t3-i", 1:n_t3)), id = "t3")
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1), id = paste0("i1-", 1:4), D = 1)
  i2 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     d = runif(2, .95, 1), id = paste0("i2-", 1:2))
  ip <- c(t1, c(i1, i2), t2, t3)
  expect_output(print(ip), "t1-i1")
  expect_false(any(grepl(pattern = "D = ", capture.output(ip))))
  expect_output(print(ip, n = 15), paste(weird_D))

  # -------------------------------------------------------------------------- #
  # A mixture of Dichotomous items without a testlet.
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1), id = paste0("i1-", 1:4))
  i3 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     id = paste0("i3-", 1:2))
  i4 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2),
                     id = paste0("i4-", 1:2))
  i5 <- itempool(b = rnorm(3), id = paste0("i5-", 1:3))

  ip <- c(i1, i3, i4, i5)
  expect_output(print(ip), "i1-1")

  # -------------------------------------------------------------------------- #
  # A mixture of Dichotmous and Polytomous items without a testlet.
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1), id = paste0("i1-", 1:4))
  i2 <- itempool(list(
    new("Item", model = "GRM", id = "GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "GRM-2", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702))
            ))
  i3 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2), c = runif(2, 0, .3),
                     id = paste0("i3-", 1:2))
  i4 <- itempool(a = rlnorm(2, 0, .3), b = rnorm(2),
                     id = paste0("i4-", 1:2))
  i5 <- itempool(b = rnorm(3), id = paste0("i5-", 1:3))
  i6 <- itempool(list(
    new("Item", model = "GPCM", id = "GPCM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GPCM", id = "GPCM-2", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702))
            ))

  ip <- c(i1, i2, i3, i4, i5, i6)
  expect_output(print(ip), "i1-1")
  expect_output(print(ip), "GRM")

  # -------------------------------------------------------------------------- #
  # A mixture of items and testlets.
  t1 <- testlet(itempool(b = rnorm(5), id = paste0("t1-i", 1:5)),
                   id = "t1")
  t2 <- testlet(itempool(a = 1:4, b = rnorm(4),
                                id = paste0("t2-i", 1:4)), id = "t2")
  i1 <- itempool(a = rlnorm(4, 0, .3), b = rnorm(4), c = runif(4, 0, .3),
                     d = runif(4, .95, 1))
  i2 <- itempool(list(
    new("Item", model = "GRM", id = "GRM-1", content = "Stress",
        parameters = list(a = 0.888, b = c(-0.837, -0.529, -0.382, 1.649),
                          D = 1.702)),
    new("Item", model = "GRM", id = "GRM-2", content = "Depression",
        parameters = list(a = 1.081, b = c(-0.692, -0.157,  0.567), D = 1.702))
            ))

  ip <- c(t1, i1, t2, i2)

  expect_output(print(ip), "id testlet model")
  expect_output(print(ip), "t2-i1")
  expect_output(print(ip), "2PL")
  expect_output(print(ip), "content")

  # -------------------------------------------------------------------------- #
  # Mixture of testlet and items with 'misc' field. If 'misc' is common it
  # should be shown when printing.
  t1 <- testlet(itempool(b = -3:-2, id = c("t1-i1", "t1-i2"), D = 1.702),
                id = "t1", misc = list(sympson_hetter_k = 0))
  t2 <- testlet(itempool(a = c(0.2, 0.2), b = 4:5, id = c("t2-i1", "t2-i2"),
                          D = 1.702), id = "t2", misc = list(sympson_hetter_k = 1))
  i1 <- item(b = -1, D = 1.702, id = "i1", misc = list(sympson_hetter_k = .2))
  i2 <- item(b = 0, D = 1.702, id = "i2", misc = list(sympson_hetter_k = .8))
  i3 <- item(b = 1, D = 1.702, id = "i3", misc = list(sympson_hetter_k = .9))
  ip <- c(t1, t2, i1, i2, i3)
  expect_output(print(ip), "b sympson_hetter_k")

  # -------------------------------------------------------------------------- #
  # 2020-04-18: When printing an item from a testlet, instead of printing the
  # item parameters, a question mark is printed like: "testlet1 3        ?".
  i1 <- item(b = 1, id = "i1", content = "A")
  i2 <- item(a = 1, b = 1, id = "i2", content = "B")
  i3 <- item(a = 1, b = 1, c = 1, id = "i3", content = "C")
  i4 <- item(a = 1, b = c(1, 2), id = "i4", content = "D", model = "GRM")
  i5 <- item(a = 1, b = c(1, 2), id = "i5", content = "E", model = "GPCM")
  t1 <- testlet(i2, i5, id = "testlet1", content = "TZM")
  ip <- c(i1, i3, t1, i4)
  expect_false(any(grepl(pattern = "\\?", capture.output(ip))))
  # Also, when printing an element, there is an error.
  expect_output(print(ip[4]))

  # -------------------------------------------------------------------------- #
  # 2020-08-10: If 'misc' or 'content' field is common among items it should
  # be printed in the header:
  ip <- generate_ip(n = 6, content = rep("Algebra", 6),
                    misc = lapply(1:6, function(x)
                      list(Form = "A1", sh = round(runif(1), 1),
                           test = 1)))
  ip[[6]]$misc$test <- NA
  expect_output(print(ip), "Content = Algebra")
  expect_output(print(ip), "Form = A1")
  expect_output(print(ip), "sh test")


})


###############################################################################@
############################# convert_model ####################################
###############################################################################@

test_that("Test convert_model", {
  # Test for item
  item <- item(a = 1, b = 1.2, c = 0.2, D = rnorm(1, 2, .2))
  expect_null(item$d)
  item_new <- convert_model(item, "4PL")
  expect_equal(item_new$d, 1)
  expect_equal(item_new$a, item$a)
  expect_equal(item_new$b, item$b)
  expect_equal(item_new$c, item$c)
  expect_equal(item_new$D, item$D)
  # Convert to Rasch
  item_new <- convert_model(item, "Rasch")
  expect_null(item_new$c)
  expect_null(item_new$a)
  expect_null(item_new$D)
  expect_equal(item_new$b, item$b)

  # -------------------------------------------------------------------------- #
  # Convert Rasch to others
  item <- item(b = 1.2, model = "Rasch")
  expect_null(item$d)

  # -------------------------------------------------------------------------- #
  # One cannot convert GRM to unidimensional IRT models.
  item <- item(a = 1, b = c(-1, 1.2), D = rnorm(1, 2, .2))
  expect_error(convert_model(ip = item, target_model = "2PL"))

  # -------------------------------------------------------------------------- #
  # Convert GPCM to GPCM2
  item <- generate_item(model = "GPCM2")
  item_new <- convert_model(item, target_model = "GPCM")
  expect_equal(item_new$b, item$b-item$d)
  expect_equal(item_new$model, "GPCM")
  theta <- rnorm(1)
  expect_equal(prob(ip = item, theta = theta),
               prob(ip = item_new, theta = theta))


  # -------------------------------------------------------------------------- #
  # Test Errors
  item <- item(a = 1, b = 2, c = .2)
  expect_error(convert_model(12, target_model = "3PL"),
               regexp = "Invalid object type.")

  expect_error(convert_model(item, target_model = "XYZ"),
               regexp = "Invalid target_model.")

  # -------------------------------------------------------------------------- #
  # Reduce itemparameters
  expect_true(is.null(convert_model(item(a = 1.2, b = .14, c = .3),
                                   target_model = "2PL")@parameters$c))
  expect_true(is.null(convert_model(item(a = 1.2, b = .14, c = .3),
                                   target_model = "1PL")@parameters$a))

  # -------------------------------------------------------------------------- #
  # Increase itemparameters
  expect_true(convert_model(item(a = 1.2, b = .14, c = .3),
                                   target_model = "4PL")@parameters$d == 1)

  expect_true(convert_model(item(b = 1.2),
                           target_model = "2PL")@parameters$a == 1)

  # -------------------------------------------------------------------------- #
  # Test Item Pool
  n <- 10 # number of items
  ip <- itempool(data.frame(a = runif(n, .5, 1.5), b = rnorm(n),
                                c = runif(n, 0,.3)), id = paste0("Item-",1:n),
                     content = sample(c("Algebra", "Arithmetic", "Geometry"),
                                     n, replace = TRUE))
  expect_equivalent(convert_model(ip, "2PL")@item_list[[n]]@parameters$b,
                    ip@item_list[[n]]@parameters$b)

  # -------------------------------------------------------------------------- #
  # Test Testlet
  t1 <- testlet(itempool(a = c(.5, 1.5, 2), b = rnorm(3),
                                c = c(.1, .2, .3)))
  expect_equal(t1[[1]]$model, "3PL")
  expect_equal(t1[[1]]$c, .1)
  expect_null(t1[[1]]$d)
  t1_new <- convert_model(t1, "2PL")
  expect_equal(t1_new[[1]]$model, "2PL")
  expect_null(t1_new[[1]]$c)
  t1_new <- convert_model(t1, "4PL")
  expect_equal(t1_new[[1]]$model, "4PL")
  expect_equal(t1_new[[1]]$d, 1)



})


###############################################################################@
############################# generate_item ####################################
###############################################################################@
test_that("Test generate_item", {
  # Test the default values:
  item <- generate_item()
  expect_is(item, "Item")
  expect_equal(item$model, "3PL")

  # -------------------------------------------------------------------------- #
  # The length of model should be 1 and it cannot be NULL.
  expect_error(generate_item(model = c("3PL", "2PL")),
               "Invalid model argument.")
  expect_error(generate_item(model = NULL), "Invalid model argument.")
  expect_error(generate_item(model = "ABC"), "Invalid model argument.")

  # -------------------------------------------------------------------------- #
  # n_categories cannot be a vector of length more than 1
  expect_error(generate_item(n_categories = 1:2),
               "Invalid number of categories")
  # n_categories cannot be NULL or NA
  expect_error(generate_item(n_categories = NULL, model = "GPCM"),
               "Invalid number of categories")
  expect_error(generate_item(n_categories = NA, model = "GPCM"),
               "Invalid number of categories")

  # -------------------------------------------------------------------------- #
  # if n_categories is a decimal it will be converted to an integer
  item <- generate_item(model = "GRM", n_categories = 4.8)
  expect_equal(item$model, "GRM")
  expect_equal(length(item$parameters$b), 3)

  # -------------------------------------------------------------------------- #
  # Can create a GPCM2 item.
  expect_is(item <- generate_item(model = "GPCM2", n_categories = 4), "Item")
  expect_equal(item$model, "GPCM2")

  # -------------------------------------------------------------------------- #
  # When se_parameters = TRUE, the function will generate se_parameters.
  item <- generate_item(model = "3PL", se_parameters = TRUE)
  expect_false(is.null(item@se_parameters))
  expect_true(item@se_parameters$a > 0.05)
  expect_true(item@se_parameters$a < 0.75)
  expect_true(item@se_parameters$b > 0.05)
  expect_true(item@se_parameters$b < 0.75)
  expect_true(item@se_parameters$c > 0.05)
  expect_true(item@se_parameters$c < 0.75)
  expect_null(item@se_parameters$D)

  # -------------------------------------------------------------------------- #
  # When se_parameters = TRUE, the function will generate se_parameters for
  # polytomous items as well.
  n_categories <- sample(3:7, 1)
  item <- generate_item(model = "GPCM", n_categories = n_categories,
                        se_parameters = TRUE)
  expect_false(is.null(item@se_parameters))
  expect_true(item@se_parameters$a > 0.05)
  expect_true(item@se_parameters$a < 0.75)
  expect_equal(length(item@parameters$b),  n_categories-1)
  expect_equal(length(item@se_parameters$b), n_categories-1)
  expect_true(all(item@se_parameters$b > 0.05))
  expect_true(all(item@se_parameters$b < 0.75))
  expect_null(item@se_parameters$D)

})



###############################################################################@
############################# generate_ip ######################################
###############################################################################@
test_that("Test generate_ip", {
  ip <- generate_ip()
  expect_is(ip, "Itempool")
  # Designate the number of items
  ip <- generate_ip(n = 12)
  expect_is(ip, "Itempool")
  expect_equal(length(ip), 12)
  # Generate item pools for other models
  ip <- generate_ip("Rasch")
  expect_true(all(ip$model == "Rasch"))
  ip <- generate_ip("1PL")
  expect_true(all(ip$model == "1PL"))
  ip <- generate_ip("2PL")
  expect_true(all(ip$model == "2PL"))
  ip <- generate_ip("4PL")
  expect_true(all(ip$model == "4PL"))
  ip <- generate_ip("GPCM")
  expect_true(all(ip$model == "GPCM"))
  ip <- generate_ip("PCM")
  expect_true(all(ip$model == "PCM"))
  ip <- generate_ip("GRM")
  expect_true(all(ip$model == "GRM"))
  ip <- generate_ip("GPCM2")
  expect_true(all(ip$model == "GPCM2"))

  # -------------------------------------------------------------------------- #
  # Function can add se_parameters
  ip <- generate_ip(n = 2, model = "2PL",
                    se_parameters = list(list(a = .2, b = .3),
                                         list(a = .4, b = .5)))
  expect_equal(ip[[1]]@se_parameters$a, 0.2)
  expect_equal(ip[[1]]@se_parameters$b, 0.3)
  expect_equal(ip[[2]]@se_parameters$a, 0.4)
  expect_equal(ip[[2]]@se_parameters$b, 0.5)

  # -------------------------------------------------------------------------- #
  # n_categories should be a vector of length 1 or length n, otherwise an error
  # will be raised
  expect_error(generate_ip(n = 4, n_categories = c(5, 6)),
               "Invalid n_categories value.")

  # -------------------------------------------------------------------------- #
  # When n_categories is a vector, the polytomous items will have that number
  # of categories.
  ip <- generate_ip(model = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
                    n_categories = c(2, 3, 6, 7, 2))

  expect_equal(ip[[2]]$model, "GRM")
  expect_equal(ip[[3]]$model, "GPCM")
  expect_equal(ip[[4]]$model, "GRM")
  expect_equal(length(ip[[2]]$parameters$b), 2)
  expect_equal(length(ip[[3]]$parameters$b), 5)
  expect_equal(length(ip[[4]]$parameters$b), 6)

  # -------------------------------------------------------------------------- #
  # When n_categories is a vector of length 1 all of the polytomous items will
  # have the same number of categories
  ip <- generate_ip(model = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
                    n_categories = 7)
  expect_equal(length(ip[[2]]$parameters$b), 6)
  expect_equal(length(ip[[3]]$parameters$b), 6)
  expect_equal(length(ip[[4]]$parameters$b), 6)

  # -------------------------------------------------------------------------- #
  # When se_parameters = TRUE, the function will generate se_parameters.
  ip <- generate_ip(model = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
                    n_categories = c(2, 3, 6, 7, 2), se_parameters = TRUE)
  expect_true(ip[[1]]@se_parameters$a > 0.05)
  expect_true(ip[[1]]@se_parameters$c < 0.75)
  expect_equal(length(ip[[4]]@se_parameters$b), 6)
})


###############################################################################@
############################# generate_testlet #################################
###############################################################################@
test_that("Test generate_ip", {
  testlet <- generate_testlet()
  expect_is(testlet, "Testlet")
  expect_equal(testlet$model, "BTM")
  # Designate the number of items
  testlet <- generate_testlet(n = 12)
  expect_is(testlet, "Testlet")
  expect_equal(length(testlet), 12)
  # Generate item pools for other models
  testlet <- generate_testlet(item_models = "Rasch")
  testlet$item_model
  expect_true(all(testlet$item_models == "Rasch"))

  # -------------------------------------------------------------------------- #
  # When n_categories is a vector, the polytomous items will have that number
  # of categories.
  testlet <- generate_testlet(
    item_models = c("3PL", "GRM", "GPCM", "GRM", "2PL"),
    n_categories = c(2, 3, 6, 7, 2))

  expect_equal(testlet[[2]]$model, "GRM")
  expect_equal(testlet[[3]]$model, "GPCM")
  expect_equal(testlet[[4]]$model, "GRM")
  expect_equal(length(testlet[[2]]$parameters$b), 2)
  expect_equal(length(testlet[[3]]$parameters$b), 5)
  expect_equal(length(testlet[[4]]$parameters$b), 6)

})
