
# library(testthat)

###############################################################################@
################################ testlet #######################################
###############################################################################@

# Test concatenation function "c" of "Itempool" class
test_that("testlet", {
  i1 <- item(b=rnorm(1))
  i2 <- item(b=rnorm(1))
  i3 <- item(b=rnorm(1), id = "mi-item-1")
  i4 <- item(b=rnorm(1), id = "myitem4")
  # A sequence of items
  expect_is(t1 <- testlet(i1, i2, i3, i4), 'Testlet')
  expect_equal(t1@model, 'BTM')
  expect_is(t1 <- testlet(i1, i2, i3, i4, id = "my_teslet1"), 'Testlet')
  expect_equal(t1@id, 'my_teslet1')
  # The order of items inthe item_list should follow: i1, i2, i3, i4
  expect_equal(t1@item_list$id, c("Item-1", "Item-2", "mi-item-1", "myitem4"))
  # A list of items
  expect_is(t1 <- testlet(list(i1, i2, i3, i4)), 'Testlet')
  expect_equal(t1@item_list[[2]]@id, 'Item-2')
  expect_is(t1 <- testlet(list(i1, i2, i3, i4), id = 'mt1'), 'Testlet')
  expect_equal(t1@id, 'mt1')
  # an Itempool
  ip <- itempool(list(i1, i2, i3, i4))
  expect_is(t1 <- testlet(ip), 'Testlet')
  expect_is(t1 <- testlet(ip, id = 'mt1'), 'Testlet')
  expect_equal(t1@id, 'mt1')
  # a data frame
  ip_dtf <- data.frame(a = runif(10, .5, 1.5), b = rnorm(10))
  expect_is(t1 <- testlet(ip_dtf, id = 'mt1'), 'Testlet')
  expect_equal(t1@id, 'mt1')

  # Testlet with 'misc' field.
  t1 <- testlet(itempool(b = rnorm(2), id = paste0("t1-i", 1:2),
                         misc = list(list(sympson_hetter_k = .8, form = "b3"),
                                     list(sympson_hetter_k = .9))),
                   id = "t1")
  expect_is(t1, "Testlet")

})

###############################################################################@
############################# $ method (Testlet) ###############################
###############################################################################@

test_that("$ method (Testlet)", {
  i1 <- item(a = 1.2, b = 0.2, c = .1, content = "Geo")
  i2 <- item(a = 1.38, b = -2.1, c = .2, content = "Geo")
  i3 <- item(a = 1.38, b = -1.1, content = "Geo")
  t1 <- testlet(c(i1, i2, i3), id = "testlet--1", content = "Alg")
  expect_equal(t1$id, "testlet--1")
  expect_equal(t1$content, "Alg")
  expect_null(t1$parameters)
  expect_is(t1$item_list, "list")
  expect_true(all(sapply(t1$item_list, is.Item)))
  expect_equivalent(t1$item_models, c("3PL", "3PL", "2PL"))

  # $max_score
  t1 <- testlet(generate_ip(model = c("2PL", "GRM", "3PL", "GRM"),
                            n_categories = c(2, 3, 2, 6)),
                id = "testlet--1", content = "Alg")
  expect_equal(t1$max_score, 9)

})

###############################################################################@
############################# $<- method (Testlet) #############################
###############################################################################@

test_that("$<- method (Testlet)", {
  i1 <- item(a = 1.2, b = 0.2, c = .1)
  i2 <- item(a = 1.38, b = -2.1, c = .2)
  i3 <- item(a = 1.38, b = -1.1)
  t1 <- testlet(c(i1, i2, i3), id = "t1")
  expect_equal(t1$item_list[[1]]$a, 1.2)
  expect_equal(t1$id, "t1")
  t1$id <- "t2"
  expect_equal(t1$id, "t2")
  t1$item_list <- convert_model(itempool(t1$item_list), "Rasch")
  expect_null(t1$item_list[[1]]$a)
})

