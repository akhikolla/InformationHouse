data_colnames = c("lon", "x", "lat", "r")

test_that("arg_check_coordnames", {
  coordnames = ~ lon + lat
  # should produce coordinate names
  expect_equal(arg_check_coordnames(coordnames, data_colnames),
               c("lon", "lat"))
  # lat3 not in data_colnames
  coordnames = ~ lon + lat3
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  coordnames = ~ lon + lat + abc
  # too many columns of data
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  coordnames = c("lon", "lat")
  # should produce coordinate names
  expect_equal(arg_check_coordnames(coordnames, data_colnames),
               c("lon", "lat"))

  # lat3 not in data_colnames
  coordnames = c("lon", "lat3")
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  coordnames = c("lon", "lat", "abc")
  # too many columns of data
  expect_error(arg_check_coordnames(coordnames, data_colnames))


  coordnames = c(1, 3)
  # should produce coordinate names
  expect_equal(arg_check_coordnames(coordnames, data_colnames),
               c("lon", "lat"))

  # column 5 is not in data_colnames
  coordnames = c(1, 5)
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  # negative indices not allowed
  coordnames = c(-1, 3)
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  # column 5 not allowed
  coordnames = c(5, 3)
  expect_error(arg_check_coordnames(coordnames, data_colnames))

  coordnames = 1:3
  # too many columns of data
  expect_error(arg_check_coordnames(coordnames, data_colnames))
})
