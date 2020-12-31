test_that("validate_nchar1 works", {
  expect_equal(validate_nchar1(c("", " ", "\u2212")), 3)
  expect_equal(validate_nchar1(c("", " ", "\u2212"), return_size = TRUE), 3)
  expect_equal(validate_nchar1(" "), 0)
  expect_equal(max_charsize("a"), 1)
})

test_that("is_space", {
  expect_true(is_space(" "))
  expect_false(is_space(""))
})

test_that("where_square_bracket_opens", {
  s <- function(x) unlist(strsplit(x, split = "", fixed = TRUE))
  expect_equal(where_square_bracket_opens(s("0123[567]9"), -1), -1)
  expect_equal(where_square_bracket_opens(s("0123[567]9"), 0), -1)
  expect_equal(where_square_bracket_opens(s("0123[567]9"), 8L), 4L)
  expect_equal(where_square_bracket_opens(s("01234567]9"), 8L), 0L)
})
