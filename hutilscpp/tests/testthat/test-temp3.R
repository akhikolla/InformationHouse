test_that("temp3 works", {
  expect_identical(do_count_logical_long(c(TRUE, FALSE, NA)), rep(1, 3))
  expect_true(do_in_int(1L, 1:4))
  expect_false(do_in_int(0L, 1:4))

  expect_true(test_intersect3_stdint())
  expect_true(test_intersect3_stdint(4))
  expect_true(test_intersect3_stdint(5))
})
