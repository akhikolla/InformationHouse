test_that("op2M agrees with do_op2M", {
  the_operators <-
    c("!=", "==", ">=",
      "<=", ">", "<",
      "%in%",
      "%between%",
      "%(between)%",
      "%]between[%")

  expect_identical(sapply(the_operators, op2M),
                   sapply(the_operators, do_op2M))
})


test_that("single_ox_x1_x2 works", {
  expect_true(test_single_ox_x1_x2(1L, op2M("!="), -1L, -1L))
  expect_true(test_single_ox_x1_x2(1L, op2M("=="), 1L, 1L))
  expect_true(test_single_ox_x1_x2(1L, op2M(">="), 0L, 0L))
  expect_true(test_single_ox_x1_x2(1L, op2M("%between%"), -1L, 2L))
  expect_true(test_single_ox_x1_x2(1L, op2M("%(between)%"), -1L, 2L))
  expect_false(test_single_ox_x1_x2(1L, op2M("%]between[%"), -1L, 2L))

  expect_true(test_single_ox_x1_x2(1, op2M("!="), -1, -1))
  expect_true(test_single_ox_x1_x2(1, op2M("=="), 1, 1))
  expect_true(test_single_ox_x1_x2(1, op2M(">="), 0, 0))
  expect_true(test_single_ox_x1_x2(1, op2M("%between%"), -1, 2))
  expect_true(test_single_ox_x1_x2(1, op2M("%(between)%"), -1, 2))
  expect_false(test_single_ox_x1_x2(1, op2M("%]between[%"), -1, 2))
})
