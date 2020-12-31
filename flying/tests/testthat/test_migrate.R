
context("minimum energy from protein")

data <- data("birds")

test_that("out of bounds protein_met throws error", {
  expect_error(migtate(data = data, method = "cmm",
                       speed_control = "constant_speed", protein_met = 1.2
                                    ))
})
