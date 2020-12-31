
context("Parameter input constant_muscle_mass")

cons <-  list(airDensity = 3)

data <- data("birds")

test_that("missing data throw error", {
  expect_error(.constant_muscle_mass(speed_control = "constant_speed",
                                     cons = cons))
})

test_that("missing constants throw error", {
  expect_error(.constant_muscle_mass(data = data, peed_control = "constant_speed"))
})


test_that("missing speed control throw error", {
  expect_error(.constant_muscle_mass(data = data, speed_control = "constant_speed",
                                    cons = cons))
})


