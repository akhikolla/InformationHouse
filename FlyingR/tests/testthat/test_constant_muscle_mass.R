
context("Parameter input constant_muscle_mass")

cons <-  list(airDensity = 3)

data <- data("birds")

test_that("missing data throw error", {
  expect_error(.constant_muscle_mass(speed_control = 1,
                                     constants = cons))
})

test_that("missing constants throw error", {
  expect_error(.constant_muscle_mass(data = data, speed_control = 1))
})


test_that("missing speed control throw error", {
  expect_error(.constant_muscle_mass(data = data, speed_control =1,
                                    constants = cons))
})


