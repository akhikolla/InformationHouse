context("stopover mass calculator")

test_that("mismatch length throw error", {
  expect_error(stopover.mass.calculator(bodyMass = c(1.1,1.2), fatMass = c(0.34),
                                        taxon = c(1, 2), duration = 24L
  ))
})


test_that("negaative duration throws error", {
  expect_error(stopover.mass.calculator(bodyMass = c(1.1,1.2), fatMass = c(0.34, 0.25),
                                        taxon = c(1, 2), duration = -24L
  ))
})


test_that("negaative duration throws error", {
  expect_error(stopover.mass.calculator(bodyMass = c(1.1,1.2), fatMass = c(0.34, 0.25),
                                        taxon = c(1, 3), duration = 24L
  ))
})
