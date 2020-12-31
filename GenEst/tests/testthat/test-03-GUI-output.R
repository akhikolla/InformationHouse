context("Check the Output-Producing Functions")

test_that("Figure size functions return numeric values", {

  expect_error(setFigW("_NOTRIGHT_"))

  data(mock)
  modSet <- pkmSet(data = mock$SE, formula_p = p ~ 1, formula_k = k ~ 1)
  expect_is(setFigW(modSet), "numeric")
  expect_equal(setFigW(modSet), 800)
  expect_is(setFigH(modSet), "numeric")
  expect_equal(setFigH(modSet), 800)

  modSet <- cpmSet(data = mock$CP, formula_l = l ~ 1, formula_s = s ~ 1,
              left = "LastPresentDecimalDays", 
              right = "FirstAbsentDecimalDays"
            )
  expect_is(setFigW(modSet), "numeric")
  expect_equal(setFigW(modSet), 800)
  expect_is(setFigH(modSet), "numeric")
  expect_equal(setFigH(modSet, "CP"), 700)
})