context("Ensemble model")

u <- v <- t(NPpc[601:813])
u.list <- list(u, u[1:2,])
v.list <- list(v, v[1:2,])

foreach::registerDoSEQ()

test_that("Ensemble model works", {
  expect_is(LDS_reconstruction(NPannual, u.list, v.list, start.year = 1800, num.restarts = 2), "list")
})


test_that("Ensemble cross validation works", {
  expect_is(cvLDS(NPannual, u.list, v.list, start.year = 1800, num.restarts = 2,
                  Z = make_Z(NPannual$Qa, nRuns = 2, contiguous = FALSE)),
            "list")
})
