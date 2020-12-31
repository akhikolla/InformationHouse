test_that("seqN_by works", {
  skip_if_not_installed("data.table")
  library(data.table)
  DM1 <- data.table(j0_k = c(7, 10, 4, 10, 7, 9, 8, 7, 10, 1, 1, 4, 7, 1, 10, 5, 8, 4, 6,
                             10, 10, 6, 6, 2, 10, 8, 5, 3, 6, 8),
                   k5 = 0.5)
  DM1[, j0_k := as.integer(j0_k)]
  setkeyv(DM1, "j0_k")
  DM2 <- copy(DM1)
  DM2[, seq_n := seq_len(.N), keyby = "j0_k"]
  DM1[, seq_n := seqN_by(j0_k)]
  expect_equal(DM2, DM1)
})
