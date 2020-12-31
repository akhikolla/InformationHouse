context("motif_info")

test_that("All positions are included",{
  expect_equal(length(intersect(unlist(lapply(1:44, function(x) motif_info(x, link = FALSE))), 1:148)), 148)
})

test_that("All positions occur exactly once",{
  expect_equal(all(table(unlist(lapply(1:44, function(x) motif_info(x, link = FALSE)))) == 1), TRUE)
})
