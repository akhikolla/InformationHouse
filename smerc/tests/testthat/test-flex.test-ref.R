# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/flex_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# flex_test_check = flex.test(coords = coords, cases = floor(nydf$cases),
#                             w = nyw, k = 10,
#                             pop = nydf$pop, nsim = 19,
#                             alpha = 0.12, longlat = TRUE)
#
# context("check flex.test with reference")
# test_that("flex.test and flex_test_ref match", {
#   for (i in seq_along(flex_test_ref$clusters)) {
#     expect_equal(flex_test_ref$clusters[[i]]$locids, flex_test_check$clusters[[i]]$locids)
#     expect_equal(flex_test_ref$clusters[[i]]$coords, flex_test_check$clusters[[i]]$centroid)
#     expect_equal(flex_test_ref$clusters[[i]]$r, flex_test_check$clusters[[i]]$r)
#     expect_equal(flex_test_ref$clusters[[i]]$max_dist, flex_test_check$clusters[[i]]$max_dist)
#     expect_equal(flex_test_ref$clusters[[i]]$pop, flex_test_check$clusters[[i]]$pop)
#     expect_equal(flex_test_ref$clusters[[i]]$cases, flex_test_check$clusters[[i]]$cases)
#     expect_equal(flex_test_ref$clusters[[i]]$ex, flex_test_check$clusters[[i]]$expected)
#     expect_equal(flex_test_ref$clusters[[i]]$smr, flex_test_check$clusters[[i]]$smr)
#     expect_equal(flex_test_ref$clusters[[i]]$rr, flex_test_check$clusters[[i]]$rr)
#     expect_equal(flex_test_ref$clusters[[i]]$loglikrat, flex_test_check$clusters[[i]]$loglikrat)
#     expect_equal(flex_test_ref$clusters[[i]]$test_statistic, flex_test_check$clusters[[i]]$test_statistic)
#     expect_equal(flex_test_ref$clusters[[i]]$pvalue, flex_test_check$clusters[[i]]$pvalue)
#     expect_equal(flex_test_ref$clusters[[i]]$w, flex_test_check$clusters[[i]]$w)
#   }
# })
