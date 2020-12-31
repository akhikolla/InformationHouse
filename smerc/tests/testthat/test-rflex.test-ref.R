# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/rflex_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# rflex_test_check = rflex.test(coords = coords, cases = floor(nydf$cases),
#                              w = nyw, k = 15,
#                              pop = nydf$pop, nsim = 19,
#                              alpha = 0.1, longlat = TRUE)
#
# context("check rflex.test with reference")
# test_that("rflex.test and rflex_test_ref match", {
#   for (i in seq_along(rflex_test_ref$clusters)) {
#     expect_equal(rflex_test_ref$clusters[[i]]$locids, rflex_test_check$clusters[[i]]$locids)
#     expect_equal(rflex_test_ref$clusters[[i]]$coords, rflex_test_check$clusters[[i]]$centroid)
#     expect_equal(rflex_test_ref$clusters[[i]]$r, rflex_test_check$clusters[[i]]$r)
#     expect_equal(rflex_test_ref$clusters[[i]]$max_dist, rflex_test_check$clusters[[i]]$max_dist)
#     expect_equal(rflex_test_ref$clusters[[i]]$pop, rflex_test_check$clusters[[i]]$pop)
#     expect_equal(rflex_test_ref$clusters[[i]]$cases, rflex_test_check$clusters[[i]]$cases)
#     expect_equal(rflex_test_ref$clusters[[i]]$ex, rflex_test_check$clusters[[i]]$expected)
#     expect_equal(rflex_test_ref$clusters[[i]]$smr, rflex_test_check$clusters[[i]]$smr)
#     expect_equal(rflex_test_ref$clusters[[i]]$rr, rflex_test_check$clusters[[i]]$rr)
#     expect_equal(rflex_test_ref$clusters[[i]]$loglikrat, rflex_test_check$clusters[[i]]$loglikrat)
#     expect_equal(rflex_test_ref$clusters[[i]]$test_statistic, rflex_test_check$clusters[[i]]$test_statistic)
#     expect_equal(rflex_test_ref$clusters[[i]]$pvalue, rflex_test_check$clusters[[i]]$pvalue)
#     expect_equal(rflex_test_ref$clusters[[i]]$w, rflex_test_check$clusters[[i]]$w)
#   }
# })
