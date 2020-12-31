# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/mlink_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# mlink_test_check = mlink.test(coords = coords, cases = floor(nydf$cases),
#                                pop = nydf$pop, w = nyw,
#                                alpha = 0.12, longlat = TRUE,
#                                nsim = 19, ubpop = 0.1, ubd = 0.2)
#
# context("check mlink.test with reference")
# test_that("mlink.test and mlink_test_ref match", {
#   for (i in seq_along(mlink_test_ref$clusters)) {
#     expect_equal(mlink_test_ref$clusters[[i]]$locids, mlink_test_check$clusters[[i]]$locids)
#     expect_equal(mlink_test_ref$clusters[[i]]$coords, mlink_test_check$clusters[[i]]$centroid)
#     expect_equal(mlink_test_ref$clusters[[i]]$r, mlink_test_check$clusters[[i]]$r)
#     expect_equal(mlink_test_ref$clusters[[i]]$max_dist, mlink_test_check$clusters[[i]]$max_dist)
#     expect_equal(mlink_test_ref$clusters[[i]]$pop, mlink_test_check$clusters[[i]]$pop)
#     expect_equal(mlink_test_ref$clusters[[i]]$cases, mlink_test_check$clusters[[i]]$cases)
#     expect_equal(mlink_test_ref$clusters[[i]]$ex, mlink_test_check$clusters[[i]]$expected)
#     expect_equal(mlink_test_ref$clusters[[i]]$smr, mlink_test_check$clusters[[i]]$smr)
#     expect_equal(mlink_test_ref$clusters[[i]]$rr, mlink_test_check$clusters[[i]]$rr)
#     expect_equal(mlink_test_ref$clusters[[i]]$loglikrat, mlink_test_check$clusters[[i]]$loglikrat)
#     expect_equal(mlink_test_ref$clusters[[i]]$test_statistic, mlink_test_check$clusters[[i]]$test_statistic)
#     expect_equal(mlink_test_ref$clusters[[i]]$pvalue, mlink_test_check$clusters[[i]]$pvalue)
#     expect_equal(mlink_test_ref$clusters[[i]]$w, mlink_test_check$clusters[[i]]$w)
#   }
# })
