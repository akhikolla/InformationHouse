# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/dc_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# dc_test_check = dc.test(coords = coords, cases = floor(nydf$cases),
#                         pop = nydf$pop, w = nyw,
#                         alpha = 0.12, longlat = TRUE,
#                         nsim = 19, ubpop = 0.1, ubd = 0.2)
#
# context("check dc.test with reference")
# test_that("dc.test and dc_test_ref match", {
#   for (i in seq_along(dc_test_ref$clusters)) {
#     expect_equal(dc_test_ref$clusters[[i]]$locids, dc_test_check$clusters[[i]]$locids)
#     expect_equal(dc_test_ref$clusters[[i]]$coords, dc_test_check$clusters[[i]]$centroid)
#     expect_equal(dc_test_ref$clusters[[i]]$r, dc_test_check$clusters[[i]]$r)
#     expect_equal(dc_test_ref$clusters[[i]]$max_dist, dc_test_check$clusters[[i]]$max_dist)
#     expect_equal(dc_test_ref$clusters[[i]]$pop, dc_test_check$clusters[[i]]$pop)
#     expect_equal(dc_test_ref$clusters[[i]]$cases, dc_test_check$clusters[[i]]$cases)
#     expect_equal(dc_test_ref$clusters[[i]]$ex, dc_test_check$clusters[[i]]$expected)
#     expect_equal(dc_test_ref$clusters[[i]]$smr, dc_test_check$clusters[[i]]$smr)
#     expect_equal(dc_test_ref$clusters[[i]]$rr, dc_test_check$clusters[[i]]$rr)
#     expect_equal(dc_test_ref$clusters[[i]]$loglikrat, dc_test_check$clusters[[i]]$loglikrat)
#     expect_equal(dc_test_ref$clusters[[i]]$test_statistic, dc_test_check$clusters[[i]]$test_statistic)
#     expect_equal(dc_test_ref$clusters[[i]]$pvalue, dc_test_check$clusters[[i]]$pvalue)
#     expect_equal(dc_test_ref$clusters[[i]]$w, dc_test_check$clusters[[i]]$w)
#   }
# })
