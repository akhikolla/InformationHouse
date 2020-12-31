# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/uls_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# uls_test_check = uls.test(coords = coords, cases = floor(nydf$cases),
#                pop = nydf$pop, w = nyw,
#                alpha = 0.05, longlat = TRUE,
#                nsim = 19, ubpop = 0.5)
#
# context("check uls.test with reference")
# test_that("uls.test and uls_test_ref match", {
#   for (i in seq_along(uls_test_ref$clusters)) {
#     expect_equal(uls_test_ref$clusters[[i]]$locids, uls_test_check$clusters[[i]]$locids)
#     expect_equal(uls_test_ref$clusters[[i]]$coords, uls_test_check$clusters[[i]]$centroid)
#     expect_equal(uls_test_ref$clusters[[i]]$r, uls_test_check$clusters[[i]]$r)
#     expect_equal(uls_test_ref$clusters[[i]]$max_dist, uls_test_check$clusters[[i]]$max_dist)
#     expect_equal(uls_test_ref$clusters[[i]]$pop, uls_test_check$clusters[[i]]$pop)
#     expect_equal(uls_test_ref$clusters[[i]]$cases, uls_test_check$clusters[[i]]$cases)
#     expect_equal(uls_test_ref$clusters[[i]]$ex, uls_test_check$clusters[[i]]$expected)
#     expect_equal(uls_test_ref$clusters[[i]]$smr, uls_test_check$clusters[[i]]$smr)
#     expect_equal(uls_test_ref$clusters[[i]]$rr, uls_test_check$clusters[[i]]$rr)
#     expect_equal(uls_test_ref$clusters[[i]]$loglikrat, uls_test_check$clusters[[i]]$loglikrat)
#     expect_equal(uls_test_ref$clusters[[i]]$test_statistic, uls_test_check$clusters[[i]]$test_statistic)
#     expect_equal(uls_test_ref$clusters[[i]]$pvalue, uls_test_check$clusters[[i]]$pvalue)
#     expect_equal(uls_test_ref$clusters[[i]]$w, uls_test_check$clusters[[i]]$w)
#   }
# })
