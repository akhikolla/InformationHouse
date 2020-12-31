# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/elliptic_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# elliptic_test_check = elliptic.test(coords = coords,
#                                      cases = floor(nydf$cases),
#                                      pop = nydf$pop, ubpop = 0.1,
#                                      nsim = 19,
#                                      alpha = 0.12)
#
# context("check elliptic.test with reference")
# test_that("elliptic.test and elliptic_test_ref match", {
#   for (i in seq_along(elliptic_test_ref$clusters)) {
#     expect_equal(elliptic_test_ref$clusters[[i]]$locids, elliptic_test_check$clusters[[i]]$locids)
#     expect_equal(elliptic_test_ref$clusters[[i]]$coords, elliptic_test_check$clusters[[i]]$centroid)
#     expect_equal(elliptic_test_ref$clusters[[i]]$r, elliptic_test_check$clusters[[i]]$r)
#     expect_equal(elliptic_test_ref$clusters[[i]]$max_dist, elliptic_test_check$clusters[[i]]$max_dist)
#     expect_equal(elliptic_test_ref$clusters[[i]]$pop, elliptic_test_check$clusters[[i]]$pop)
#     expect_equal(elliptic_test_ref$clusters[[i]]$cases, elliptic_test_check$clusters[[i]]$cases)
#     expect_equal(elliptic_test_ref$clusters[[i]]$ex, elliptic_test_check$clusters[[i]]$expected)
#     expect_equal(elliptic_test_ref$clusters[[i]]$smr, elliptic_test_check$clusters[[i]]$smr)
#     expect_equal(elliptic_test_ref$clusters[[i]]$rr, elliptic_test_check$clusters[[i]]$rr)
#     expect_equal(elliptic_test_ref$clusters[[i]]$loglikrat, elliptic_test_check$clusters[[i]]$loglikrat)
#     expect_equal(elliptic_test_ref$clusters[[i]]$test_statistic, elliptic_test_check$clusters[[i]]$test_statistic)
#     expect_equal(elliptic_test_ref$clusters[[i]]$pvalue, elliptic_test_check$clusters[[i]]$pvalue)
#     expect_equal(elliptic_test_ref$clusters[[i]]$w, elliptic_test_check$clusters[[i]]$w)
#     expect_equal(elliptic_test_ref$clusters[[i]]$semiminor_axis, elliptic_test_check$clusters[[i]]$semiminor_axis)
#     expect_equal(elliptic_test_ref$clusters[[i]]$semimajor_axis, elliptic_test_check$clusters[[i]]$semimajor_axis)
#     expect_equal(elliptic_test_ref$clusters[[i]]$angle, elliptic_test_check$clusters[[i]]$angle)
#     expect_equal(elliptic_test_ref$clusters[[i]]$shape, elliptic_test_check$clusters[[i]]$shape)
#   }
# })
