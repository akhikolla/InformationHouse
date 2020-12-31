# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/cepp_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# cepp_test_check = cepp.test(coords = coords,
#                             cases = floor(nydf$cases),
#                             pop = nydf$pop,
#                             nstar = 1000, alpha = 1)
#
# context("check cepp.test with reference")
# test_that("cepp.test and cepp_test_ref match", {
#   for (i in seq_along(cepp_test_ref$clusters)) {
#     expect_equal(cepp_test_ref$clusters[[i]]$locids, cepp_test_check$clusters[[i]]$locids)
#     expect_equal(cepp_test_ref$clusters[[i]]$coords, cepp_test_check$clusters[[i]]$centroid)
#     expect_equal(cepp_test_ref$clusters[[i]]$r, cepp_test_check$clusters[[i]]$r)
#     expect_equal(cepp_test_ref$clusters[[i]]$max_dist, cepp_test_check$clusters[[i]]$max_dist)
#     expect_equal(cepp_test_ref$clusters[[i]]$pop, cepp_test_check$clusters[[i]]$pop)
#     expect_equal(cepp_test_ref$clusters[[i]]$cases, cepp_test_check$clusters[[i]]$cases)
#     expect_equal(cepp_test_ref$clusters[[i]]$ex, cepp_test_check$clusters[[i]]$expected)
#     expect_equal(cepp_test_ref$clusters[[i]]$smr, cepp_test_check$clusters[[i]]$smr)
#     expect_equal(cepp_test_ref$clusters[[i]]$rr, cepp_test_check$clusters[[i]]$rr)
#     # expect_equal(cepp_test_ref$clusters[[i]]$loglikrat, cepp_test_check$clusters[[i]]$loglikrat)
#     expect_equal(cepp_test_ref$clusters[[i]]$test_statistic, cepp_test_check$clusters[[i]]$test_statistic)
#     expect_equal(cepp_test_ref$clusters[[i]]$pvalue, cepp_test_check$clusters[[i]]$pvalue)
#     expect_equal(cepp_test_ref$clusters[[i]]$w, cepp_test_check$clusters[[i]]$w)
#   }
# })
