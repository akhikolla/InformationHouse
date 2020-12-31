# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/edmst_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# edmst_test_check = edmst.test(coords = coords, cases = floor(nydf$cases),
#                               pop = nydf$pop, w = nyw,
#                               alpha = 0.12, longlat = TRUE,
#                               nsim = 19, ubpop = 0.1, ubd = 0.2)
#
# context("check edmst.test with reference")
# test_that("edmst.test and edmst_test_ref match", {
#   for (i in seq_along(edmst_test_ref$clusters)) {
#     expect_equal(edmst_test_ref$clusters[[i]]$locids, edmst_test_check$clusters[[i]]$locids)
#     expect_equal(edmst_test_ref$clusters[[i]]$coords, edmst_test_check$clusters[[i]]$centroid)
#     expect_equal(edmst_test_ref$clusters[[i]]$r, edmst_test_check$clusters[[i]]$r)
#     expect_equal(edmst_test_ref$clusters[[i]]$max_dist, edmst_test_check$clusters[[i]]$max_dist)
#     expect_equal(edmst_test_ref$clusters[[i]]$pop, edmst_test_check$clusters[[i]]$pop)
#     expect_equal(edmst_test_ref$clusters[[i]]$cases, edmst_test_check$clusters[[i]]$cases)
#     expect_equal(edmst_test_ref$clusters[[i]]$ex, edmst_test_check$clusters[[i]]$expected)
#     expect_equal(edmst_test_ref$clusters[[i]]$smr, edmst_test_check$clusters[[i]]$smr)
#     expect_equal(edmst_test_ref$clusters[[i]]$rr, edmst_test_check$clusters[[i]]$rr)
#     expect_equal(edmst_test_ref$clusters[[i]]$loglikrat, edmst_test_check$clusters[[i]]$loglikrat)
#     expect_equal(edmst_test_ref$clusters[[i]]$test_statistic, edmst_test_check$clusters[[i]]$test_statistic)
#     expect_equal(edmst_test_ref$clusters[[i]]$pvalue, edmst_test_check$clusters[[i]]$pvalue)
#     expect_equal(edmst_test_ref$clusters[[i]]$w, edmst_test_check$clusters[[i]]$w)
#   }
# })
