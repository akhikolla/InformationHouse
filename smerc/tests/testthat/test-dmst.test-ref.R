# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/dmst_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# dmst_test_check = dmst.test(coords = coords, cases = floor(nydf$cases),
#                             pop = nydf$pop, w = nyw,
#                             alpha = 0.12, longlat = TRUE,
#                             nsim = 19, ubpop = 0.1, ubd = 0.2)
#
# context("check dmst.test with reference")
# test_that("dmst.test and dmst_test_ref match", {
#   for (i in seq_along(dmst_test_ref$clusters)) {
#     expect_equal(dmst_test_ref$clusters[[i]]$locids, dmst_test_check$clusters[[i]]$locids)
#     expect_equal(dmst_test_ref$clusters[[i]]$coords, dmst_test_check$clusters[[i]]$centroid)
#     expect_equal(dmst_test_ref$clusters[[i]]$r, dmst_test_check$clusters[[i]]$r)
#     expect_equal(dmst_test_ref$clusters[[i]]$max_dist, dmst_test_check$clusters[[i]]$max_dist)
#     expect_equal(dmst_test_ref$clusters[[i]]$pop, dmst_test_check$clusters[[i]]$pop)
#     expect_equal(dmst_test_ref$clusters[[i]]$cases, dmst_test_check$clusters[[i]]$cases)
#     expect_equal(dmst_test_ref$clusters[[i]]$ex, dmst_test_check$clusters[[i]]$expected)
#     expect_equal(dmst_test_ref$clusters[[i]]$smr, dmst_test_check$clusters[[i]]$smr)
#     expect_equal(dmst_test_ref$clusters[[i]]$rr, dmst_test_check$clusters[[i]]$rr)
#     expect_equal(dmst_test_ref$clusters[[i]]$loglikrat, dmst_test_check$clusters[[i]]$loglikrat)
#     expect_equal(dmst_test_ref$clusters[[i]]$test_statistic, dmst_test_check$clusters[[i]]$test_statistic)
#     expect_equal(dmst_test_ref$clusters[[i]]$pvalue, dmst_test_check$clusters[[i]]$pvalue)
#     expect_equal(dmst_test_ref$clusters[[i]]$w, dmst_test_check$clusters[[i]]$w)
#   }
# })
