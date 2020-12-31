# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/mlf_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# mlf_test_check = mlf.test(coords = coords, cases = floor(nydf$cases),
#                           pop = nydf$pop, w = nyw,
#                           alpha = 0.12, longlat = TRUE,
#                           nsim = 19, ubpop = 0.1, ubd = 0.5)
#
# context("check mlf.test with reference")
# test_that("mlf.test and mlf_test_ref match", {
#   for (i in seq_along(mlf_test_ref$clusters)) {
#     expect_equal(mlf_test_ref$clusters[[i]]$locids, mlf_test_check$clusters[[i]]$locids)
#     # expect_equal(mlf_test_ref$clusters[[i]]$coords, mlf_test_check$clusters[[i]]$centroid)
#     # expect_equal(mlf_test_ref$clusters[[i]]$r, mlf_test_check$clusters[[i]]$r)
#     # expect_equal(mlf_test_ref$clusters[[i]]$max_dist, mlf_test_check$clusters[[i]]$max_dist)
#     expect_equal(mlf_test_ref$clusters[[i]]$pop, mlf_test_check$clusters[[i]]$pop)
#     expect_equal(mlf_test_ref$clusters[[i]]$cases, mlf_test_check$clusters[[i]]$cases)
#     expect_equal(mlf_test_ref$clusters[[i]]$ex, mlf_test_check$clusters[[i]]$expected)
#     expect_equal(mlf_test_ref$clusters[[i]]$smr, mlf_test_check$clusters[[i]]$smr)
#     expect_equal(mlf_test_ref$clusters[[i]]$rr, mlf_test_check$clusters[[i]]$rr)
#     expect_equal(mlf_test_ref$clusters[[i]]$loglikrat, mlf_test_check$clusters[[i]]$loglikrat)
#     # expect_equal(mlf_test_ref$clusters[[i]]$test_statistic, mlf_test_check$clusters[[i]]$test_statistic)
#     expect_equal(mlf_test_ref$clusters[[i]]$pvalue, mlf_test_check$clusters[[i]]$pvalue)
#     expect_equal(mlf_test_ref$clusters[[i]]$w, mlf_test_check$clusters[[i]]$w)
#   }
# })
