# fpath = system.file("testdata",  package = "smerc")
# fname = paste(fpath, "/fast_test_ref.rda", sep = "")
# load(fname)
#
# # test updated methods
# data(nydf)
# data(nyw)
# coords = with(nydf, cbind(longitude, latitude))
# set.seed(1)
# fast_test_check = fast.test(coords = coords, cases = floor(nydf$cases),
#                             pop = nydf$pop,
#                             alpha = 0.1, longlat = TRUE,
#                             nsim = 19, ubpop = 0.5)
#
# context("check fast.test with reference")
# test_that("fast.test and fast_test_ref match", {
#   for (i in seq_along(fast_test_ref$clusters)) {
#     expect_equal(fast_test_ref$clusters[[i]]$locids, fast_test_check$clusters[[i]]$locids)
#     expect_equal(fast_test_ref$clusters[[i]]$coords, fast_test_check$clusters[[i]]$centroid)
#     expect_equal(fast_test_ref$clusters[[i]]$r, fast_test_check$clusters[[i]]$r)
#     expect_equal(fast_test_ref$clusters[[i]]$max_dist, fast_test_check$clusters[[i]]$max_dist)
#     expect_equal(fast_test_ref$clusters[[i]]$pop, fast_test_check$clusters[[i]]$pop)
#     expect_equal(fast_test_ref$clusters[[i]]$cases, fast_test_check$clusters[[i]]$cases)
#     expect_equal(fast_test_ref$clusters[[i]]$ex, fast_test_check$clusters[[i]]$expected)
#     expect_equal(fast_test_ref$clusters[[i]]$smr, fast_test_check$clusters[[i]]$smr)
#     expect_equal(fast_test_ref$clusters[[i]]$rr, fast_test_check$clusters[[i]]$rr)
#     expect_equal(fast_test_ref$clusters[[i]]$loglikrat, fast_test_check$clusters[[i]]$loglikrat)
#     expect_equal(fast_test_ref$clusters[[i]]$test_statistic, fast_test_check$clusters[[i]]$test_statistic)
#     expect_equal(fast_test_ref$clusters[[i]]$pvalue, fast_test_check$clusters[[i]]$pvalue)
#     expect_equal(fast_test_ref$clusters[[i]]$w, fast_test_check$clusters[[i]]$w)
#   }
# })
