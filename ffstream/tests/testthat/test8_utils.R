context("Test  8: utils, mostly the computation of pvalues")

test_that("convert a pvalue to make sure significant pvalues are small: ", {
        expect_equal(convertPvalueToCorrectSideR(0.95), 0.05)
        expect_equal(convertPvalueToCorrectSideR(0.01), 0.01)
        })


test_that("computation of one-sided pvalue: ", {
        expect_equal(computeOneSidedPvalueR(1, 0, 1), 0.025, tolerance=1e-5)
        })


test_that("combine two pvalues into single pvalue: ", {
        expect_equal(combineTwoOneSidedPvaluesR(0.025, 0.4), 0.05)
        })



test_that("computation of two-sided pvalue: ", {
        expect_equal(computeTwoSidedPvalueR(1, 0, 1), 0.025, tolerance=1e-5)
        expect_equal(computeTwoSidedPvalueR(0, 0, 1), 0.025, tolerance=1e-5)
        })



test_that("computation of normal cdf: ", {
        q1 <- 1.95
        #accurate to 1.5e-7
       expect_equal(computeStdNormCdf(q1), pnorm(q1), tolerance=2e-7) 
       
        q2 <- 2.00
       expect_equal(computeStdNormCdf(q2), pnorm(q2), tolerance=2e-7) 
        
        })



test_that("computation of pvalue for AFF: ", {
        expect_equal(makeTwoSidedPvalueOneSidedR(0.95), 0.1)
        expect_equal(makeTwoSidedPvalueOneSidedR(0.05), 0.1)
        
        })



